#! /usr/bin/env python
from peewee import Model,SqliteDatabase,TextField,IntegerField,DoubleField,ForeignKeyField,BooleanField,TimeField,DateTimeField
import json
from pathlib import Path
from functools import reduce,partial
from bs4 import BeautifulSoup
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
#import POS as dcsPOS
try:
    from .FAOR import FAOR as dcsFAOR
except (ModuleNotFoundError,ImportError):
    from FAOR import FAOR as dcsFAOR
try:
    from . import MIS as dcsMIS
except (ModuleNotFoundError,ImportError):
    import MIS as dcsMIS
#from .POS import POS as dcsPOS
from fifi import SCT as dcsSCT
from collections import deque
from pandas import read_html, concat, DataFrame
import numpy as np
import re
import datetime

DB_TYPE_MAP = {'int':IntegerField,'float':DoubleField,'str':TextField,'bool':BooleanField,
               'foreign':ForeignKeyField,'time':TimeField,'datetime':DateTimeField}
STR_TYPE_MAP = {'int':int,'float':float,'str':str,'bool':bool}
HTMLPARSER = 'lxml'

POS_SPLIT_RE = re.compile(r'\#-*\sTarget:\s.*\s-*\#\n')
POS_CONFIG_RE = re.compile(r'\#\d\d\_([^,\|]*\|){15}\n')
GUIDE_RE = re.compile(r'\#\s(?P<Cat>.*)\,\s(?P<Cam>FPI-TO|WFI|FFI)\,\sRadius=(?P<Radius>[\d\.]*)\,\s(?P<MagData>.*)\n(?P<Data>.*)')


def _as_pandas(rows):
    if isinstance(rows,dict):
        columns = rows.keys()
    else:
        columns = rows[0].keys()
    return DataFrame.from_records(rows,columns=columns)

def _as_table(rows):
    if isinstance(rows,dict):
        columns = rows.keys()
    else:
        columns = rows[0].keys()
    #print(rows[0])
    #print(columns)
    try:
        table = Table(data=rows,names=columns)
    except ValueError:
        rows = [{k:row.get(k,None) for k in columns} for row in rows]
        table = Table(data=rows,names=columns)
        
    return table

def _as_json(rows):
    return json.dumps(rows,default=str)

def convert_dtypes(row, unitdict):
    unitdict = {k:STR_TYPE_MAP.get(v,str) for k,v in unitdict.items()}
    row = {k:unitdict[k](v) if k in unitdict and v is not None else v for k,v in row.items()}
    return row

def generate_field(key, units, options):
    '''Generate database Field objects for each key'''
    # get field type, default to string
    ftype = units.get(key,'str')
    
    # convert ftypes to sql fields
    field = DB_TYPE_MAP[ftype]
    opt = options.get('_ALL_').copy()
    opt.update(options.get(key,{}))  # override options for field

    # connect foreign key to model
    if ftype == 'foreign':
        model = globals()[opt['model']]
        opt['model'] = model
        
    return field(**opt)


def get_keydict(root,keys,exclude=None,as_tag=False):
    if exclude:
        if isinstance(exclude,str):
            exclude = (exclude,)
        keys = list(filter(lambda x:x not in exclude, keys))
    
    try:
        if as_tag:
            return {key:getattr(root,key.lower()) for key in keys}
        else:
            return {key:getattr(root,key.lower()).text for key in keys}
    except AttributeError:
        vals = (root.find(key.lower()) for key in keys)
        if not as_tag:
            vals = (val.text if val else None for val in vals)
        return {key:val for key,val in zip(keys,vals)}

def get_attrdict(tags,data_attr,ext='_'):
    """Expand out _start and _end and get attrs"""
    try:
        attrdict = {ext.join((k,attr)):tag[attr] for k,tag in tags.items() for attr in data_attr[k]}
    except TypeError:
        # likely a non-science leg
        attrdict = {ext.join((k,attr)):None for k,tag in tags.items() for attr in data_attr[k]}
    return attrdict

def proc_waypoint(tag,keys,data_units):
    """Get waypoint as dict"""
    waypt = dict()
    for k in keys:
        v = tag.get(k.lower())
        v = v if v != 'N/A' else None
        v = data_units.get(k,str)(v) if v is not None else None
        waypt[k] = v
        
    return waypt

def get_waypoints(leg,keys,data_units):
    """Get waypoint data for leg"""
    try:
        waypts = leg.waypoints.find_all('waypoint')
    except AttributeError:
        return None
    data_units = {k:STR_TYPE_MAP.get(v,str) for k,v in data_units.items()}
    waypt_func = partial(proc_waypoint,keys=keys,data_units=data_units)
    '''
    waypts = map
    waypts = map(lambda tag: {k:tag.get(k.lower()) for k in keys},waypts)
    waypts = {k:v if v != 'N/A' for k,v in wa
    print(list(waypts)[0])
    '''
    waypts = tuple(map(waypt_func,waypts))
    return waypts
        
def make_position(request):
    try:
        ra =  float(request.target.position.lon.string)
        dec = float(request.target.position.lat.string)
    except AttributeError:
        # target likely ephemeris
        return None
    try:
        equ = request.target.position.coordSystem.equinoxDesc.string
    except AttributeError:
        equ = 'J2000'
    coord = SkyCoord(ra=ra,dec=dec,equinox=equ,unit=(u.deg,u.deg))
    return coord

def combine_AOR_data(name,position,rkeys,dkeys,dithers,maps,blkdict,comments,meta):
    row = meta.copy()
    row.update(rkeys)
    row.update(dkeys)
    row.update(dithers)
    row.update(maps)
    row['target'] = name
    try:
        row['RA'],row['DEC'] = position.to_string('hmsdms',sep=':').split()
    except AttributeError:
        row['RA'] = None
        row['DEC'] = None
    row['ObsBlkID'] = blkdict.get(row['aorID'],None)
    if comments:
        row['ObsBlkComment'] = comments.get(row['ObsBlkID'])

    if row['InstrumentName'] == 'FORCAST':
        # THIS IS NEW AND POSSIBLY DANGEROUS
        if row.get('TotalTime') and row.get('Repeat'):
            row['TotalTime'] = float(row['TotalTime']) * float(row['Repeat'])

    return row

def combine_GUIDE_data(target,aorid,dkeys,akeys,mkeys,blkdict,meta):
    row = dkeys
    row.update(akeys)
    row.update(mkeys)
    row.update(meta)
    row['aorID'] = aorid
    row['target'] = target
    row['planID'] = '_'.join(aorid.split('_')[0:-1])
    row['pkey'] = ': '.join((aorid,row['Name']))
    row['ObsBlkID'] = blkdict.get(row['aorID'],None)
    return row

def combine_MIS_data(legnum,dkeys,attrs,waypts,meta):
    row = meta.copy()
    row['Leg'] = legnum
    row.update(dkeys)
    row.update(attrs)
    row['WAYPTS'] = json.dumps(waypts) if waypts else None

    #row['ObsBlk'] = row['ObsBlkID']
    row['planID'] = row['ObsPlanID']
    row['fkey'] = 'Leg%s_%s' % (legnum,row['FlightPlan'])
    try:
        row['Comment'] = row['Comment'].replace('<br>','\n')
    except AttributeError:
        row['Comment'] = None
        
    return row

def combine_FAOR_data(config,run,meta):
    row = meta.copy()
    row.update(config)
    row.update(run)
    row['ckey'] = '%s_%s' % (row['fkey'],row['AORID'])

    return row

def combine_AORSEARCH_data(row):
    aorid = row['AORID'].replace('ViewPlan','').strip()
    row['AORID'] = aorid
    row['planID'] = '_'.join(aorid.split('_')[0:-1])

    flightplans = row['FlightPlanIDs']
    if flightplans:
        flightplans = ', '.join(flightplans.replace('[OB]','').split())

    '''
    if isinstance(flightplans,str):
        flightplans = ', '.join(flightplans.replace('[OB]','').split())
    elif np.isnan(flightplans):
        flightplans = ''
    else:
        flighplans = ''
    '''
    row['FlightPlanIDs'] = flightplans

    return row

def replace_keys(row, keymap):
    for k,v in keymap.items():
        newval = row.pop(k)
        if newval:
            row[v] = newval
    return row

def get_fifi_dithers(request):
    if request.instrument.data.instrumentname.text != 'FIFI-LS':
        return {'deltaRaV':None,'deltaDecW':None}

    try:
        offsets = request.instrument.ditheroffsets.find_all('ditheroffset')
    except:
        return {'deltaRaV':None,'deltaDecW':None}
    
    if not offsets:
        return {'deltaRaV':None,'deltaDecW':None}

    offsets = ((float(offset.deltarav.text),float(offset.deltadecw.text)) for offset in offsets)
    dra,ddec = zip(*offsets)
    dra,ddec = json.dumps(dra), json.dumps(ddec)
    return {'deltaRaV':dra,'deltaDecW':ddec}

def get_great_map(request):
    if request.instrument.data.instrumentname.text != 'GREAT':
        return {'MapPos':None}

    try:
        maps = request.instrument.mappos.find_all('indexedlocation')
    except:
        return {'MapPos':None}
    
    if not maps:
        return {'MapPos':None}

    maps = [[m['label'],int(m.index.text),float(m.ra.text),float(m.dec.text)] for m in maps]
    maps = json.dumps(maps)
    
    #offsets = ((float(offset.deltarav.text),float(offset.deltadecw.text)) for offset in offsets)
    #dra,ddec = zip(*offsets)
    #dra,ddec = json.dumps(dra), json.dumps(ddec)
    return {'MapPos':maps}
    

def AOR_to_rows(filename, aorcfg, convert_dtype=False):
    """Converts AOR xml files to rows for DB ingestion"""
    with open(filename,'r') as f:
        try:
            aor = BeautifulSoup(f,HTMLPARSER).body.aors.list
        except AttributeError:
            return None

    meta = get_keydict(aor,json.loads(aorcfg['meta_keys']))
    try:
        pi = aor.investigator.attrs
        inst = pi['institution']
        if pi['honorific'] in ('','NONE',None):
            pi['honorific'] = ''
        pi = ' '.join((pi['honorific'],pi['firstname'],pi['lastname'])).strip()
        pi = ', '.join((pi,inst))
    except:
        pi = ''
    meta['PI'] = pi

    # copy ProposalID to planID
    meta['planID'] = meta['ProposalID']

    # save filename
    meta['FILENAME'] = str(Path(filename).resolve())

    # save timestamp
    stats = Path(filename).stat()
    ts = stats.st_mtime if stats.st_mtime > stats.st_ctime else stats.st_ctime
    meta['TIMESTAMP'] = ts

    # Get requests, and pull xml keys
    requests = aor.vector.find_all('request')
    req_func = partial(get_keydict,keys=json.loads(aorcfg['request_keys']))
    rkeys = map(req_func,requests)

    # Get target name and position info from request
    names = (r.target.find('name').text for r in requests)
    positions = (make_position(r) for r in requests)

    # Get instrument config info from request.data
    keys = json.loads(aorcfg['data_keys'])
    fkeys = json.loads(aorcfg['FORCAST_keys'])
    fikeys = json.loads(aorcfg['FIFI_keys'])
    grkeys = json.loads(aorcfg['GREAT_keys'])
    keys += fkeys.keys()
    keys += fikeys
    #keys += grkeys

    data_func = partial(get_keydict,keys=keys)
    dkeys = map(data_func,(r.data for r in requests))

    # replace forcast specific keys
    rep_func = partial(replace_keys,keymap=fkeys)
    dkeys = map(rep_func,dkeys)

    # get fifi dithers
    dithers = map(get_fifi_dithers,requests)

    # get great maps
    maps = map(get_great_map,requests)

    # Get obsblk info from request.obsplanobsblockinfolist
    try:
        obsplans = aor.obsplanobsblockinfolist.find_all('obsblockinfo')
    except AttributeError:
        obsplans = None

    try:
        comments = {obsplan.obsblockid.text:obsplan.comment.text for obsplan in obsplans}
    except (AttributeError,TypeError):
        comments = None


    if obsplans is not None:
        # mapping of obsblk:[aorid,aorid]...
        obsblkdict = {obs.obsblockid.text:map(lambda x:x.text,obs.find_all('aorid')) for obs in obsplans}

        # reverse dictionary, to have aorids mapped to obsblocks
        blkdict = {aorid:block for block,aorlist in obsblkdict.items() for aorid in aorlist}
    else:
        # likely a calibrator
        aorlist = (r.data.aorid.text for r in requests)
        blkdict = {aorid:None for aorid in aorlist}

    # combine all xml data into row
    row_func = partial(combine_AOR_data,blkdict=blkdict,comments=comments,meta=meta)
    rows = map(row_func,names,positions,rkeys,dkeys,dithers,maps)

    if convert_dtype:
        units = json.loads(aorcfg['data_units'])
        con_func = partial(convert_dtypes,unitdict=units)
        rows = map(con_func,rows)
        
    return list(rows)


def MIS_to_rows(filename, miscfg):
    """Converts MIS xml files to rows for DB ingestion"""
    with open(filename,'r') as f:
        mis = BeautifulSoup(f,HTMLPARSER).body.flightplan

    if mis is None or Path(filename).suffix == '.mis':
        # fallback to legacy format
        tab = Table.read(filename,format='mis-tab')
        rows = dcsMIS.MIS_table_to_DB(tab,miscfg)
        return rows
        
    #wps = mis.find_all('waypoint')
    #[w.extract() for w in wps]

    # Get flight info
    meta = {"FlightPlan":mis['id'],"Series":'_'.join(mis['id'].split('_')[0:-1])}

    # get additional flight header info
    mkeys = json.loads(miscfg['meta_keys'])
    exkeys = ('FlightPlan','Leg')  # ignore these
    for mkey in mkeys:
        if mkey in exkeys:
            continue
        if 'Sun' in mkey:
            try:
                meta[mkey] = mis.sunsetrise.select_one(mkey).text
            except AttributeError:
                meta[mkey] = None
        else:
            meta[mkey] = mis.get(mkey.lower())

    # format datetime keys
    dunits = json.loads(miscfg['data_units'])
    for k,v in dunits.items():
        if v == 'datetime':
            if meta[k] is None:
                continue
            try:
                dt = datetime.datetime.strptime(meta[k], '%Y-%m-%dT%H:%M:%SZ')
            except ValueError:
                dt = datetime.datetime.strptime(meta[k], '%Y-%b-%d %H:%M:%S %Z')
            meta[k] = dt.strftime('%Y-%m-%d %H:%M:%S')

    # save filename
    meta['FILENAME'] = str(Path(filename).resolve())

    # set name (e.g. GAVIN)
    meta['FlightName'] = mis['id'].split('_')[-1]

    # save timestamp
    stats = Path(filename).stat()
    ts = stats.st_mtime if stats.st_mtime > stats.st_ctime else stats.st_ctime
    meta['TIMESTAMP'] = ts

    # Get legs, and pull xml keys
    legs = mis.legs.find_all('leg')
    legnums = (leg['id'] for leg in legs)

    data_func = partial(get_keydict,keys=json.loads(miscfg['data_keys']))
    dkeys = map(data_func,legs)

    # Get keys with attributes
    data_attr = json.loads(miscfg['data_attr'])
    attr_func = partial(get_keydict,keys=data_attr.keys(),as_tag=True)
    attrs = map(attr_func,legs)

    # Now get attributes (start,end)
    attr_func = partial(get_attrdict,data_attr=data_attr)
    attrs = map(attr_func,attrs)

    # Get waypoints
    waypt_func = partial(get_waypoints,
                         keys=json.loads(miscfg['waypts_keys']),
                         data_units=json.loads(miscfg['waypts_data_units']))
    waypts = list(map(waypt_func,legs))
    
    # combine all xml data into row
    row_func = partial(combine_MIS_data,meta=meta)
    rows = map(row_func,legnums,dkeys,attrs,waypts)

    return list(rows)

def FAOR_to_rows(filename, faorcfg,faor=None):
    """Converts FAOR file to rows for DB ingestion"""
    if faor is None:
        faor = dcsFAOR.read(filename)
        
    # Get meta data
    meta = {key:faor.preamble[key] for key in json.loads(faorcfg['meta_keys'])}
    meta['FlightName'] = meta['FlightPlan'].split('_')[-1]
    meta['FILENAME'] = str(Path(filename).resolve())
    meta['planID'] = '_'.join(faor.config[0]['AORID'].split('_')[0:-1])
    meta['fkey'] = 'Leg%s_%s' % (meta['Leg'],meta['FlightPlan'])

    # save timestamp
    stats = Path(filename).stat()
    ts = stats.st_mtime if stats.st_mtime > stats.st_ctime else stats.st_ctime
    meta['TIMESTAMP'] = ts

    # Get config data from each config block
    config = map(lambda c: {k:c.get(k,None) for k in json.loads(faorcfg['config_keys'])}, faor.config)

    # Get run data from each run block
    run = map(lambda r: {k:r.get(k,None) for k in json.loads(faorcfg['run_keys'])}, faor.run)

    config = list(config)

    comments = faor.get_comments_as_dicts()

    # get dithscale and comment data
    for fc,c,m in zip(faor.config,config,comments):
        c.update(m)

    # combine all xml data into row
    row_func = partial(combine_FAOR_data,meta=meta)
    rows = map(row_func,config,run)

    return list(rows)


def SCT_to_rows(filename, sctcfg, sct=None):
    """Converts SCT file to rows for DB ingestion"""
    if sct is None:
        sct = dcsSCT.read_sctfile(filename)

    # save timestamp
    stats = Path(filename).stat()
    ts = stats.st_mtime if stats.st_mtime > stats.st_ctime else stats.st_ctime
    sct['TIMESTAMP'] = ts
    return [sct]
        

def AORSEARCH_to_frame(filename, aorsearchcfg):
    """Convert AORSEARCH result to pandas frame"""
    with open(filename,'r') as f:
        soup = BeautifulSoup(f.read(),HTMLPARSER)
    ths = soup.find_all('th')
    # this table should be unique
    try:
        htable = list(filter(lambda x: 'NAIF_ID' in x.text,ths))[0]
    except IndexError:
        # no aorids found
        return None

    htable = htable.parent.parent
    table = read_html(str(htable))[0]
    table = table[json.loads(aorsearchcfg['data_keys'])]
    table['FILENAME'] = str(filename)
    return table
    

def AORSEARCH_to_rows(filenames, aorsearchcfg):
    """Converts AORSEARCH result pages for DB ingestion"""
    if isinstance(filenames,str):
        filenames = [filenames]

    frame_func = partial(AORSEARCH_to_frame, aorsearchcfg=aorsearchcfg)
    rows = map(frame_func, filenames)
    rows = filter(lambda x: x is not None, rows)
    rows = concat(rows)
    rows.fillna('', inplace=True)
    rows = rows.to_dict('records')

    #row_func = partial(combine_AORSEARCH_data)
    rows = map(combine_AORSEARCH_data,rows)
    return list(rows)

def POS_to_rows(filename, poscfg):

    data_keys = json.loads(poscfg['data_keys'])
    pos_keys = json.loads(poscfg['pos_keys'])
    
    with open(filename,'r') as f:
        text = f.read()

    rows = deque()
    target_blocks = POS_SPLIT_RE.split(text)
    for target in target_blocks[1:]:
        # get config block
        config = POS_CONFIG_RE.search(target).group()
        config = config[1:].strip().split('|')[:-1]
        config = {k:v for k,v in zip(data_keys,config)}
        config['FILENAME'] = filename
        config['planID'] = '_'.join(config['AORID'].split('_')[:-1])

        # get all nods
        nods = re.findall('%s.*J2000.*\n'%config['Target'],target)
        for idx,nod in enumerate(nods):
            nod = nod.split()
            if nod[-1][0] == '#':
                # no PM
                nod = {k:v for k,v in zip(pos_keys[:-2],nod)}
                nod['PMRA'] = ''
                nod['PMDEC'] = ''
            else:
                nod = {k:v for k,v in zip(pos_keys,nod)}
            nod['pkey'] = '%s: %s' % (config['AORID'],nod['POSName'])
            if 'Nod' in nod['POSName']:
                nod['isNod'] = True
            else:
                nod['isNod'] = False
            nod.update(config)
            
            nods[idx] = nod
        rows.extend(nods)

    return list(rows)

def GUIDE_to_rows(filename, guidecfg):
    """Converts AOR xml files to rows for DB ingestion"""
    with open(filename,'r') as f:
        try:
            aor = BeautifulSoup(f,HTMLPARSER).body.aors.list
        except AttributeError:
            return None

    # save filename
    meta = {'FILENAME':str(Path(filename).resolve())}
    # save timestamp
    stats = Path(filename).stat()
    ts = stats.st_mtime if stats.st_mtime > stats.st_ctime else stats.st_ctime
    meta['TIMESTAMP'] = ts

    requests = aor.vector.find_all('request')
    gstars = [star for r in requests for star in r.find_all('guidestar')]
    if not gstars:
        return None
    data_keys = json.loads(guidecfg['data_keys'])
    gdata_keys = ['GuideStar%s'%key for key in data_keys]
    
    data_func = partial(get_keydict,keys=gdata_keys)
    dkeys = map(data_func,gstars)
    dkeys = (dict(zip(data_keys,d.values())) for d in dkeys)

    # Get guidestar attributes
    data_attr = json.loads(guidecfg['attr_keys'])
    akeys = map(lambda tag:{a:tag[a.lower()] if tag and a.lower() in tag.attrs else None for a in data_attr},gstars)

    # Get target name from request
    targets = (g.parent.parent.find('name').text for g in gstars)
    aorids = (g.parent.parent.parent.instrument.data.aorid.text for g in gstars)
    
    # Get obsblk info from request.obsplanobsblockinfolist
    try:
        obsplans = aor.obsplanobsblockinfolist.find_all('obsblockinfo')
    except AttributeError:
        obsplans = None

    if obsplans is not None:
        # mapping of obsblk:[aorid,aorid]...
        obsblkdict = {obs.obsblockid.text:map(lambda x:x.text,obs.find_all('aorid')) for obs in obsplans}

        # reverse dictionary, to have aorids mapped to obsblocks
        blkdict = {aorid:block for block,aorlist in obsblkdict.items() for aorid in aorlist}
    else:
        # likely a calibrator
        aorlist = (r.data.aorid.text for r in requests)
        blkdict = {aorid:None for aorid in aorlist}

    # get mags
    mks = json.loads(guidecfg['mag_keys'])
    mtags = [{m:g.find('guidestar%s'%m.lower()) for m in mks} for g in gstars]
    mkeys = ({k:'%s_%s'%(v['type'],v.text) if ((v and ('type' in v.attrs)) and (np.float(v.text) < 900)) else None for k,v in m.items()} \
             for m in mtags)
    '''
    mkeys = deque()
    for m in mtags:
        try:
            mk = {k:'%s_%s'%(v['type'],v.text) if (v and ('type' in v.attrs)) and (np.float(v.text) > 900) else None for k,v in m.items()}
        except (KeyError,AttributeError):
            print(m)
            mk = {k:None for k in m.keys()}
        mkeys.append(mk)
    #mkeys = ({k:'%s_%s'%(v['type'],v.text) if np.float(v.text) != 999 else None for k,v in m.items()} for m in mtags)
    '''

    row_func = partial(combine_GUIDE_data,blkdict=blkdict,meta=meta)
    rows = map(row_func,targets,aorids,dkeys,akeys,mkeys)

    return list(rows)
    

def insert_rows(db, rows, cls):
    """Insert rows to database"""
    with db.atomic():
        cls.insert_many(rows).execute()
def replace_rows(db, rows, cls):
    """Upssert rows to database"""
    with db.atomic():
        cls.replace_many(rows).execute()


def ModelFactory(name, config, db, register=True):
    '''Generate peewee Model on the fly'''
    
    # first make meta class
    Meta = type('Meta', (object,), {'database':db,
                                    'table_name':config[name]['table_name']})

    # then generate all fields from config keys
    keys = config[name]['keys'].split('+')
    keys = map(json.loads, keys)
    keys = reduce(lambda x,y:x+y, keys)

    units = json.loads(config[name]['data_units'])
    options = json.loads(config[name]['data_options'])
    
    field_func = partial(generate_field,
                         units=units,
                         options=options)
    fields = {key:field_func(key) for key in keys}
    fields['Meta'] = Meta

    # finally, generate class
    cls = type(name, (Model,), fields)
    if register:
        # register globally
        globals()[name] = cls

    # add converter function
    cls.to_rows = CONVERTER_FUNCS.get(name,NotImplementedError())

    # add insert functions
    cls.insert_rows = partial(insert_rows,cls=cls)
    cls.replace_rows = partial(replace_rows,cls=cls)

    cls.MODEL_NAME = name

    # add rendering functions
    cls.as_pandas = _as_pandas
    #cls.as_json = lambda x: json.dumps(x)
    cls.as_json = _as_json
    cls.as_table = _as_table
    
    return cls


CONVERTER_FUNCS = {'AOR':AOR_to_rows,
                   'MIS':MIS_to_rows,
                   'FAOR':FAOR_to_rows,
                   'SCT':SCT_to_rows,
                   'POS':POS_to_rows,
                   'GUIDE':GUIDE_to_rows,
                   'AORSEARCH':AORSEARCH_to_rows}

    

if __name__ == '__main__':
    from configparser import ConfigParser
    c = ConfigParser()
    c.read('DBmodels.cfg')


    sctfile = '../fifi/202002_FI/LAWRENCE/VS44_FIFI-LS_VS44_08_0071_1.sct'
    rows = SCT_to_rows(sctfile,c['SCT'])
    print(rows)
    exit()

    aorfile = '/home/msgordo1/Downloads/07_0140.aor'
    rows = AOR_to_rows(aorfile,c['AOR'])
    print(rows[0])
    exit()
    
    misfile = '/home/msgordo1/Downloads/202001_HA_JONAS_WX12.misxml'
    rows = MIS_to_rows(misfile,c['MIS'])
    print(rows)

    exit()

    aorfile = '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/c8e3b8507f078fef10b677a8813bb058'
    aors = AOR_to_rows(aorfile,c['AOR'])
    for aor in aors:
        if aor['aorID'] == '07_0049_60':
            print(aor)

    exit()

    misfile = '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/2045115f3d5f18458abc0eb97fcbf123'

    mis = MIS_to_rows(misfile,c['MIS'])
    exit()
    
    aorfile = '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/e4289dd22b7aee455660b29a4b5e431e'
    #'../test/07_0130.aor'
    guide = GUIDE_to_rows(aorfile,c['GUIDE'])
    print(guide)
    
    exit()
    
    #posfile = '../test/201910_FO_GIMLI/201910_FO_GIMLI_V838_Mon.pos'
    posfile = '../test/07_0225.pos'
    #pos = dcsPOS.read_POS_file('../test/201910_FO_GIMLI/201910_FO_GIMLI_V838_Mon.pos',guide=True)
    pos = POS_to_rows(posfile,c['POS'])
    #print(pos)

    guide = GUIDE_to_rows(posfile,c['GUIDE'])
    print(guide)
    exit()


    #faor = dcsFAOR.read('../test/Leg13__90_0062_Alpha_Cet.faor')

    #AOR_to_rows('../test/07_0193.aor',c['AOR'])

    #MIS_to_rows('../test/201909_HA_FABIO.misxml',c['MIS'])
    #rows = AORSEARCH_to_rows(['/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/2a06a89e01eb342562a6c9b578c771bb',
    #                          '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/192432211ebc8069df5c3f84247b7d3b',
    #                          '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/7cb0d37a6d387c20148bb1c48f76b7e8',
    #                          '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/4a012c6e452ab87cbe8117ee118c603e'],
    #                         c['AORSEARCH'])
    #rows = AORSEARCH_to_rows(['../test/AORSEARCH1.html','../test/AORSEARCH2.html','../test/AORSEARCH3.html','../test/AORSEARCH4.html'],c['AORSEARCH'])
    #rows = AORSEARCH_to_rows(['/home/gordon/.astropy/cache/DCS/astropy/download/py3/9eebcae65a2a06f1c1d810d87b776a49',
    #                          '/home/gordon/.astropy/cache/DCS/astropy/download/py3/02e4ff2f922f1021029f5154cdfdf465',
    #                          '/home/gordon/.astropy/cache/DCS/astropy/download/py3/74fb34e0c12c46acd946027e694d1472',
    #                          '/home/gordon/.astropy/cache/DCS/astropy/download/py3/a72797468df001ca71dea3a3598b70c5'],
    #                         c['AORSEARCH'])
    #print(rows)
