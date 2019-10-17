#! /usr/bin/env python
from peewee import Model,SqliteDatabase,TextField,IntegerField,DoubleField,ForeignKeyField
import json
from pathlib import Path
from functools import reduce,partial
from bs4 import BeautifulSoup
from astropy.coordinates import SkyCoord
import astropy.units as u
from FAOR import FAOR as dcsFAOR
#from .FAOR import FAOR as dcsFAOR
from collections import deque
from pandas import read_html, concat
import numpy as np

DB_TYPE_MAP = {'int':IntegerField,'float':DoubleField,'str':TextField,'foreign':ForeignKeyField}

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
    attrdict = {ext.join((k,attr)):tag[attr] for k,tag in tags.items() for attr in data_attr[k]}
    # this works, I promise
    return attrdict
        
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

def combine_AOR_data(name,position,rkeys,dkeys,blkdict,meta):
    row = meta.copy()
    row.update(rkeys)
    row.update(dkeys)
    row['target'] = name
    row['RA'],row['DEC'] = position.to_string('hmsdms',sep=':').split()
    row['ObsBlk'] = blkdict.get(row['aorID'],None)
    return row

def combine_MIS_data(legnum,dkeys,attrs,meta):
    row = meta.copy()
    row['Leg'] = legnum
    row.update(dkeys)
    row.update(attrs)

    row['ObsBlk'] = row['ObsBlkID']
    row['planID'] = row['ObsPlanID']
    row['fkey'] = 'Leg%s_%s' % (legnum,row['FlightPlan'])
    
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
    

def AOR_to_rows(filename, aorcfg):
    """Converts AOR xml files to rows for DB ingestion"""
    with open(filename,'r') as f:
        aor = BeautifulSoup(f,'lxml').body.aors.list

    meta = get_keydict(aor,json.loads(aorcfg['meta_keys']))
    try:
        pi = aor.investigator.attrs
        inst = pi['institution']
        pi = ' '.join((pi['honorific'],pi['firstname'],pi['lastname']))
        pi = ', '.join((pi,inst))
    except:
        pi = ''
    meta['PI'] = pi

    # copy ProposalID to planID
    meta['planID'] = meta['ProposalID']

    # save filename
    meta['FILENAME'] = str(Path(filename).resolve())

    # Get requests, and pull xml keys
    requests = aor.vector.find_all('request')
    req_func = partial(get_keydict,keys=json.loads(aorcfg['request_keys']))
    rkeys = map(req_func,requests)

    # Get target name and position info from request
    names = (r.target.find('name').text for r in requests)
    positions = (make_position(r) for r in requests)

    # Get instrument config info from request.data
    data_func = partial(get_keydict,keys=json.loads(aorcfg['data_keys']))
    dkeys = map(data_func,(r.data for r in requests))

    # Get obsblk info from request.obsplanobsblockinfolist
    obsplans = aor.obsplanobsblockinfolist.find_all('obsblockinfo')

    # mapping of obsblk:[aorid,aorid]...
    obsblkdict = {obs.obsblockid.text:map(lambda x:x.text,obs.find_all('aorid')) for obs in obsplans}

    # reverse dictionary, to have aorids mapped to obsblocks
    blkdict = {aorid:block for block,aorlist in obsblkdict.items() for aorid in aorlist}

    # combine all xml data into row
    row_func = partial(combine_AOR_data,blkdict=blkdict,meta=meta)
    rows = map(row_func,names,positions,rkeys,dkeys)

    return list(rows)


def MIS_to_rows(filename, miscfg):
    """Converts MIS xml files to rows for DB ingestion"""
    with open(filename,'r') as f:
        mis = BeautifulSoup(f,'lxml').body.flightplan

    # Get flight info
    meta = {"FlightPlan":mis['id'],"Series":'_'.join(mis['id'].split('_')[0:-1])}

    # save filename
    meta['FILENAME'] = str(Path(filename).resolve())

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

    # combine all xml data into row
    row_func = partial(combine_MIS_data,meta=meta)
    rows = map(row_func,legnums,dkeys,attrs)

    return list(rows)

def FAOR_to_rows(filename, faorcfg):
    """Converts FAOR file to rows for DB ingestion"""
    faor = dcsFAOR.read(filename)

    # Get meta data
    meta = {key:faor.preamble[key] for key in json.loads(faorcfg['meta_keys'])}
    meta['FILENAME'] = str(Path(filename).resolve())
    meta['planID'] = '_'.join(faor.config[0]['AORID'].split('_')[0:-1])
    meta['fkey'] = 'Leg%s_%s' % (meta['Leg'],meta['Flight'])

    # Get config data from each config block
    config = map(lambda c: {k:c.get(k,None) for k in json.loads(faorcfg['config_keys'])}, faor.config)

    # Get run data from each run block
    run = map(lambda r: {k:r.get(k,None) for k in json.loads(faorcfg['run_keys'])}, faor.run)

    # combine all xml data into row
    row_func = partial(combine_FAOR_data,meta=meta)
    rows = map(row_func,config,run)

    return list(rows)

def AORSEARCH_to_frame(filename, aorsearchcfg):
    """Convert AORSEARCH result to pandas frame"""
    with open(filename,'r') as f:
        soup = BeautifulSoup(f.read(),'lxml')
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

    print(name)
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
        
    return cls


CONVERTER_FUNCS = {'AOR':AOR_to_rows,
                   'MIS':MIS_to_rows,
                   'FAOR':FAOR_to_rows,
                   'AORSEARCH':AORSEARCH_to_rows}

if __name__ == '__main__':
    from configparser import ConfigParser
    c = ConfigParser()
    c.read('models.cfg')


    #faor = dcsFAOR.read('../test/Leg13__90_0062_Alpha_Cet.faor')

    #AOR_to_rows('../test/07_0193.aor',c['AOR'])

    #MIS_to_rows('../test/201909_HA_FABIO.misxml',c['MIS'])
    #rows = AORSEARCH_to_rows(['/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/2a06a89e01eb342562a6c9b578c771bb',
    #                          '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/192432211ebc8069df5c3f84247b7d3b',
    #                          '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/7cb0d37a6d387c20148bb1c48f76b7e8',
    #                          '/home/msgordo1/.astropy/cache/DCS/astropy/download/py3/4a012c6e452ab87cbe8117ee118c603e'],
    #                         c['AORSEARCH'])
    rows = AORSEARCH_to_rows(['../test/AORSEARCH1.html','../test/AORSEARCH2.html','../test/AORSEARCH3.html','../test/AORSEARCH4.html'],
                             c['AORSEARCH'])
    #rows = AORSEARCH_to_rows(['/home/gordon/.astropy/cache/DCS/astropy/download/py3/9eebcae65a2a06f1c1d810d87b776a49',
    #                          '/home/gordon/.astropy/cache/DCS/astropy/download/py3/02e4ff2f922f1021029f5154cdfdf465',
    #                          '/home/gordon/.astropy/cache/DCS/astropy/download/py3/74fb34e0c12c46acd946027e694d1472',
    #                          '/home/gordon/.astropy/cache/DCS/astropy/download/py3/a72797468df001ca71dea3a3598b70c5'],
    #                         c['AORSEARCH'])
    print(rows)
