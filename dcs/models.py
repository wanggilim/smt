#! /usr/bin/env python
from peewee import Model,SqliteDatabase,TextField,IntegerField,DoubleField
import json
from pathlib import Path
from functools import reduce,partial
from bs4 import BeautifulSoup
from astropy.coordinates import SkyCoord
import astropy.units as u

DB_TYPE_MAP = {'int':IntegerField,'float':DoubleField,'str':TextField}

def generate_field(key, units, options):
    '''Generate database Field objects for each key'''
    # get field type, default to string
    ftype = units.get(key,'str')
    
    # convert ftypes to sql fields
    field = DB_TYPE_MAP[ftype]
    opt = options.get('_ALL_').copy()
    opt.update(options.get(key,{}))  # override options for field
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

def insert_rows(db, rows, cls):
    """Insert rows to database"""
    with db.atomic():
        cls.insert_many(rows).execute()
def replace_rows(db, rows, cls):
    """Upssert rows to database"""
    with db.atomic():
        cls.replace_many(rows).execute()

def ModelFactory(name, config, db, register=False):
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
        
    return cls


CONVERTER_FUNCS = {'AOR':AOR_to_rows,
                   'MIS':MIS_to_rows}

if __name__ == '__main__':
    from configparser import ConfigParser
    c = ConfigParser()
    c.read('models.cfg')


    #AOR_to_rows('../test/07_0193.aor',c['AOR'])

    MIS_to_rows('../test/201909_HA_FABIO.misxml',c['MIS'])
