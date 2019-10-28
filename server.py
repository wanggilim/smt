#! /usr/bin/env python
from configparser import ConfigParser
import argparse
import eel
import sqlite3
from bs4 import BeautifulSoup
from functools import partial
from astropy.coordinates import SkyCoord
import astropy.units as u
from collections import OrderedDict,defaultdict
from pathlib import Path
from peewee import Model,SqliteDatabase,TextField,IntegerField,DoubleField
import json
from dcs import DBmodels, DCS

AOR_META_KEYS = ('ProposalID','ProposalTitle','Category','PI')
REQUEST_KEYS = ('title',)
#TARGET_KEYS = ('name','position')
#INVESTIGATOR_ATTRS = ('FirstName','LastName','Honorific','Institution')
DATA_KEYS = ('aorID','InstrumentName',
             'InstrumentSpectralElement1','InstrumentSpectralElement2',
             'order','TotalTime','Repeat',
             'ObsPlanConfig','ObsPlanMode',
             'ScanType','ScanCoordSys','ScanRate','ScanTime',
             'ScanAmplitudeEL','ScanAmplitudeXEL',
             'ChopType','NodType',
             'ChopThrow','ChopAngle','ChopAngleCoordinate',
             'NodThrow','NodAngle','NodAngleCoordinate','NodTime',
             'DitherPattern','DitherCoord','DitherScale',
             'target','RA','DEC')
DATA_UNITS = {'Repeat':int,'TotalTime':float,'Repeat':int,
              'ScanRate':float,'ScanTime':float,
              'ScanAmplitudeEL':float,'ScanAmplitudeXEL':float,
              'ChopThrow':float,'ChopAngle':float,
              'NodThrow':float,'NodAngle':float,'NodTime':float,
              'DitherScale':float,'order':int}
DATA_OPTIONS = {'aorID':{'primary_key':True,'null':False,'unique':True},
                '_ALL_':{'null':True}}


ROW_KEYS = DATA_KEYS + AOR_META_KEYS + REQUEST_KEYS

DB_TYPE_DICT = {"int":IntegerField,"float":DoubleField,"str":TextField}

DB_FILE = 'smt.db'
AOR_TABLE_NAME = 'aor_tbl'
'''
class MetaAOR(type,Model):
    """Dynamically create class variables from AOR keys"""
    def __new__(mcs, clsname, bases, dictionary):
        print(clsname)
        for key in ROW_KEYS:
            # get field type, default to string
            ftype = DATA_UNITS.get(key,str)

            # convert ftypes to sql fields
            field = DB_TYPE_DICT[ftype]
            options = DATA_OPTIONS.get(key,{})

            dictionary[key] = field(**options)
            #setattr(AOR,key,field(**options))
        return type.__new__(mcs, clsname, bases, dictionary)

class AOR(metaclass=MetaAOR):
    class Meta:
        database = SqliteDatabase(DB_FILE)
        table_name = AOR_TABLE_NAME
'''

def _generate_field(key):
    '''Generate database Field objects for each key'''
    # get field type, default to string
    ftype = DATA_UNITS.get(key,str)
    
    # convert ftypes to sql fields
    field = DB_TYPE_DICT[ftype]
    options = DATA_OPTIONS.get('_ALL_').copy()
    options.update(DATA_OPTIONS.get(key,{}))  # override options for field
    return field(**options)

'''
class AORTEST(Model):
    # Dynamically create class variables from AOR keys
    locals().update({key:_generate_field(key) for key in ROW_KEYS})
    
    class Meta:
        database = SqliteDatabase(DB_FILE)
        table_name = AOR_TABLE_NAME
'''


def jsonify(func):
    def wrapped(*args,**kwargs):
        return json.dumps(func(*args,**kwargs))
    return wrapped


@eel.expose
def get_mission(stuff):
    print(stuff)
    
def initialize_database(mcfg, models=('AOR','MIS','FAOR','AORSEARCH')):
    """Initialize database, create models dynamically, and register models globally"""
    db_file = mcfg['DEFAULT']['db_file']
    db = SqliteDatabase(db_file)

    # generate models
    #   do not register FAOR in dcs namespace
    #register = (False if model == 'FAOR' else True for model in models)
    #mods = [dcs.models.ModelFactory(name, mcfg, db, reg) for name,reg in zip(models,register)]
    mods = [DBmodels.ModelFactory(name, mcfg, db) for name in models]

    # bind database to all models
    #db.bind([AOR])

    # create database
    with db:
        db.create_tables(mods,safe=True)

    # register models globally
    for name,model in zip(models,mods):
        globals()[name] = model
    return db

def create_AOR_tbl(db):
    with db:
        db.create_tables([AOR],safe=True)

def create_aor_table(con,name='aor_table',okeys=ROW_KEYS,dtypes=DATA_UNITS):
    dunits = [DB_TYPE_DICT[DATA_UNITS[key]] if key in DATA_UNITS else DB_TYPE_DICT[str] \
              for key in okeys]
    cols = ', '.join(['%s %s' % (key,dunit) for key,dunit in zip(okeys,dunits)])
    stmt = 'CREATE TABLE %s(%s);' % (name,cols)

    print(stmt)

    try:
        with con:
            con.execute(stmt)
    except sqlite3.OperationalError:
        #table exists?
        pass
    return con

"""
def insert_rows(con,rows,name='aor_table',okeys=ROW_KEYS):
    '''
    cols = ', '.join(okeys)
    qs = ', '.join(['?']*len(okeys))
    stmt = 'INSERT INTO %s(%s) VALUES (%s);' % (name,cols,qs)

    print(stmt)
    print(rows[0])

    with con:
        con.executemany(stmt, rows)

    '''
    with con:
        for row in rows:
            cols,vals = zip(*row.items())
            cols = ', '.join(cols)
            qs = ', '.join(['?']*len(vals))
            stmt = 'INSERT INTO %s(%s) VALUES (%s);' % (name,cols,qs)
            print(stmt)
            print(vals)
            con.execute(stmt,vals)
"""
def insert_rows2(db,rows):
    with db.atomic():
        AOR.insert_many(rows).execute()
def replace_rows2(db,rows):
    with db.atomic():
        AOR.replace_many(rows).execute()


def get_keydict(root,keys,exclude=None):
    if exclude:
        if isinstance(exclude,str):
            exclude = (exclude,)
        keys = list(filter(lambda x:x not in exclude, keys))
    
    try:
        return {key:getattr(root,key.lower()).text for key in keys}
    except AttributeError:
        vals = (root.find(key.lower()) for key in keys)
        vals = (val.text if val else None for val in vals)
        return {key:val for key,val in zip(keys,vals)}
        
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

def combine(name,position,rkeys,dkeys,meta,okeys=ROW_KEYS):
    row = meta.copy()
    row.update(rkeys)
    row.update(dkeys)
    row['target'] = name
    row['RA'],row['DEC'] = position.to_string('hmsdms',sep=':').split()


    #row = [(k,row[k]) for k in ROW_KEYS]
    #row = OrderedDict(row)
    #row = [row[k] for k in okeys]
    return row

def apply_types(row):
    for k,v in row.items():
        if v is None:
            continue
        row[k] = DATA_UNITS[k](v) if k in DATA_UNITS else str(v)
    return row

def make_ordered(row,okeys=ROW_KEYS):
    rows = (row[k] for k in okeys)
    return tuple(rows)

def validate(row):
    return {k:v for k,v in row.items() if v is not None}
        
def aor_to_rows(filename, aorcfg):
    with open(filename,'r') as f:
        aor = BeautifulSoup(f,'lxml').body.aors.list

    #meta = {key:aor.find(key.lower()).text for key in AOR_META_KEYS}
    #meta = {key:getattr(aor,key.lower()).text for key in AOR_META_KEYS}
    meta = get_keydict(aor,json.loads(aorcfg['meta_keys']))
    try:
        pi = aor.investigator.attrs
        inst = pi['institution']
        pi = ' '.join((pi['honorific'],pi['firstname'],pi['lastname']))
        pi = ', '.join((pi,inst))
    except:
        pi = ''
    meta['PI'] = pi

    meta['planID'] = meta['ProposalID']

    requests = aor.vector.find_all('request')
    req_func = partial(get_keydict,keys=json.loads(aorcfg['request_keys']))
    rkeys = map(req_func,requests)

    names = (r.target.find('name').text for r in requests)
    positions = (make_position(r) for r in requests)

    data_func = partial(get_keydict,keys=json.loads(aorcfg['data_keys']))
    dkeys = map(data_func,(r.data for r in requests))

    row_func = partial(combine,meta=meta)

    rows = map(row_func,names,positions,rkeys,dkeys)
    #rows = map(apply_types,rows)
    #rows = map(make_ordered,rows)
    #rows = map(validate,rows)
    return list(rows)

@eel.expose
def getAOR(aorID):
    if len(aorID.split('_')) == 2:
        # this is a planID
        return getAORs_by_planID(aorID)
    try:
        aor = AOR.get(AOR.aorID==aorID)
    except:
        # check DCS
        planID = '_'.join(aorID.split('_')[:-1])
        cfile = dcs.getObsPlan(planID)
        if cfile is None:
            return None
        return getAOR(aorID)
        
    return json.dumps(aor.__data__)

@eel.expose
def getAORs_by_planID(planID):
    aors = AOR.select().where(AOR.planID==planID)
    aors = (aor.__data__ for aor in aors)
    aors = list(aors)
    if not aors:
        # check DCS
        cfile = dcs.getObsPlan(planID)
        if cfile is None:
            return None
        return getAORs_by_planID(planID)
    
    return json.dumps(list(aors))

@eel.expose
def getAORs_by_ObsBlk(ObsBlk):
    aors = AOR.select().where(AOR.ObsBlk==ObsBlk).order_by(AOR.order,AOR.aorID)
    aors = (aor.__data__ for aor in aors)
    aors = list(aors)
    if not aors:
        # check DCS
        planID = '_'.join(ObsBlk.split('_')[1:-1])
        cfile = dcs.getObsPlan(planID)
        if cfile is None:
            return None
        return getAORs_by_ObsBlk(ObsBlk)
    
    return json.dumps(list(aors))

@eel.expose
def getFlightSeries(series):
    #flights = MIS.select().where(MIS.Series==series)
    #flights = {flight.FlightPlan for flight in flights}
    return dcs.getFlightSeries(series)
    


@eel.expose
def getFlightPlan(FlightPlan):
    legs = MIS.select().where(MIS.FlightPlan==FlightPlan).order_by(MIS.Leg)
    legs = (leg.__data__ for leg in legs)
    if not legs:
        cfile = dcs.getFlightPlan(FlightPlan)
        if cfile is None:
            return None
        
    return json.dumps(list(legs))

@eel.expose
def getFAORs_by_aorID(aorID):
    faors = FAOR.select().where(FAOR.AORID==aorID)
    faors = (faor.__data__ for faor in faors)
    return json.dumps(list(faors))

@eel.expose
def getFAORs_by_target(target):
    faors = FAOR.select().where(FAOR.TARGET==target)
    faors = (faor.__data__ for faor in faors)
    return json.dumps(list(faors))


def register_models(dcs):
    """Register active DB models in DCS instance globally"""
    active_models = dcs.get_active_models()
    if active_models:
        for name,model in active_models.items():
            globals()[name] = model

    # finally, register dcs
    globals()['dcs'] = dcs

    return active_models
    

def main():
    parser = argparse.ArgumentParser(description='Run server and startup SOFIA Mission Toolbox')
    parser.add_argument('-scfg',type=str,default='server.cfg',help='Server config file (default=server.cfg)')
    parser.add_argument('-mcfg',type=str,default='dcs/DBmodels.cfg',help='Model config file (default=dcs/DBmodels.cfg)')
    parser.add_argument('-r','--refresh-cache',
                        dest='refresh_cache',action='store_true',
                        help='Force update from DCS')
    args = parser.parse_args()

    # Read in options
    '''
    config = ConfigParser()
    config.read(args.scfg)
    eelcfg = config['EEL']  # eel config section
    
    eel.init(eelcfg['webdir'])
    eel.start(eelcfg['start'],
              host=eelcfg['host'],
              port=eelcfg['port'],
              mode=eelcfg['mode'],
              jinja=eelcfg['jinja'])
    '''


    mcfg = ConfigParser()
    mcfg.read(args.mcfg)

    dcs = DCS.DCS(modelcfg=mcfg,refresh_cache=args.refresh_cache)
    register_models(dcs)

    #dcs._force_db_sync()
    
    #rows = POS.to_rows('test/07_0225.pos',mcfg['POS'])
    #print('Digest POS')
    #POS.replace_rows(dcs.db,rows)

    #rows = GUIDE.to_rows('test/07_0225.pos',mcfg['GUIDE'])
    #print('Digest GUIDE')
    #GUIDE.replace_rows(dcs.db,rows)
    res = dcs.getPOS('07_0193',guide=True)
    #print(res)
    
    exit()
    
    #dcs._force_db_sync()
    #exit()
    #db = dcs.db
    #register_models(dcs)

    res = dcs.getAORs(['OB_07_0149_02','06_0182'], raw=True)
    #res2 = dcs.getAORs(['OB_06_0021_01'],as_json=True)
    print(res)
    #print()
    #print(res2)
    #print(str(res)==str(res2))
    exit()

    '''
    #db = initialize_database(mcfg)
    
    #exit()
    #mconfig = ConfigParser()
    #mconfig.read(args.mcfg)

    #db_file = mconfig['DEFAULT']['db_file']
    #db_file = Path(args.mcfg).parent/db_file

    #db = init_db(db_file)

    #config = ConfigParser()
    #config.read(args.cfg)
    '''

    '''
    # TESTING BLOCK1
    aorcfg = mcfg['AOR']
    rows = AOR.to_rows('test/07_0130.aor',aorcfg)
    print('Digest AOR')
    AOR.replace_rows(db,rows)

    miscfg = mcfg['MIS']
    rows = MIS.to_rows('test/201909_HA_FABIO.misxml',miscfg)
    print('Digest MIS')
    MIS.replace_rows(db,rows)

    faorcfg = mcfg['FAOR']
    rows = FAOR.to_rows('test/Leg13__90_0062_Alpha_Cet.faor',faorcfg)
    print('Digest FAOR')
    FAOR.replace_rows(db,rows)

    aorsearchcfg = mcfg['AORSEARCH']
    rows = AORSEARCH.to_rows(['test/AORSEARCH1.html','test/AORSEARCH2.html',
                              'test/AORSEARCH3.html','test/AORSEARCH4.html'],
                             aorsearchcfg)
    print('Digest AORSEARCH')
    AORSEARCH.replace_rows(db,rows)
    #'''

    # TESTING BLOCK2
    res = getAOR('07_0149_1')
    print(res)
    res = getAORs_by_planID('07_0225')
    print(res)
    res = getFlightPlan('201909_HA_FABIO')
    print(res)
    res = getAORs_by_ObsBlk('OB_07_0225_08')
    print(res)
    res = getFAORs_by_aorID('90_0062_1')
    print(res)
    res = getFAORs_by_target('Alpha Cet')
    print(res)
    res = getFlightSeries('201909_HA')
    print(res)
    #'''
    
    

if __name__ == '__main__':
    main()
    '''
    db = init_db()
    
    #create_aor_table(con)
    create_AOR_tbl(db)

    rows = aor_to_rows('test/07_0130.aor')
    replace_rows(db,rows)
    '''
