#! /usr/bin/env python
import argparse
from peewee import Model,SqliteDatabase,TextField,IntegerField,DoubleField
from configparser import ConfigParser
import json
from pathlib import Path
from functools import reduce,partial

DB_TYPE_MAP = {'int':IntegerField,'float':DoubleField,'str':TextField}
config = ConfigParser()

def generate_field(key, units, options):
    '''Generate database Field objects for each key'''
    # get field type, default to string
    ftype = units.get(key,'str')
    
    # convert ftypes to sql fields
    field = DB_TYPE_MAP[ftype]
    opt = options.get('_ALL_').copy()
    opt.update(options.get(key,{}))  # override options for field
    return field(**opt)

def make_table_name(model_class):
    model_name = model_class.__name__
    return model_name.lower() + '_tbl'
    #return config[model_name]['table_name']


class BaseModel(Model):
    def __new__(cls, name, bases, attrs):
        # Dynamically create class variables from AOR keys in config file
        print(cls)
        model_name = model_class.__name__
        # load all keys as json lists.  mapreduce to combine lists
        keys = config[model_name]['keys'].split('+')
        keys = map(json.loads, keys)
        keys = reduce(lambda x,y:x+y, keys)
        locals().update({key:_generate_field(key) for key in keys})

        # invoke database and table at runtime
        #database = SqliteDatabase(config['DEFAULT']['db_file'])
        #table_name = config['AOR']['table_name']
        #new_attrs = {'database':database,'table_name':table_name}

        #if 'Meta' in attrs:
        #    attrs['Meta'].update(new_attrs)
        #else:
        #    attrs['Meta'] = new_attrs
        
        return super(Model,cls).__new__(cls, name, bases, attrs)

    
    class Meta:
        table_function = make_table_name
        

class AORTEST(Model):
    """AOR database model for sqlite3"""
    def __new__(cls, name, bases, attrs):
        # Dynamically create class variables from AOR keys in config file
        aor_meta_keys = config['AOR']['aor_meta_keys']
        print('WHAT')
        locals().update({key:_generate_field(key) for key in ROW_KEYS})

        # invoke database and table at runtime
        database = SqliteDatabase(config['DEFAULT']['db_file'])
        #table_name = config['AOR']['table_name']
        #new_attrs = {'database':database,'table_name':table_name}

        if 'Meta' in attrs:
            attrs['Meta'].update(new_attrs)
        else:
            attrs['Meta'] = new_attrs
        
        return super(Model,cls).__new__(cls, name, bases, attrs)
        #database = SqliteDatabase(config['DEFAULT']['db_file'])
        #table_name = config['AOR']['table_name']
        #print(config)

        # override Meta class at runtime
        #model._meta.set_database(database)
        #model._meta.set_table_name(table_name)
        #setattr(model._meta, 'database', database)
        #setattr(model._meta, 'table_name', table_name)
    
        #return model


    #class Model:
    #    database = DATABASE
    '''    
    class Meta:
        #database = SqliteDatabase(config['DEFAULT']['db_file'])
        #table_name = config['AOR']['table_name']
        #database = SqliteDatabase(None)
    '''

def ModelFactory(name, config, db):
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
    #globals()[name] = cls # register globally
    return cls
    
    
def main():
    parser = argparse.ArgumentParser(description='Convert DCS xml to database model')
    parser.add_argument('-cfg',type=str,default='models.cfg',help='Config file with model keywords (default=models.cfg).')

    args = parser.parse_args()
    config.read(args.cfg)
    db = SqliteDatabase(config['DEFAULT']['db_file'])
    m = ModelFactory('AOR', config, db)

    # load all keys as json lists.  mapreduce to combine lists
    #keys = config['AOR']['keys'].split('+')
    #keys = map(json.loads, keys)
    #keys = reduce(lambda x,y:x+y, keys)

    #AOR = type("AOR", (Model,)

    exit()

    db = SqliteDatabase(config['DEFAULT']['db_file'])

    # bind to model
    db.bind([AOR])
    
    #db_filename = SqliteDatabase(config['DEFAULT']['db_file'])
    #db = DATABASE
    #db.initialize(SqliteDatabase(config['DEFAULT']['db_file']))
    
    with db:
        db.create_tables([AOR],safe=True)

    

    
if __name__ == '__main__':
    # TESTING
    main()
#else:
#    cfgfile = Path(__file__).parent/'models.cfg'
#    config.read(cfgfile)
#    db = SqliteDatabase(config['DEFAULT']['db_file'])
#    # bind to model
#    db.bind([AOR])
