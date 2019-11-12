#! /usr/bin/env python
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table,join,Column
from astropy.io import registry
import numpy as np
import re
from collections import deque
from pandas import DataFrame

KEY_RE = re.compile('Following.*\n\#(.*)')
#LINES_RE = re.compile('\#-*\sTarget:\s(?P<Target>.*)\s-*\#\n(?P<Params>.*\s*)(?P=Target)(?P<Coord>\s.*)\n')
LINES_RE = re.compile('\#-*\sTarget:\s(?P<Target>.*)\s-*\#\n((?P<Param>#.*\n)*)(?P<Coord>.*\n)')
GUIDE_RE = re.compile('\#\s(?P<Cat>.*)\,\s(?P<Cam>FPI-TO|WFI|FFI)\,\sRadius=(?P<Radius>[\d\.]*)\,\s(?P<MagData>.*)\n(?P<Data>.*)')


def proc_line(line):
    '''Process line from POS file'''

    target,params0,params1,coord = line
    
    if not params0:
        return None


    params0 = params0.strip().split('\n')
    params1 = params1.strip().split('\n')
    params = set()

    # in case we get copies
    if len(params0) > 1:
        for p in params0:
            params.add(p)
    else:
        params.add(params0[0])
        
    if len(params1) > 1:
        for p in params1:
            params.add(p)
    
    else:
        params.add(params1[0])

    #params = [p.replace('||4_POINT||','|4_POINT|') for p in params]
    #params = [p.replace('||||','||') for p in params]
    
    coord = coord.strip().split()
    pm = [float(coord[-2]),float(coord[-1])]*u.mas/u.yr
    coord = SkyCoord(ra=coord[1],dec=coord[2],equinox=coord[3],unit=(u.hourangle,u.deg))

    return params,[coord]*len(params),[pm]*len(params)

def proc_guide_line(line):
    '''Process guide line from POS file'''
    catalog, camera, radius, magdata, data = line

    radius = u.Quantity(float(radius),u.arcmin)
    mags = magdata.split(',')
    mags = (mag.split('=') for mag in mags)
    #mags = (mag for mag in mags if not 'NOTSET' in mag)
    mags = {k.strip():float(v) for k,v in mags}

    data = {k:v for k,v in zip(('Star','RA','DEC','Equinox','PM_RA','PM_DEC'),
                               data.split())}
    coord = SkyCoord(ra=data['RA'],dec=data['DEC'],
                     equinox=data['Equinox'],unit=(u.hourangle,u.deg))
    data['RA'],data['DEC'] = coord.to_string('hmsdms',sep=':',precision=2).split()
    data['PM_RA'] =  u.Quantity(float(data['PM_RA']), u.mas/u.yr)
    data['PM_DEC'] = u.Quantity(float(data['PM_DEC']),u.mas/u.yr)

    data['Catalog'] = catalog
    data['Camera'] = camera
    data['Radius'] = radius
    data.update(**mags)

    return data

def make_table(params,coords,keys,pms):
    '''Generate table from params'''
    params = [p.split('|')[:-1] for p in params]

    tab = Table(rows=params,names=keys)
    tab['AORID'] = [x.replace('#','') for x in tab['AORID']]

    # change 90_0086_6 to 90_0086_06
    aorids = deque()
    for idx,aorid in enumerate(tab['AORID']):
        suffix = aorid.split('_')[-1]
        if len(suffix) == 1:
            aorid = '_'.join(aorid.split('_')[0:-1])+'_0%s'%suffix
            #tab['AORID'][idx] = aorid
        aorids.append(aorid)
    tab.replace_column('AORID',Column(list(aorids),name='AORID'))


    coords = [coord.to_string('hmsdms',sep=':',precision=2).split() for coord in coords]
    ra,dec = zip(*coords)
    tab.add_column(Column(ra,name='RA'))
    tab.add_column(Column(dec,name='DEC'))

    pm_ra,pm_dec = zip(*pms)
    tab.add_column(Column(u.Quantity(pm_ra),name='PM_RA'))
    tab.add_column(Column(u.Quantity(pm_dec),name='PM_DEC'))
    
    return tab

def read_POS_file(filename, guide=False):
    '''Parse POS file into table.'''
    with open(filename,'r') as f:
        pos = f.read()
        keys = KEY_RE.findall(pos)[0].split('|')
        if keys[-1] == '':
            del keys[-1]
        keys.insert(1,'Target')

        # each contains (Target,Params,Coord)
        lines = LINES_RE.findall(pos)

        if guide:
            # each contains (Cat,Cam,Radius,MagData,Data)
            guidelines = GUIDE_RE.findall(pos)
            guidelines = [proc_guide_line(line) for line in guidelines]
            guidetab = DataFrame.from_dict(guidelines)
            guidetab = Table.from_pandas(guidetab)
            guidetab.remove_column('NOTSET Mag')
            for col in guidetab.colnames:
                if ' Mag' in col:
                    guidetab[col].unit = u.mag
                if 'PM_' in col:
                    guidetab[col].unit = u.mas/u.yr
            guidetab['Radius'] = u.arcmin
        else:
            guidetab = None

    lines = [proc_line(line) for line in lines if line[1]]
    params,coords,pms = zip(*lines)

    params = [item for sublist in params for item in sublist]
    coords = [item for sublist in coords for item in sublist]
    pms = [item for sublist in pms for item in sublist]

    tab = make_table(params,coords,keys,pms)
    tab.guide = guidetab
    return tab

# Register Table class readers
registry.register_reader('pos-tab', Table, read_POS_file)
