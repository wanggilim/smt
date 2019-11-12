#! /usr/bin/env python
import re
from astropy.table import Table,Column
from astropy.io import registry
import numpy as np
from collections import OrderedDict, deque
from astropy.coordinates import SkyCoord
import astropy.units as u

LEG_RE = re.compile('(Leg\s\d.+?)\n\n\n', re.S)
LEG_APPROACH_RE = re.compile('(Leg\s\d*\s\(Approach.+?)\n\n', re.S)

LEG_TAB_RE = re.compile('(UTC\s*MHdg.+?)(?:(?:Leg)|\Z)', re.S)

LEG_NAME_RE = re.compile('(Leg\s\d.+?)\s(?:Start)', re.S)
LEG_NUM_RE = re.compile('Leg\s(\d+)\s')

THDG_RATE_RE = re.compile('THdg:.*(rate:)')
def rate_repl(match):
    return(match.group(0).replace('rate','rate2'))

TOP_META_RE = re.compile('Filename:\s(.*)\sSaved:\s(.*)')
MISSION_META_RE = re.compile('Mission\sSummary.*(Flight.+?)Leg\s1', re.S)

SUA_RE = re.compile('\*\**\s?Potential\sSUA.*\nZone.*', re.S)

#LEG_START_RE = re.compile('Start:\s(.*?)(?:\n\n)', re.S)

META_KEYS = ['Flight Plan ID','Start','Leg Dur','Alt.','ObspID','Blk','Priority','Obs Dur',
             'Target','RA','Dec','Equinox','Elev','ROF','rate','FPI','Moon Angle',
             'Moon Illum','THdg','rate2','Init Hdg','Sun Az Delta','Runway','End Lat','End Lon',
             'End lat','End lon',
             'Sunrise','Sunrise Az','Airport','NAIF ID','Sunset','Sunset Az',
             'Comment','DCS comments','Note','neighbors',
             'Legs','Mach','Takeoff','Obs Time','Flt Time','Landing']

SPLIT_CHECK = 'XXXXXXX'

def convert_column_dtypes(table,
                          scols=['ObsBlk','AOR','RA','DEC','Target']):
    '''Ensure that certain columns are processed with appropriate dtypes'''
    for col in scols:
        if col not in table.colnames:
            continue
        table.replace_column(col,Column(table[col],dtype=str))
        
        for row in table:
            if row[col] == 'None':
                row[col] = ''
        
    return table


def extract_keys(leg,keys=META_KEYS):
    '''Extract keys from each leg'''
    params = OrderedDict()
    for key in keys:
        # get all values after key
        kcolon = '%s:'%key
        if kcolon in leg:
            #startcol = leg.find(key) + len(key) + 1 # 1 for colon
            startcol = leg.find(kcolon) + len(kcolon)
            match = leg[startcol:]
        else:
            continue

        # get positions of every other key in match
        cols = np.array([match.find('%s:'%k) for k in keys],dtype=float)
        #cols = np.array([match.find(k) for k in keys],dtype=float)
        # set -1 to nan and find first key after val
        cols[cols==-1] = np.nan
        try:
            ksplit = keys[np.nanargmin(cols)]
        except ValueError:
            # all nans, therefore this is the last key
            params[key] = match.strip()
            continue
        # split on this key
        match = match.split(ksplit)[0]
        params[key] = match.strip()

    return params
    

def read_UTC_tab(utctext):
    '''Reads utc text table into Table object'''

    utclines = utctext.split('\n')

    # Add comment column to headerline
    utclines[0] = '%sComment' % utclines[0]

    # delete blank lines
    utclines = [line for line in utclines if line]

    # get start of columns in fixed width format
    col_starts = [utclines[0].index(name) for name in utclines[0].split()]

    return Table.read(utclines,format='ascii.fixed_width',col_starts=col_starts)

    
def get_legs(filename):
    '''Divide text into legs'''
    with open(filename,'r') as f:
        text = f.read()

        # rate: keyword is specified twice, which brings things
        #  replace it
        text = THDG_RATE_RE.sub(rate_repl,text)
        
        legs = LEG_RE.findall(text)         # both leg preamble and utc tab
        utc_tabs = LEG_TAB_RE.findall(text) # just utc tab

        # remove possible SUAs
        for idx,tab in enumerate(utc_tabs):
            sua = SUA_RE.findall(tab)
            if sua:
                utc_tabs[idx] = tab.replace(sua[0],'')

        # get toplevel metadata
        top = TOP_META_RE.findall(text)[0]
        top = OrderedDict(zip(('Filename','Saved'),top))
        mission = MISSION_META_RE.findall(text)[0]

        

    # remove utc tabs from legs
    legs = [l.replace(u.strip(),'') for l,u in zip(legs,utc_tabs)]
    
    # make sure comments are attached
    for idx,pair in enumerate(zip(legs,utc_tabs)):
        leg,utc = pair
        if 'Comment' in leg:
            continue

        # check for comments by replacing leg header and utc tab text with XXXXXXX
        rep = text.replace(leg.strip(), SPLIT_CHECK)
        rep = rep.replace(utc, SPLIT_CHECK)

        splitrep = rep.split(SPLIT_CHECK)[1].strip()
        if splitrep:
            # comment is present and not captured by LEG_RE
            legs[idx] = '\n'.join((legs[idx], splitrep))

    # convert utc_tabs to tables
    utc_tabs = [Table.read(tab,format='utc-tab') for tab in utc_tabs]

    # get leg metadata
    leg_names = [LEG_NAME_RE.findall(leg)[0].strip() for leg in legs]

    leg_meta = [extract_keys(leg) for leg in legs]
    
    for meta,leg_name in zip(leg_meta,leg_names):
        # add leg name to meta data
        meta['LegName'] = leg_name
        meta['Leg'] = int(LEG_NUM_RE.findall(leg_name)[0])
        meta.move_to_end('LegName',last=False)
        meta.move_to_end('Leg',last=False)

    # finish top level metadata
    summary = extract_keys(mission)
    summary.update(**top)
    summary['summary'] = mission
    
    # add metadata to tables
    for tab,meta,raw in zip(utc_tabs,leg_meta,legs):
        tab.meta = meta
        tab.raw = raw
        tab.meta['summary'] = summary

    return summary, utc_tabs


def read_MIS_file(filename):
    '''Parse MIS file into table'''
    meta, legs = get_legs(filename)

    colnames = ['Leg','LegName','Blk','ObspID','Target','RA','Dec','Start','Obs Dur','Alt.']#,'MIS']
    rows = deque()
    for leg in legs:
        # extract keys from each leg's metadata
        row = [leg.meta.get(col) for col in colnames]

        # get AOR from Blk
        #idx = colnames.index('Blk')
        #aorid = row[idx].split('OB_')[1] if row[idx] else None
        #row.append(aorid)
        #row[idx] = row[idx].split('OB_')[1] if row[idx] else None
        
        #row = [leg.meta.get(col) for col in colnames[:-1]]
        #row.append(leg)
        rows.append(row)

    table = Table(rows=list(rows),names=colnames)
    table.rename_column('Blk','ObsBlk')
    table.rename_column('ObspID','AOR')
    table.rename_column('Dec','DEC')

    # remove rows with empry leg
    idx = np.where(table['Leg'] == None)
    if idx:
        table.remove_rows(idx[0])

    # attach metadata
    table.meta = meta

    # change ra, dec columns
    for row in table:
        if not row['RA']:
            continue
        #print(row)
        coord = SkyCoord(ra=row['RA'],dec=row['DEC'],unit=(u.hourangle,u.deg))
        ra,dec = coord.to_string('hmsdms',sep=':',precision=2).split()
        row['RA'], row['DEC'] = ra, dec

    # Convert columns to appropriate dtypes
    table = convert_column_dtypes(table)

    # add row metadata
    for idx,leg in enumerate(legs):
        try:
            table.meta['Leg%i'%leg.meta['Leg']] = leg.meta
        except KeyError:
            table.meta['Leg%i'%(idx+1)] = leg.meta

    return table
    
# Register Table class readers
registry.register_reader('utc-tab', Table, read_UTC_tab)
registry.register_reader('mis-tab', Table, read_MIS_file)



if __name__ == '__main__':
    tab = Table.read('mis/201807_HA_IAGO.mis',format='mis-tab')
    tab.pprint()

