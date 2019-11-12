#! /usr/bin/env python
import argparse
from bs4 import BeautifulSoup, element
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table,join,Column
from astropy.io import registry
import numpy as np
import DCS
import xmltodict

# allowed FORCAST imaging types for identifying if an AOR table is all FORCAST
FORCAST_IMG_MODES = ('IMAGING','IMAGING_SWC','IMAGING_LWC','IMAGING_DUAL',
                     'GRISM','GRISM_SWC','GRISM_LWC','GRISM_DUAL','GRISM_XD',
                     'ACQUISITION')

def is_integer(n):
    '''check if float can be safely cast to int'''
    return n % 1 == 0

def get_first_tag(elem,tag,as_string=False):
    '''Return first tag value from elem'''
    val = elem.find(tag.lower())

    if val is None:
        return None
    
    if val.string:
        if as_string:
            return val.string
        
        try:
            num = float(val.string)
            if is_integer(num):
                return int(num)
            else:
                return num
        except ValueError:
            return val.string
    else:
        return [s for s in val.stripped_strings]

def get_all_tags(elem,tag):
    '''Return all tag values'''
    vals = elem.find_all(tag.lower())
    return vals

def convert_column_dtypes(table,
                          dcols=['#Iter','Nod Time','Chop Throw','Chop Angle',
                                 'Scan Dur','Dith Offset','Scan Rate'],
                          scols=['AORID','RA','DEC','Mode','Dith Unit','Scan Amp'],
                          icols=['#Iter','Order']):
    '''Ensure that certain columns are processed with appropriate dtypes'''
    for col in dcols:
        if col not in table.colnames:
            continue
        try:
            table.replace_column(col,Column(table[col],dtype=float))
        except ValueError:
            # leave as string values
            continue

    for col in icols:
        if col not in table.colnames:
            continue
        try:
            table.replace_column(col,Column(table[col],dtype=int))
        except TypeError:
            continue


    for col in scols:
        if col not in table.colnames:
            continue
        table.replace_column(col,Column(table[col],dtype=str))
        for row in table:
            if row[col] == 'None':
                row[col] = ''
        
    return table

def get_preamble(aor):
    '''Returns preamble properties (e.g. ProposalID, Title, PI)'''

    tags = ['ProposalID', 'ProposalTitle', 'Category','ScienceKeywords']

    preamble = {tag:str(get_first_tag(aor,tag,as_string=True)) for tag in tags}

    investigators = list(get_all_tags(aor,'Investigator'))

    for idx,inv in enumerate(investigators):
        # convert tag objects to dicts
        investigators[idx] = {k:v for k,v in inv.attrs.items() if v and v != 'NONE'}
    preamble['Investigators'] = investigators

    try:
        preamble['PI'] = investigators[0]
    except IndexError:
        pi = str(get_first_tag(aor,'proposalpi',as_string=True))
        preamble['PI'] = {'firstname':'','lastname':pi,'honorific':'','institution':''}
                
    return preamble


def get_comments(aor):
    '''Returns obsblock comments'''
    blks = get_all_tags(aor,'ObsBlockInfo')
    comments = {}
    for blk in blks:
        comment = blk.comment
        if comment and comment.string:
            comments[blk.obsblockid.string] = comment.string

    return comments

def get_requests(aor):
    '''Returns request objects'''
    requests = get_all_tags(aor,'Request')

    for request in requests:
        request.aorid = request.instrument.aorid.string

        try:
            ra =  float(request.target.position.lon.string)
            dec = float(request.target.position.lat.string)
        except AttributeError:
            # target likely ephemeris
            request.coord = None
            continue
            
        try:
            equ = request.target.position.coordSystem.equinoxDesc.string
        except AttributeError:
            equ = 'J2000'
        request.coord = SkyCoord(ra=ra,dec=dec,equinox=equ,unit=(u.deg,u.deg))

    return requests

def get_obs_mode(request):
    '''Get observation mode from request'''
    if get_first_tag(request,'InstrumentName') == 'FORCAST':
        # check if acquisition
        obsmode = get_first_tag(request,'ObsPlanConfig')

        if obsmode == 'ACQUISITION':
            return obsmode
        elif obsmode == 'GRISM_XD':
            return 'GRISM_XD'
        else:
            return get_first_tag(request,'InstrumentConfiguration')
    
    if get_first_tag(request,'ObsPlanConfig') == 'POLARIZATION':
        return 'Polarimetry'
    if get_first_tag(request,'ScanType') == 'Lissajous':
        return 'Lissajous'
    if get_first_tag(request,'ScanType') == 'Box':
        return 'Box'
    return None

def get_filter(request):
    '''Get filters and boresights from request'''

    if get_first_tag(request,'InstrumentName') == 'FORCAST':
        # FORCAST
        mode = get_first_tag(request,'InstrumentConfiguration')
        if mode in ('IMAGING_DUAL','GRISM_DUAL'):
            filt1 = get_first_tag(request,'InstrumentSpectralElement1').split('_')[1]
            filt2 = get_first_tag(request,'InstrumentSpectralElement2').split('_')[1]
            filt = '%s/%s' % (filt1,filt2)
        elif mode in ('IMAGING_SWC','GRISM_SWC','GRISM_XD'):
            filt = get_first_tag(request,'InstrumentSpectralElement1').split('_')[1]
        elif mode in ('IMAGING_LWC','GRISM_LWC'):
            filt = get_first_tag(request,'InstrumentSpectralElement2').split('_')[1]
        else:
            # likely ACQ
            obsmode = get_first_tag(request,'ObsPlanConfig')
            if obsmode in ('ACQUISITION','GRISM_XD'):
                filt1 = get_first_tag(request,'InstrumentSpectralElement1').split('_')[1]
                try:
                    filt2 = get_first_tag(request,'InstrumentSpectralElement2').split('_')[1]
                except IndexError:
                    filt2 = ''
                filt = filt1 if filt1 != 'OPEN' else filt2
            else:
                filt = ''
    else:
        # HAWC+
        try:
            filt = get_first_tag(request,'InstrumentSpectralElement1').split('_')[1]
        except IndexError:
            return ''

        '''
        bore = '%s_main' % filt.lower()
        if filt == 'C':
            bore = 'bc_main'
        return '%s/%s' % (filt,bore)
        '''


    return filt

def get_chop_coord_type(request):
    '''Get SIRF or ERF from coord type'''
    typescan = get_first_tag(request,'ChopAngleCoordinate')
    if typescan in ['Sky','J2000']:
        return 'ERF'
    elif typescan == 'Array':
        return 'SIRF'
    elif typescan == 'HORIZON':
        return 'HOR'
    else:
        raise ValueError('ChopAngleCoordinate type %s not understood' % typescan)
    

def make_request_table(requests,include_meta=False):
    obs_overview = [(r.aorid,r.title.string,*r.coord.to_string('hmsdms',sep=':',precision=2).split()) if r.coord else (r.aorid,r.title.string,'','') for r in requests]

    obs_overview = Table(rows=obs_overview, names=['AORID','Name','RA','DEC'])

    modes = [get_obs_mode(r) for r in requests]
    filts = [get_filter(r) for r in requests]
    aorids = [r.aorid.string for r in requests]
    angles = [get_first_tag(r,'RollAngle') for r in requests]
    angles = [float(a) if a not in ('',None) else 0 for a in angles]

    obs_overview['_AORID'] = aorids

    # change 90_0086_6 to 90_0086_06
    for idx,aorid in enumerate(aorids):
        suffix = aorid.split('_')[-1]
        if len(suffix) == 1:
            aorids[idx] = '_'.join(aorid.split('_')[0:-1])+'_0%s'%suffix

    obs_overview.replace_column('AORID',aorids)

    obsnames = [r.title.string for r in requests]
    numiters = [get_first_tag(r,'Repeat') for r in requests]

    '''
    if 'Box' in modes or all([mode == 'Polarimetry' for mode in modes]):
        nodtimes = [get_first_tag(r,'NodTime') for r in requests]
        chopthrow = [get_first_tag(r,'ChopThrow') for r in requests]
        chopangle = [get_first_tag(r,'ChopAngle') for r in requests]
        choptype = [get_chop_coord_type(r) for r in requests]

        rows = zip(modes, filts, aorids, obsnames, numiters,
                   nodtimes, chopthrow, chopangle, choptype)
        names = ['Mode','Band/Bore','AORID','Name','#Iter',
                 'Nod Time','Chop Throw','Chop Angle','Sys']


    elif all([mode == 'Lissajous' for mode in modes]):
        scandur = [get_first_tag(r,'TotalTime') for r in requests]
        scanamp = [get_first_tag(r,'ScanAmplitudeEL') for r in requests]
        
        rows = zip(modes, filts, aorids, obsnames, numiters,
                   scandur, scanamp)
        names = ['Mode','Band/Bore','AORID','Name','#Iter',
                 'Scan Dur','Scan Amp']
    else:
        # mixed modes in AOR, get all properties
        nodtimes = [get_first_tag(r,'NodTime') for r in requests]
        chopthrow = [get_first_tag(r,'ChopThrow') for r in requests]
        chopangle = [get_first_tag(r,'ChopAngle') for r in requests]
        choptype = [get_chop_coord_type(r) for r in requests]
        scandur = [get_first_tag(r,'TotalTime') for r in requests]
        scanamp = [get_first_tag(r,'ScanAmplitudeEL') for r in requests]

        rows = zip(modes, filts, aorids, obsnames, numiters,
                   nodtimes, chopthrow, chopangle, scandur, scanamp, choptype)
        names = ['Mode','Band/Bore','AORID','Name','#Iter',
                 'Nod Time','Chop Throw','Chop Angle','Scan Dur','Scan Amp','Sys']
    '''
    nodtimes = [get_first_tag(r,'NodTime') for r in requests]
    chopthrow = [get_first_tag(r,'ChopThrow') for r in requests]
    chopangle = [get_first_tag(r,'ChopAngle') for r in requests]
    choptype = [get_chop_coord_type(r) for r in requests]
    #scandur = [get_first_tag(r,'TotalTime') for r in requests]
    #scandur = [get_first_tag(r,'ScanTime') for r in requests]

    orders = [get_first_tag(r,'Order') for r in requests]

    scandur = []
    scanamp = []
    scanrate = []
    dithunit = []
    dithoffset = []
    for r in requests:
        if get_first_tag(r,'InstrumentName') == 'FORCAST':
            repeats = get_first_tag(r,'repeat')
            time = get_first_tag(r,'TotalTime')
            try:
                repeats = int(repeats)
                time *= repeats
            except (ValueError,TypeError):
                pass
            scandur.append(time)
            sa = get_first_tag(r,'NodType')
            if sa == 'Nod_Match_Chop':
                sa = 'NMC'
            elif sa == 'C2_N_C2':
                sa = 'C2NC2'
            else:
                sa = sa
            scanamp.append(sa)  # for forcast, scanamp is the nodtype

            du = get_first_tag(r,'DitherPattern')
            du = du[0] if du and du != 'None' else ''                      
            dithunit.append(du)  # for forcast, dithunit is the pattern

            do = get_first_tag(r,'ditherOffsets')
            if do:
                do = np.max([float(x) for x in do])
            dithoffset.append(do)  # for forcast, dithoffset is the scale

            scanrate.append(None)
        else:
            scandur.append(get_first_tag(r,'ScanTime'))
            #scanamp.append(get_first_tag(r,'ScanAmplitudeEL'))
            try:
                el = np.float(get_first_tag(r,'ScanAmplitudeEL'))
                xel = np.float(get_first_tag(r,'ScanAmplitudeXEL'))
                if el.is_integer():
                    el = int(el)
                if xel.is_integer():
                    xel = int(xel)

                # if same el and xel, don't add /
                if el == xel:
                    samp = str(el)
                else:
                    samp = '%s/%s'%(el,xel)
            except TypeError:
                samp = ''
            scanamp.append(samp)
            scanrate.append(get_first_tag(r,'ScanRate'))
            dithunit.append(get_first_tag(r,'DitherOffsetUnit'))
            dithoffset.append(get_first_tag(r,'DitherScale'))
    
    rows = zip(modes, filts, aorids, obsnames, numiters,
               nodtimes, chopthrow, chopangle, scandur, scanamp, scanrate, choptype,
               dithoffset,dithunit,orders,angles)
    names = ['Mode','Band','AORID','Name','#Iter',
             'Nod Time','Chop Throw','Chop Angle','Scan Dur','Scan Amp','Scan Rate','Sys',
             'Dith Offset','Dith Unit','Order','Angle']

    obs_details = Table(rows=list(rows),names=names)
    
    obs_table = join(obs_details,obs_overview,keys=['AORID','Name'])

    # if any c2nc2, add nod angles
    if any([sa == 'C2NC2' for sa in scanamp]):
        nodangles = [get_first_tag(r,'NodAngle') for r in requests]
        nodthrow = [get_first_tag(r,'NodThrow') for r in requests]
        nodangles = [float(x) if x not in (None,'None','','--') else None for x in nodangles]
        nodthrow = [float(x) if x not in (None,'None','','--') else None for x in nodthrow]
        obs_table.add_column(Column(nodangles,name='Nod Angle',dtype=float))
        obs_table.add_column(Column(nodthrow,name='Nod Throw',dtype=float))
        obs_table['Nod Angle'].unit = u.deg
        obs_table['Nod Throw'].unit = u.arcsec
        

    # if FORCAST grism is present, add slit if necessary
    if any(np.in1d(obs_table['Mode'],['GRISM_SWC','GRISM_LWC','GRISM_DUAL','ACQUISITION','GRISM_XD'])):
        # add slit column
        slits = [get_first_tag(r,'Slit',as_string=True) for r in requests]
        slits = [s if s not in ('',None,'None','NONE') else '' for s in slits]
        aorids = [r.aorid.string for r in requests]
        slit_tab = Table(rows=list(zip(aorids,slits)),names=('_AORID','Slit'))
        obs_table = join(obs_table,slit_tab,join_type='left',keys=['_AORID'])
        obs_names = list(obs_table.colnames.copy())
        idx = obs_names.index('Slit')
        del obs_names[idx]
        obs_names.insert(2,'Slit')
        obs_table = obs_table[obs_names]
        
        
        #print(requests[0])
        #exit()
        #slits = [r.instrument.data.slit.string for r in requests]
        #exit()
        #slits = [s if s not in ('',None,'None','NONE') else '' for s in slits]
        #obs_table.add_column(Column(slits,name='Slit'),index=2)

    # sort by aorid
    #sort_order = np.argsort([int(aorid.split('_')[-1]) for aorid in obs_table['AORID']])
    #obs_table = obs_table[sort_order]

    # fix column dtypes
    obs_table = convert_column_dtypes(obs_table)

    if include_meta:
        # add metadata for each row
        for row,request in zip(obs_table,requests):
            meta = xmltodict.parse(str(request))['request']
            obs_table.meta[row['AORID']] = meta

    return obs_table

def get_obs_blockids(aor,aorids,d=None):
    '''Get associated obsblockid for given aorid'''
    blocks = get_all_tags(aor,'ObsBlockInfo')

    # for each obsblock, get aorids
    blockdict = {block.obsblockid.string:[a.string \
                                          for a in block.aorlist.find_all('aorid')] \
                 for block in blocks}

    # reverse dictionary, to have aorids mapped to obsblocks
    aordict = {aorid:block for block,aorlist in blockdict.items() for aorid in aorlist}

    '''
    if not aordict:
        # dict is empty--likely a legacy aor with no obsblocks
        obsblock_col = Column(['OB_%s'%aorid for aorid in aorids],name='ObsBlk')
        return obsblock_col
    '''

    if not aordict:
        # dict is empty--likely a legacy aor with no obsblocks
        #   must get obsblocksearch to determine aorids in obsblock
        if d is None:
            d = DCS.DCS()

        # have to remove the zero (_03 -> _3)
        aorstr = [x.split('_') for x in aorids]
        for idx,aor in enumerate(aorstr):
            if aor[-1][0] == '0':
                aorstr[idx] = [aor[0],aor[1],aor[-1][1]]
            aorstr[idx] = '_'.join(aorstr[idx])
        aorsearch = ','.join(aorstr)
        aorsearch = d.getAORSearch(aorids=aorsearch)

        # get indices where aorstr matches since the table isn't sorted
        idx = [np.where(aorsearch['AORID'] == x) for x in aorstr]
        obsblocks = [str(aorsearch['ObsBlockID'][x][0]).strip() for x in idx]
        obsblock_col = Column(obsblocks,name='ObsBlk')
        return obsblock_col

            

    # change 90_0086_6 to 90_0086_06
    for aorid,block in aordict.copy().items():
        suffix = aorid.split('_')[-1]
        if len(suffix) == 1:
            del aordict[aorid]
            aorid = '_'.join(aorid.split('_')[0:-1])+'_0%s'%suffix
            aordict[aorid] = block

    # generate column of associated aorids
    aoridlist = [aordict.get(aorid,'') for aorid in aorids]
    obsblock_col = Column(aoridlist,name='ObsBlk')
    return obsblock_col


def read_AOR_file(filename,obsblock=None,d=None,include_meta=False):
    '''Parse AOR file into table. If obsblock, return only those aorids
       that are in the specified obsblock'''
    with open(filename,'r') as f:
        aor = BeautifulSoup(f,'lxml')

    preamble = get_preamble(aor)
    requests = get_requests(aor)

    comments = get_comments(aor)
    
    obs_table = make_request_table(requests,include_meta=include_meta)

    # attach obsblockids to table
    obs_blocks = get_obs_blockids(aor,obs_table['AORID'],d=d)
    obs_table.add_column(obs_blocks,index=obs_table.index_column('AORID')+1)

    if comments:
        # add obsblk comments
        obsblk,cmt = zip(*[(k,v) for k,v in comments.items()])
        comtable = Table()
        comtable['ObsBlk'] = obsblk
        comtable['ObsBlkComment'] = cmt
        obs_table = join(obs_table,comtable,join_type='left',keys=['ObsBlk'])

    obs_table.meta.update(**preamble)

    if obsblock:
        obs_table = obs_table[obs_table['ObsBlk'] == obsblock]
        #obs_table.comments = comments

    # check if all forcast
    if all(np.in1d(obs_table['Mode'],['',*FORCAST_IMG_MODES])):
        # forcast-specific table
        obs_table.rename_column('Scan Dur','Total Exp')
        obs_table.rename_column('Scan Amp','Obs Type')
        obs_table.rename_column('Nod Time','Nod Dwell')
        obs_table.rename_column('Dith Unit','Num Dith')
        obs_table.remove_column('#Iter')

    obs_table.meta['CFILENAME'] = filename
    #obs_table.meta['ObsBlkComments'] = comments
    
    return obs_table

# Register Table class reader
try:
    registry.register_reader('aor-tab', Table, read_AOR_file)
except registry.IORegistryError:
    pass

def main():
    parser = argparse.ArgumentParser(description='Translate HAWC+ AOR to dossier')
    parser.add_argument('aor',type=str,help='AOR file')

    args = parser.parse_args()

    aor = Table.read(args.aor,format='aor-tab')

    aor.pprint()
    

if __name__ == '__main__':
    main()
