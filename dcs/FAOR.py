#! /usr/bin/env python
import argparse
import xml.etree.ElementTree as ET
from copy import deepcopy
from collections import OrderedDict, deque, defaultdict
from pathlib import Path
import logging
from itertools import takewhile,dropwhile
import sys
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.utils.console import ProgressBar
from functools import partial
from bs4 import BeautifulSoup

XMLPARSER = 'lxml-xml'

KEYS = ['config','run','NODDWELL','REPEATS','loop','endloop','STOP','rewind']
CKEYS = ('FlightPlan','Leg','Leg Dur','Obs Dur','Tot Dur','Start','ROF')

PIX_SCALE = 0.77 # arcsec/pix

# setup logging
log = logging.getLogger(__file__)

class FAOR(object):
    '''Functional AOR for FORCAST'''
    def __init__(self,root,preamble,config,run):
        self.root = root
        self.preamble = preamble
        self.config = config
        self.run = run

    def __repr__(self):
        return 'FAOR(root=%s,\npreamble=%r\n,config=%r\n,run=%r)'%(self.root,self.preamble,self.config,self.run)

    @staticmethod
    def _read_run_block(lines,start):
        '''Read run block from FAOR file'''
        run = deque(takewhile(lambda line: not (('config' in line) or ('#config' in line) or ('EOF' in line)), lines[start+1:]))

        if '# RUN' not in lines[start]:
            run.appendleft(lines[start])
        if not any(['config' in r for r in run]):
            return None

        # backup until no more comment lines
        idx = start - 1
        while idx >= 0:
            if lines[idx][0].lstrip() == '#':
                if '# RUN' not in lines[idx]:
                    run.appendleft(lines[idx])
                idx -= 1
            else:
                break

        # remove lines after config with comments
        run = list(run)
        valid = True
        for idx,line in enumerate(run):
            if 'config' in line:
                valid = False
            if not valid:
                if any([key in line for key in KEYS]):
                    continue
                try:
                    if line.lstrip()[0] == '#':
                        run[idx] = ''
                except IndexError:
                    run[idx] = ''

        run = [line.strip() for line in run]
        run = list(filter(lambda x:x,run))

        return run

    @staticmethod
    def _proc_run_block(run):
        '''Process run block into dictionary'''
        comments = list(takewhile(lambda line: line[0].lstrip() == '#' and 'config' not in line, run))
        comments = ''.join(filter(lambda line:line.strip() != '#' and '# RUN' not in line,comments)) # remove blank comments and join
        r = list(dropwhile(lambda line: line[0].lstrip() == '#' and 'config' not in line, run))
        r = [line.split() for line in r]
        r = filter(lambda x: x, r)
        r = [pair if len(pair) == 2 else (pair[0],'') for pair in r]
        r = OrderedDict(r)
        r['COMMENT'] = comments
        r.move_to_end('COMMENT',last=False)
        if len(r['COMMENT']) > 2 and r['COMMENT'][0:2] == '# ':
            r['COMMENT'] = r['COMMENT'][2:]
        r['COMMENT'] = r['COMMENT'].replace('#','\n#')

        return r

    @classmethod
    def read(cls,aorfile,aorids=None):
        '''Read FAOR from file'''
        with open(aorfile,'r') as f:
            lines = f.readlines()

        # get until first aorid
        preamblelines = takewhile(lambda line: 'AORID' not in line, lines)
        preamble = [line.split() for line in preamblelines]
        preamble = filter(lambda x: x, preamble)
        preamble = [pair if len(pair) == 2 else (pair[0],' '.join(pair[1:])) for pair in preamble]
        preamble = OrderedDict(filter(lambda x: x, preamble))

        # get dur keywords and add to preamble
        idcs = list(filter(lambda x: any([key in x[1] for key in CKEYS]),enumerate(lines)))
        if idcs:
            rem,params = zip(*idcs)
            lines = list(np.delete(np.array(lines),rem))
            for param in params:
                if param[0] == '#':
                    param = param[1:]
                key,val = [x.strip() for x in param.split(': ')]
                preamble[key] = val

        # get until run blocks (or config keyword)
        configlines = list(takewhile(lambda line: not (('# RUN' in line) or ('config' in line) or ('#config' in line)), lines))
        configstartidx = [idx for idx,val in enumerate(configlines) if 'AORID' in val]

        configs = [list(takewhile(lambda line: line.strip() != '',configlines[i:])) for i in configstartidx]
        for idx,config in enumerate(configs):
            cfg = [line.split() for line in config]
            cfg = [pair if len(pair) == 2 else (pair[0],' '.join(pair[1:])) for pair in cfg] # dithers have 2 nums per line
            configs[idx] = OrderedDict(filter(lambda x: x,cfg))

        # run blocks start at # RUN BLOCKS or first 'config'
        runlines = list(dropwhile(lambda line: not (('# RUN' in line) or ('config' in line) or ('#config' in line)), lines))
        runstartidx = [idx for idx,val in enumerate(runlines) if '# RUN' in val or 'config' in val or '#config' in val]

        runs = [FAOR._read_run_block(runlines,start) for start in runstartidx]
        runs = list(filter(lambda x: x is not None,runs))
        runs = [FAOR._proc_run_block(run) for run in runs]

        aor = '_'.join(configs[0]['AORID'].split('_')[0:2])
        root = '_'.join((aor,preamble['TARGET']))
        faor = cls(root,preamble,configs,runs)

        if aorids is not None:
            # only keep specific aorids in list (usually from an obsblk)
            faor.keep_aorids(aorids)
        return faor

    def keep_aorids(self,aorids,comment=False):
        '''Keep config and run blocks with aorids specified.
        If comment, keep run blocks but comment them out.  Default is to delete them.
        '''
        if isinstance(aorids,str):
            aorids = [aorids]
        remove = deque()
        for i,cfg in enumerate(self.config):
            if cfg['AORID'] not in aorids:
                remove.append(i)

        if comment:
            for i in remove:
                for k,v in self.run[i].copy().items():
                    del self.run[i][k]
                    self.run[i]['#%s'%k] = v
        else:
            self.config = list(np.delete(self.config,list(remove)))
            self.run = list(np.delete(self.run,list(remove)))

    def comment_aorids(self,aorids):
        if isinstance(aorids,str):
            aorids = [aorids]
        keep = deque()
        for i,cfg in enumerate(self.config):
            if cfg['AORID'] not in aorids:
                keep.append(i)
        self.keep_aorids(list(keep),comment=True)

    def is_empty(self):
        if self.run and len(self.run) > 0 and any(['run' in run for run in self.run]):
            # the last case checks to make sure at least one run block is active
            return False
        else:
            return True

    def add_loop(self,aorid,numloops):
        '''Add loops to run block'''
        idx = self.index(aorid)
        run = list(self.run[idx].items()) # make aorid run block into tuple
        keys = list(self.run[idx].keys())
        # add endloop after noddwell or repeats
        if 'REPEATS' in keys:
            run.insert(keys.index('REPEATS')+1,('endloop',''))
        else:
            run.insert(keys.index('NODDWELL')+1,('endloop',''))

        # add loop before run cmd
        run.insert(keys.index('run'),('loop',numloops))

        self.run[idx] = OrderedDict(run)

    def modify_run(self,planner):
        '''Modify run blocks with planner table'''
        for row in planner:
            aorid = row['AORID']
            if aorid in self:
                idx = self.index(aorid)
                self.run[idx]['NODDWELL'] = row['noddwell']
                if row['repeats'] > 0:
                    self.run[idx]['REPEATS'] = row['repeats']
                    if self.config[idx]['STYLE'] == 'C2_N_C2':
                        # disallow repeats in C2NC2
                        self.run[idx].pop('REPEATS',None)
                else:
                    # disallow repeats in C2NC2
                    self.run[idx].pop('REPEATS',None)
                if row['rewind'] > 0:
                    self.run[idx]['rewind'] = row['rewind']
                if row['rewind'] == -1:
                    if 'rewind' in self.run[idx]:
                        del self.run[idx]['rewind']
                if row['rewind'] == -2:
                    self.run[idx]['rewind'] = ''
                if row['loops'] > 1:
                    self.add_loop(aorid,row['loops'])

                if 'STOP' in row.colnames and row['STOP'] == True:
                    self.run[idx]['STOP'] = 'STOP'

                # if stop and rewind (blank), make stop first
                if 'STOP' in self.run[idx] and 'rewind' in self.run[idx]:
                    if self.run[idx]['rewind'] == '':
                        self.run[idx].move_to_end('rewind')

                # modify ditherpattern
                if 'DitherPattern' in planner.meta:
                    numdiths = planner.meta['DitherPattern'][aorid]
                    orignumdiths = int(self.config[idx]['DITHER'])
                    if numdiths != orignumdiths:
                        try:
                            scale = np.max([np.abs(float(x)) for x in self.config[idx]['DITHER2'].split()])
                        except KeyError:
                            scale = 10
                        # requested numdiths is different than configured, so need to replace them
                        self.config[idx] = FAOR.modify_dithers(self.config[idx],numdiths,scale)

                for key in CKEYS:
                    if key in planner.meta:
                        self.preamble[key] = planner.meta[key]

                self.run[idx]['COMMENT'] = self.run[idx]['COMMENT'] + '\n#   %s' % row['Comments']

    def __contains__(self,aorid):
        '''Check if aorid is in config'''
        return any([cfg['AORID'] == aorid for cfg in self.config])

    def index(self,aorid):
        '''Return index of aorid in config'''
        return [cfg['AORID'] for cfg in self.config].index(aorid)

    @staticmethod
    def build_dither(num,scale=10):
        '''Build dither pattern'''
        scale = int(scale)
        if num == 2:
            pattern = [[0,0],[scale,scale]]
        elif num == 3:
            pattern = [[0,0],[0,scale],[scale,0]]
        elif num == 4:
            pattern = [[-scale,scale],[-scale,-scale],[scale,-scale],[scale,scale]]
        elif num == 5:
            pattern = [[0,0],[scale,scale],[scale,-scale],[-scale,scale],[-scale,-scale]]
        elif num == 6:
            pattern = [[scale,scale],[scale,0],[scale,-scale],[-scale,scale],[-scale,0],[-scale,scale]]
        elif num == 7:
            pattern = [[0,0],[scale,scale],[scale,0],[scale,-scale],[-scale,scale],[-scale,0],[-scale,scale]]
        elif num == 8:
            pattern = [[scale,scale],[scale,0],[scale,-scale],[0,scale],[0,-scale],[-scale,scale],[-scale,0],[-scale,-scale]]
        elif num == 9:
            pattern = [[0,0],[scale,scale],[scale,0],[scale,-scale],[0,scale],[0,-scale],[-scale,scale],[-scale,0],[-scale,-scale]]
        elif num == 10:
            pattern = [[0,0],[scale,scale],[scale,0],[scale,-scale],[0,scale],[0,-scale],[-scale,scale],[-scale,0],[-scale,-scale],
                       [int(scale/2),int(scale/2)]]
        else:
            raise NotImplementedError('Only dither patterns with 2-10 offsets implemented')

        return pattern

    @staticmethod
    def modify_dithers(cfg,num,scale):
        '''modify configuration to requested number of dithers'''
        # get keys with dither and clear them out
        for key in cfg.copy():
            if 'DITHER' in key:
                del cfg[key]

        if num in (0,1):
            return cfg
        pattern = FAOR.build_dither(num,scale)
        cfg['DITHER'] = str(num)
        for idx,p in enumerate(pattern,start=1):
            p = [str(x) for x in p]
            cfg['DITHER%i'%idx] = '   '.join(p)
        return cfg
        
    def get_comments_as_tables(self):
        '''Return comments in run blocks as tables'''
        tabs = deque()
        for ii,run in enumerate(self.run):
            if run.get('COMMENT'):
                lines = run.get('COMMENT').split('\n')
                # remove # character
                lines = [line[1:] if line[0] == '#' else line for line in lines]

                # first line is name
                meta = {'name':lines[0]}

                # second line contains obs details
                dets = lines[1].split(',')
                dets = {k:v.strip() for k,v in zip(('Mode','Obs Type','Dithers','TLSPN'),dets)}
                
                dets['Dithers'] = int(dets['Dithers'].split()[0])
                ROFrate,span = dets['TLSPN'].split('@')
                ROFrate = u.Quantity(ROFrate,u.deg/u.minute)
                span = span.strip().split()[:2]
                dets['TLSPN'] = u.Quantity(span[0],span[1])
                dets['ROFrate'] = ROFrate

                # get dither scale
                if 'DITHER2' in self.config[ii]:
                    scale = [np.abs(np.rint(float(x))) for x in self.config[ii]['DITHER2'].split()]
                    dets['DITHSCALE'] = int(np.max(scale))
                else:
                    dets['DITHSCALE'] = None

                meta.update(**dets)

                # third line is names
                names = lines[2].split()
                # fifth line contains unit header
                units = lines[3].split()
                # six line is data
                data = lines[5].split()

                table = Table(rows=[data],names=names,dtype=[float]*len(names),meta=meta)
                for col,unit in zip(table.colnames,units):
                    table[col].unit = unit

                table.meta.update(self.preamble)

                tabs.append(table)

            else:
                tabs.append(None)
        return list(tabs)

    def get_comments_as_dicts(self):
        tabs = self.get_comments_as_tables()
        dicts = [{col:t[col][0] for col in t.colnames} for t in tabs]
        for d,t in zip(dicts,tabs):
            d.update(t.meta)
        return dicts
                

        
    def write(self,outfile,ljust=15):
        '''Write out FAOR to file'''
        log.info('Writing FAORs to %s' % outfile)

        with open(outfile,'w') as f:
            # write preamble, skipping empty key:value pairs
            for k,v in self.preamble.items():
                if k in CKEYS:
                    continue
                if v is not None and not v.startswith('\n') and not v == '':
                    f.write('%s%s\n' % (k.ljust(ljust), v))
            f.write('\n')

            # Write configs
            for cfg in self.config:
                for k,v in cfg.items():
                    if k in ('COMMENT','TYPE'):  # skip these
                        continue
                    f.write('%s%s\n' % (k.ljust(ljust), v))
                f.write('\n')

            # Write run blocks
            f.write('# RUN BLOCKS\n')
            for key in CKEYS:
                if key in self.preamble:
                    f.write('#   %s: %s\n' % (key,self.preamble[key]))
            f.write('\n')
                
            for run in self.run:
                for k,v in run.items():
                    if k in ('COMMENT','#COMMENT'):
                        if v:
                            f.write('# %s\n' % v)
                    elif k in ('STOP','#STOP'):
                        f.write('\n%s\n'%k)
                    else:
                        if v:
                            f.write('%s%s\n' % (k.ljust(ljust), v))
                        else:
                            # dont leftpad
                            f.write('%s\n' % k)
                f.write('\n')
            f.write('EOF\n')


def get_keydict(root,keys,exclude=None,as_tag=False):
    if exclude:
        if isinstance(exclude,str):
            exclude = (exclude,)
        keys = list(filter(lambda x:x not in exclude, keys))
    
    try:
        if as_tag:
            return {key:getattr(root,key) for key in keys}
        else:
            return {key:getattr(root,key).text for key in keys}
    except AttributeError:
        vals = (root.find(key) for key in keys)
        if not as_tag:
            vals = (val.text if val else None for val in vals)
        return {key:val for key,val in zip(keys,vals)}
    

# FORCAST faor <- USPOT mapping for Instrument Block
Imap = [  
    ['AORID', 'aorID'],  
    ['INSTMODE', 'InstrumentConfiguration'], 
    ['ORDER', 'order'], 
    ['TOTALEXP', 'TotalTime'], 
    ['CYCLES', 'Repeat'], 
    ['STYLE', 'NodType'], 
    ['CHPTHROW', 'ChopThrow'], 
    ['CHPANG', 'ChopAngle'], 
    ['CHPCOORD', 'ChopAngleCoordinate'], 
    ['NODTHROW', 'NodThrow'], 
    ['NODANG', 'NodAngle'], 
    ['NODCOORD', 'NodAngleCoordinate'], 
    ['DITHER', 'DitherPattern'], 
    ['DITHCOORD', 'DitherCoordinate'], 
    ['SWC', 'InstrumentSpectralElement1'], 
    ['LWC', 'InstrumentSpectralElement2'], 
    ['SLIT', 'Slit'], 
    ['BORESITE', 'undef'], 
    ['IRSRCTYPE', 'IRSourceType']]
    
# FORCAST faor <- USPOT mapping for Target Block
Tmap = [ 
    ['TARGET', 'name'],  
    ['RA', 'lon'],  
    ['DEC', 'lat'],  
    ['PMRA', 'lonPm'],  
    ['PMDEC', 'latPm'],  
    ['COORDSYS', 'equinoxDesc'],  
    ['EPOCH', 'epoch'],  
    # add these to Target Info block
    ['SOURCETYPE', 'SourceType'],  
    ['VMAG', 'VisibleMagnitude'],  
    ['VWAVE', 'VisibleWavelength'],  
    ['IRFLUX', 'IRFlux'],  
    ['IRUNIT', 'IRFluxUnit'],  
    ['IRWAVE', 'IRWavelength']]

#KEYS_LIST = [x for x in zip(*(Tmap+Imap))]
TMAP = dict([t[::-1] for t in Tmap])
IMAP = dict([i[::-1] for i in Imap])
KEYS_MAP = TMAP.copy()
KEYS_MAP.update(IMAP)


def fix_dither(request):
    instrument = request.instrument
    data = instrument.data


    # rename dither keywords            
    dithcoord = data.find('DitherCoordinate')
    if dithcoord and dithcoord.string == 'Array':
        for dRA, dDec in zip(instrument.find_all('deltaRaV'),instrument.find_all('deltaDecW')):
            foo = float(dRA.string)/PIX_SCALE
            dRA.string = str('%.1f' % foo)  #float(dRA.text)/0.76)
            foo = float(dDec.string)/PIX_SCALE
            dDec.string = str('%.1f' % foo)  #float(dDec.text)/0.76)

    # if DITHER exists AND there are no ditherOffsets -> DITHER = '0'
    if data.find('DitherPattern') is not None and \
       instrument.find('ditherOffsets') is None:
        data.find('DitherPattern').string = '0'
        data.find('DitherPattern').name = 'DITHER'
    # otherwise, DITHER equals the number of ditherOffsets found.
    # rename DitherOffset as DITHERN (matches number of dither points)
    else:
        if data.find('DitherPattern') is None:
            dither = BeautifulSoup('<DitherPattern></DitherPattern>',XMLPARSER)
            data.append(dither)
        dpoints = instrument.find('ditherOffsets')
        if dpoints is None:
            data.find('DitherPattern').string = '0'
            data.find('DitherPattern').name = 'DITHER'
        else:
            dpoints = dpoints.find_all('DitherOffset')
            data.find('DitherPattern').string = str(len(dpoints))
            data.find('DitherPattern').name = 'DITHER'
            for idx,tag in enumerate(dpoints):
                tag.name = 'DITHER'+str(idx+1)

    return request

def get_target_data(request):

    request = fix_dither(request)

    instrument = request.instrument
    data = instrument.data
    cfg = {v:data.find(k) for k,v in IMAP.items()}
    cfg = {k:v.string if v else None for k,v in cfg.items()}

    cfg['DITHER'] = data.find('DITHER').string

    # get and reorder dithers
    if cfg['DITHER'] not in ('0','',0,None):
        dkeys = ['DITHER%i'%(i+1) for i in range(int(cfg['DITHER']))]
        for k in dkeys:
            off = instrument.ditherOffsets.find(k)
            dRA = off.find('deltaRaV').string
            dDEC = off.find('deltaDecW').string
            cfg[k] = dRA.ljust(8) + dDEC
        lkeys = list(cfg.keys())
        lk1,lk2 = lkeys[0:lkeys.index('DITHCOORD')+1],lkeys[lkeys.index('DITHCOORD'):-1]
        lkeys = lk1+dkeys+lk2
        cfg = {k:cfg[k] for k in lkeys}

        
    # OPEN (or Slit) value for SWC, LWC, SLIT translates to
    # NONE value for SWC, LWC, SLIT
    if cfg['SWC'] == 'OPEN':
        cfg['SWC'] = 'NONE'
    if cfg['LWC'] == 'OPEN':
        cfg['LWC'] = 'NONE'
    if cfg['SLIT'] in ('OPEN','Slit'):
        cfg['SLIT'] = 'NONE'

    # INSTMODE value of IMG_DUAL, IMG_SWC, IMG_LWC  translates to
    # IMAGING_SWC, IMAGING_LWC, IMAGING_DUAL, GRISM_SWC, GRISM_LWC,
    # GRISM_DUAL, GRISM_XD
    instmode = cfg['INSTMODE']
    if instmode == 'IMG_DUAL':
        cfg['INSTMODE'] = 'IMAGING_DUAL'
    if instmode == 'IMG_SWC':
        cfg['INSTMODE'] = 'IMAGING_SWC'
    if instmode == 'IMG_LWC':
        cfg['INSTMODE'] = 'IMAGING_LWC'
    if instmode == 'InstrumentConfiguration':     # GRISM LS
        if cfg['SWC'] == 'NONE':
            cfg['INSTMODE'] = 'GRISM_LWC'
        elif cfg['LWC'] == 'NONE':
            cfg['INSTMODE'] = 'GRISM_SWC'
        else:
            cfg['INSTMODE'] = 'GRISM_DUAL'
        if cfg['SWC'] == 'FOR_XG063':
            cfg['INSTMODE'] = 'GRISM_XD'
        if cfg['SWC'] == 'FOR_XG111':
            cfg['INSTMODE'] = 'GRISM_XD'

    # some STYLE values need to be translated
    if cfg['STYLE'] == 'NPC_CHOP_ALONG_SLIT':
        cfg['STYLE'] = 'Nod_Perp_Chop_CAS'
    if cfg['STYLE'] == 'NPC_NOD_ALONG_SLIT':
        cfg['STYLE'] = 'Nod_Perp_Chop_NAS'

    # BORESITE value determination:
    # SLIT = NONE -> BORESITE = IMAGE
    # SLIT = FOR_LS* -> BORESITE = LSLIT
    # SLIT = FOR_SS24 -> BORESITE = SSLIT
    if cfg['SLIT'] == 'NONE':
        cfg['BORESITE'] = 'IMAGE'
    elif cfg['SLIT'] == 'FOR_SS24':
        cfg['BORESITE'] = 'SSLIT'
    else:
        cfg['BORESITE'] = 'LSLIT'

    # IRSCRTYPE value determination (available in SSpot v2.4.1 and up)
    # only IMAGING_* have no IRSourceType defined in SSpot, so set value
    # to 'Unknown'
    if cfg['IRSRCTYPE'] in(None,'IRSourceType'):
        cfg['IRSRCTYPE'] = 'Unknown'

    # Acquisition AORs: Cycles is always 1, DITHCOORD is always Array,
    # BORESITE is always LSLIT
    if data.ObsPlanConfig.string == 'ACQUISITION':
        cfg['CYCLES'] = '1'
        cfg['DITHCOORD'] = 'Array'
        cfg['BORESITE'] = 'LSLIT'

    # add TYPE key for post-processing FAOR
    cfg['TYPE'] = data.ObsPlanConfig.string

    target = get_keydict(request.target, keys=TMAP.keys())
    target = {v:target[k] for k,v in TMAP.items()}
    extra_data = get_keydict(data, keys=TMAP.keys())
    extra_data = {v:extra_data[k] for k,v in TMAP.items()}

    if request.target['class'] in ('TargetMovingSingle','SofiaTargetMovingSingle'):
        target['RA'] = '0.0'
        target['DEC'] = '0.0'
        target['PMRA'] = 'UNKNOWN'
        target['PMDEC'] = 'UNKNOWN'
        target['COORDSYS'] = 'J2000'
        target['EPOCH'] = '2000.0'
    # combine extra_data with target only if target info is None
    target.update({k:v for k,v in extra_data.items() if target[k] is None})
    
    return target, cfg

def replace_badchar(string):
    """ replace some reserved characters with '_' 
    """
    string = string.replace('/', '_')
    string = string.replace('\\', '_')
    string = string.replace(':', '_')
    string = string.replace(' ', '_')
    return string


def FO_rename_tags(requests):
    """ rename USPOT tags to FORCAST-required tag (parameter) names
        do not ADD or REMOVE any tag/parameter pairs here - that is done in
        cleanAOR
        Input: XML
        Output: XML
    """
    # create list oldkey (USPOT tags) and newkey (FORCAST tags)
    for r in requests:
        for k,v in KEYS_MAP.items():
            try:
                r.instrument.data.find(k).name = v
            except AttributeError:
                continue
    
    for r in requests:
        # If the DITHCOORD keyword is Array, then the DITHERN RA and DEC
        # (deltaRaV and deltaDecW) values are divided by 0.77
        instrument = r.instrument
        data = instrument.data
        dithcoord = data.find('DITHCOORD')
        if dithcoord and dithcoord.string == 'Array':
            for dRA, dDec in zip(instrument.find_all('deltaRaV'),instrument.find_all('deltaDecW')):
                foo = float(dRA.string)/PIX_SCALE
                dRA.string = str('%.1f' % foo)  #float(dRA.text)/0.76)
                foo = float(dDec.string)/PIX_SCALE
                dDec.string = str('%.1f' % foo)  #float(dDec.text)/0.76)

        # if DITHER exists AND there are no ditherOffsets -> DITHER = '0'
        if data.find('DITHER') is not None and \
           instrument.find('ditherOffsets') is None:
            data.find('DITHER').string = '0'
        # otherwise, DITHER equals the number of ditherOffsets found.
        # rename DitherOffset as DITHERN (matches number of dither points)
        else:
            if data.find('DITHER') is None:
                dither = BeautifulSoup('<DITHER></DITHER>',XMLPARSER)
                data.append(dither)
            dpoints = instrument.find('ditherOffsets')
            if dpoints is None:
                data.find('DITHER').string = '0'
            else:
                dpoints = dpoints.find_all('DitherOffset')
                data.find('DITHER').string = str(len(dpoints))
                for idx,tag in enumerate(dpoints):
                    tag.name = 'DITHER'+str(idx+1)

    return requests


def FO_clean_aor(requests):
    """ remove unwanted info extracted from USpot *.aor file;
        reorder tag-value pairs to match FORCAST faor

        input: xml element containing AORs only (no Proposal info, target list)
        output: array [[tag1,val1], ... , [tagM,valM]] needed for FORCAST FAOR
    """

    instr_func = partial(get_keydict,keys=IMAP.values())
    instr = list(map(instr_func,requests))

    # add dithers
    for idx,tup in enumerate(zip(instr,requests)):
        d,r = tup
        if d['DITHER'] not in ('0','',0,None):
            dkeys = ['DITHER%i'%(i+1) for i in range(int(d['DITHER']))]
            for k in dkeys:
                off = r.instrument.ditherOffsets.find(k)
                dRA = off.find('deltaRaV').string
                dDEC = off.find('deltaDecW').string
                d[k] = dRA.ljust(8) + dDEC
            lkeys = list(d.keys())
            lk1,lk2 = lkeys[0:lkeys.index('DITHCOORD')+1],lkeys[lkeys.index('DITHCOORD'):-1]
            lkeys = lk1+dkeys+lk2
            instr[idx] = {k:d[k] for k in lkeys}

        d = instr[idx]
        # OPEN (or Slit) value for SWC, LWC, SLIT translates to
        # NONE value for SWC, LWC, SLIT
        if d['SWC'] == 'OPEN':
            d['SWC'] = 'NONE'
        if d['LWC'] == 'OPEN':
            d['LWC'] = 'NONE'
        if d['SLIT'] in ('OPEN','Slit'):
            d['SLIT'] = 'NONE'

        # INSTMODE value of IMG_DUAL, IMG_SWC, IMG_LWC  translates to
        # IMAGING_SWC, IMAGING_LWC, IMAGING_DUAL, GRISM_SWC, GRISM_LWC,
        # GRISM_DUAL, GRISM_XD
        instmode = d['INSTMODE']
        if instmode == 'IMG_DUAL':
            d['INSTMODE'] = 'IMAGING_DUAL'
        if instmode == 'IMG_SWC':
            d['INSTMODE'] = 'IMAGING_SWC'
        if instmode == 'IMG_LWC':
            d['INSTMODE'] = 'IMAGING_LWC'
        if instmode == 'InstrumentConfiguration':     # GRISM LS
            if d['SWC'] == 'NONE':
                d['INSTMODE'] = 'GRISM_LWC'
            elif d['LWC'] == 'NONE':
                d['INSTMODE'] = 'GRISM_SWC'
            else:
                d['INSTMODE'] = 'GRISM_DUAL'
            if d['SWC'] == 'FOR_XG063':
                d['INSTMODE'] = 'GRISM_XD'
            if d['SWC'] == 'FOR_XG111':
                d['INSTMODE'] = 'GRISM_XD'

        # some STYLE values need to be translated
        if d['STYLE'] == 'NPC_CHOP_ALONG_SLIT':
            d['STYLE'] = 'Nod_Perp_Chop_CAS'
        if d['STYLE'] == 'NPC_NOD_ALONG_SLIT':
            d['STYLE'] = 'Nod_Perp_Chop_NAS'

        # BORESITE value determination:
        # SLIT = NONE -> BORESITE = IMAGE
        # SLIT = FOR_LS* -> BORESITE = LSLIT
        # SLIT = FOR_SS24 -> BORESITE = SSLIT
        if d['SLIT'] == 'NONE':
            d['BORESITE'] = 'IMAGE'
        elif d['SLIT'] == 'FOR_SS24':
            d['BORESITE'] = 'SSLIT'
        else:
            d['BORESITE'] = 'LSLIT'

        # IRSCRTYPE value determination (available in USpot v2.4.1 and up)
        # only IMAGING_* have no IRSourceType defined in USpot, so set value
        # to 'Unknown'
        if d['IRSRCTYPE'] in ('IRSourceType',None):
            d['IRSRCTYPE'] = 'Unknown'

        # Acquisition AORs: Cycles is always 1, DITHCOORD is always Array,
        # BORESITE is always LSLIT
        if r.instrument.data.ObsPlanConfig.string == 'ACQUISITION':
            d['CYCLES'] = '1'
            d['DITHCOORD'] = 'Array'
            d['BORESITE'] = 'LSLIT'

    # Populate Target Block w/ values from AORs
    target_func = partial(get_keydict,keys=TMAP.keys())
    extra_func = partial(get_keydict,keys=TMAP.values())

    targets = map(target_func,(r.target for r in requests))
    targets = [{v:t[k] for k,v in TMAP.items()} for t in targets]
    extra_data = map(extra_func,(r.instrument.data for r in requests))

    for t,r,e in zip(targets,requests,extra_data):
        if r.target['class'] in ('TargetMovingSingle','SofiaTargetMovingSingle'):
            # non-sidereal
            t['RA'] = '0.0'
            t['DEC'] = '0.0'
            t['PMRA'] = 'UNKNOWN'
            t['PMDEC'] = 'UNKNOWN'
            t['COORDSYS'] = 'J2000'
            t['EPOCH'] = '2000.0'
        # combine extra_data with target only if target info is None
        t.update({k:v for k,v in e.items() if t[k] is None})

    return targets[0], instr


def basic_run_block(cfg):
    '''Generate default run block from FAOR cfg block'''

    aorid = cfg['AORID']

    if cfg['INSTMODE'] in ('IMAGING_DUAL','GRISM_DUAL'):
        mode = '%s/%s' % (cfg['SWC'][-4:],cfg['LWC'][-4:]) # F077/F315
    elif cfg['INSTMODE'] in ('IMAGING_SWC','GRISM_SWC'):
        mode = cfg['SWC'][-4:]
    elif cfg['INSTMODE'] in ('IMAGING_LWC','GRISM_LWC'):
        mode = cfg['LWC'][-4:]
    else:
        mode = ''

    # add slit if present
    if cfg['SLIT'] not in ('NONE','None',None):
        mode = '%s, %s' % (mode,cfg['SLIT'][-4:])

    # add boresite
    mode = '%s | %s boresite' % (mode, cfg['BORESITE'])


    if cfg['TYPE'] == 'ACQUISITION':
        mode = '%s - Acq - %s' % (aorid,mode)
        run = [('COMMENT',mode),
               ('config',aorid),
               ('run',aorid),
               ('NODDWELL','10'),
               ('REPEATS','2'),
               ('STOP','STOP')]
    else:
        mode = '%s - %s' % (aorid,mode)
        run = [('COMMENT',mode),
               ('config',aorid),
               ('run',aorid),
               ('NODDWELL','30'),
               ('REPEATS','1'),
               ('rewind','')]

    run = OrderedDict(run)
    
    if cfg['STYLE'] == 'C2_N_C2':
        # remove repeats
        run.pop('REPEATS',None)
    
    return run

def look_ahead(configs, runs):
    '''Look ahead and add stops if the next config switches boresight'''
    boresites = [cfg['BORESITE'] for cfg in configs]
    for idx,bore in enumerate(boresites):
        if idx+1 == len(boresites):
            continue
        if 'SLIT' in bore and boresites[idx+1] == 'IMAGE':
            runs[idx]['STOP'] = 'STOP'

    return runs

def proc_aorid(combo, requests, PropID):
    '''Process aorids in configs'''

    inst = combo[1]
    if inst != 'FORCAST':
        #raise NotImplementedError('Only FORCAST is supported')
        return None,None,None,None


    # remove non-relevant AORs
    vector = list(filter(lambda r:(r.target.find('name').string,r.instrument.data.InstrumentName.string) == combo,
                         requests))

    # create root, PropID_Target_Inst
    root = PropID + '_' + combo[0] # + '_' + combo[1]  ### don't need _FORCAST
    # replace some reserved characters with '_'
    root = replace_badchar(root)
    log.info('Processing AORID %s' % root)


    data = [(int(r.instrument.data.order.string),r) for r in vector]
    data.sort(key=lambda x: x[0])  # sort by priority only
    # insert the last item (second item; element) from each tuple
    vector = [d[1] for d in data]

    # get the AORID of each AOR
    aorID_ObsConfig = dict((r.instrument.data.aorID.string,r.instrument.data.ObsPlanConfig.string) for r in vector)

    vector = FO_rename_tags(vector)

    preamble, configs = FO_clean_aor(vector)

    log.info('Generated preamble for %s' % preamble['TARGET'])
    log.info('Generated configs for %s' % ', '.join([cfg['AORID'] for cfg in configs]))

    # add config type to dict for post-processing
    for cfg in configs:
        if cfg['AORID'] in aorID_ObsConfig:
            cfg['TYPE'] = aorID_ObsConfig[cfg['AORID']]

    # generate basic run blocks
    runs = [basic_run_block(cfg) for cfg in configs]
    runs = look_ahead(configs,runs)
    log.info('Generated run blocks for %s' % ', '.join([cfg['AORID'] for cfg in configs]))
    
    return root, preamble, configs, runs


def AOR_to_FAORdicts(aorfile, aorids=None, comment=False):
    '''Read FORCAST AORs from .aor file
    If comment, keep run blocks but comment them out.  Default is to delete them.
    '''

    # parse input file
    log.info('Parsing %s' % aorfile)

    with open(aorfile,'r') as f:
        soup = BeautifulSoup(f.read(),XMLPARSER)

    requests = soup.find_all('Request')

    # get proposal id from header
    PropID = soup.AORs.list.ProposalInfo.ProposalID
    try:
        PropID = PropID.string
    except AttributeError:
        PropID = None
    if PropID in (None,''):
        PropID = "00_0000"    # indicates no PropID


    # group by targets
    groups = defaultdict(deque)
    for r in requests:
        target = r.target.find('name').string
        instrument = r.instrument.data.InstrumentName.string
        if instrument != 'FORCAST':
            continue

        # extract instrument data
        preamble, cfg = get_target_data(r)
        if preamble is None:
            continue
        groups[target].append((preamble,cfg))

    faors = deque()
    for _,tup in groups.items():
        preambles,cfgs = zip(*tup)
        preamble = preambles[0]
        cfgs = sorted(cfgs,key=lambda cfg:int(cfg['ORDER']))
        root = replace_badchar('_'.join((PropID,preamble['TARGET'])))

        # generate basic run blocks
        runs = [basic_run_block(cfg) for cfg in cfgs]
        runs = look_ahead(cfgs,runs)
        faors.append((root,preamble,cfgs,runs))

    faors = (FAOR(root,preamble,config,run) for root,preamble,config,run in faors)
    faors = list(filter(lambda faor: faor.preamble is not None, faors))

    # only keep faors for selected aorids
    if aorids is not None:
        if isinstance(aorids,str):
            aorids = [aorids]
        # keep only aorids in config and run specified
        [faor.keep_aorids(aorids,comment=comment) for faor in faors]

    return faors


def main():
    parser = argparse.ArgumentParser(description='AOR to FAOR translator')
    parser.add_argument('aors',nargs='+',help='.aor file(s) to process')
    parser.add_argument('-o',metavar='outdir',type=str,default='.',help="Output directory (default='.')")

    args = parser.parse_args()
    odir = Path(args.o)

    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(message)s')

    for aorfile in args.aors:
        # Generate aors
        faors = AOR_to_FAORdicts(aorfile)

        # Write out basic FAORs
        for faor in faors:
            # Make output directory
            outdir = odir/Path('_'.join(faor.root.split('_')[0:2])) # PropID
            outdir.mkdir(exist_ok=True,parents=True)

            outfile = str((outdir/faor.root).with_suffix('.basic.faor')).replace(',','_')
            faor.write(outfile)

if __name__ == '__main__':
    main()
