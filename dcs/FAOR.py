#! /usr/bin/env python
import argparse
import xml.etree.ElementTree as ET
from copy import deepcopy
from collections import OrderedDict, deque
from pathlib import Path
import logging
from itertools import takewhile,dropwhile
import sys
import numpy as np
import astropy.units as u
from astropy.table import Table

KEYS = ['config','run','NODDWELL','REPEATS','loop','endloop','STOP','rewind']
CKEYS = ('Flight','Leg','Leg Dur','Obs Dur','Tot Dur','Start','ROF')

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


        '''
        for idx,run in enumerate(runs):
            comments = list(takewhile(lambda line: line[0].lstrip() == '#' and 'config' not in line, run))
            comments = ''.join(filter(lambda line:line.strip() != '#' and '# RUN' not in line,comments)) # remove blank comments and join
            r = list(dropwhile(lambda line: line[0].lstrip() == '#' and 'config' not in line, run))
            r = [line.split() for line in r]
            r = filter(lambda x: x, r)
            r = [pair if len(pair) == 2 else (pair[0],'') for pair in r]
            r = OrderedDict(r)
            r['COMMENT'] = comments
            r.move_to_end('COMMENT',last=False)
            runs[idx] = r

        print(runs)
        print(len(runs))
        exit()
        '''
        
        '''
        runs = []
        skip = False
        
        for j,idx in enumerate(runstartidx):
            # starting at each config line, go until next config line
            if skip:
                skip = False
                continue
            
            run = list(takewhile(lambda line: not (('config' in line) or ('#config' in line) or ('EOF' in line)),runlines[idx+1:]))
            if any(['config' in run]):
                # valid run
                runs.append(run)
            else:
                # go until next runstart
                try:
                    run = runlines[idx:runstartidx[j+1]-1]
                    if any(['config' in run]):
                        runs.append(run)
                        skip = True
                        continue
                    else:
                        # likely '# RUN BLOCK'
                        run = runlines[idx:runstartidx[j+2]-1]
                        runs.append(run)
                        skip = True
                        continue
                except IndexError:
                    # this is the last one
                    run = list(takewhile(lambda line: 'EOF' not in line, runlines[idx:]))
                runs.append(run)

        for idx,run in enumerate(runs):
            comments = list(takewhile(lambda line: line[0].lstrip() == '#' and 'config' not in line, run))
            comments = ''.join(filter(lambda line:line.strip() != '#' and '# RUN' not in line,comments)) # remove blank comments and join
            r = list(dropwhile(lambda line: line[0].lstrip() == '#' and 'config' not in line, run))
            r = [line.split() for line in r]
            r = filter(lambda x: x, r)
            r = [pair if len(pair) == 2 else (pair[0],'') for pair in r]
            r = OrderedDict(r)
            r['COMMENT'] = comments
            r.move_to_end('COMMENT',last=False)
            runs[idx] = r

        # duplicates can come through 
        _,unique_runs = np.unique([run.get('#config') if run.get('#config') else run['config'] for run in runs],return_index=True)
        runs = list(np.array(runs)[sorted(unique_runs)])
        [print(x) for x in runs]
        print(len(runs))
        exit()
        '''

        aor = '_'.join(configs[0]['AORID'].split('_')[0:2])
        root = '_'.join((aor,preamble['TARGET']))
        faor = cls(root,preamble,configs,runs)

        if aorids is not None:
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


def replace_badchar(string):
    """ replace some reserved characters with '_' 
    see http://en.wikipedia.org/wiki/Filename for a list
    """
    string = string.replace('/', '_')
    string = string.replace('\\', '_')
    string = string.replace(':', '_')
    string = string.replace(' ', '_')
    return string

def FO_rename_tags(vector):
    """ rename USPOT tags to FORCAST-required tag (parameter) names
        do not ADD or REMOVE any tag/parameter pairs here - that is done in
        cleanAOR
        Input: XML
        Output: XML
    """
    # create list oldkey (USPOT tags) and newkey (FORCAST tags)
    newkey,oldkey = [list(x) for x in zip(*(Tmap+Imap))]

    # rename oldkey to newkey
    for idx, tag in enumerate(oldkey):
        for element in vector.iter(tag):
            element.tag = newkey[idx]

    requests = vector.findall('Request')

    for item in requests:
        # If the DITHCOORD keyowrd is Array, then the DITHERN RA and DEC
        # (deltaRaV and deltaDecW) values are divided by 0.77
        if item.findtext('instrument/data/DITHCOORD') == 'Array':
            for dRA, dDec in zip(item.iter('deltaRaV'), item.iter('deltaDecW')):
                foo = float(dRA.text)/0.77
                dRA.text = str('%.1f' % foo)  #float(dRA.text)/0.76)
                foo = float(dDec.text)/0.77
                dDec.text = str('%.1f' % foo)  #float(dDec.text)/0.76)

        # if DITHER exists AND there are no ditherOffsets -> DITHER = '0'
        if item.find('instrument/data/DITHER') is not None and \
        item.find('instrument/ditherOffsets') is None:
            item.find('instrument/data/DITHER').text = '0'
        # otherwise, DITHER equals the number of ditherOffsets found.
        # rename DitherOffset as DITHERN (matches number of dither points)
        else:
            for dpoints in item.iter('ditherOffsets'):
                item.find('instrument/data/DITHER').text = str(len(dpoints))
                for i in range(len(dpoints)):
                    dpoints.find('DitherOffset').tag = 'DITHER'+str(i+1)

    return vector


def FO_clean_aor(vector):
    """ remove unwanted info extracted from SSpot *.aor file;
        reorder tag-value pairs to match FORCAST faor

        input: xml element containing AORs only (no Proposal info, target list)
        output: array [[tag1,val1], ... , [tagM,valM]] needed for FORCAST FAOR
    """

    Target = {}    # dictionary
    Instrument = {}
    for idx,m in enumerate(vector):  # create dictionary with length
        Instrument[idx] = {}              # equal to number of AORs

    # Populate Inst Block w/ values from AORs (add DITHERN if DITHER > 0)
    newkey,oldkey = [list(x) for x in zip(*Imap)]

    for inst in Instrument:    # fill dictionary w/ dummy values
        for n, o in zip(newkey, oldkey):
            Instrument[inst][n] = o

    # replace dummy w/ real values, add DITHERN RA DEC
    requests = vector.findall('Request')
    for (iter, item) in enumerate(requests):     # loop over AORs
        for n in newkey:
            if item.find('instrument/data/' + n) is not None:
                #print item.find('instrument/data/'+n).text
                # replace dummy w/ real value
                Instrument[iter][n] = item.find('instrument/data/' + n).text
        # some Instruments lack DitherPattern tag in the AOR
        if item.find('.//DITHER') is None:
            Instrument[iter]['DITHER'] = '0'
        # for non-zero DITHER
        elif int(item.find('instrument/data/DITHER').text) > 0:
            for i in range(int(item.find('instrument/data/DITHER').text)):
                dRA = item.find('instrument/ditherOffsets/DITHER' +
                      str(i+1) + '/deltaRaV').text
                dDEC = item.find('instrument/ditherOffsets/DITHER' +
                      str(i+1) + '/deltaDecW').text
                Instrument[iter]['DITHER' + str(i+1)] = dRA.ljust(8) + dDEC

        # OPEN (or Slit) value for SWC, LWC, SLIT translates to
        # NONE value for SWC, LWC, SLIT
        if Instrument[iter]['SWC'] == 'OPEN':
            Instrument[iter]['SWC'] = 'NONE'
        if Instrument[iter]['LWC'] == 'OPEN':
            Instrument[iter]['LWC'] = 'NONE'
        if (Instrument[iter]['SLIT'] == 'OPEN' or
                Instrument[iter]['SLIT'] == 'Slit'):
            Instrument[iter]['SLIT'] = 'NONE'

        # INSTMODE value of IMG_DUAL, IMG_SWC, IMG_LWC  translates to
        # IMAGING_SWC, IMAGING_LWC, IMAGING_DUAL, GRISM_SWC, GRISM_LWC,
        # GRISM_DUAL, GRISM_XD
        instmode = Instrument[iter]['INSTMODE']
        if instmode == 'IMG_DUAL':
            Instrument[iter]['INSTMODE'] = 'IMAGING_DUAL'
        if instmode == 'IMG_SWC':
            Instrument[iter]['INSTMODE'] = 'IMAGING_SWC'
        if instmode == 'IMG_LWC':
            Instrument[iter]['INSTMODE'] = 'IMAGING_LWC'
        if instmode == 'InstrumentConfiguration':     # GRISM LS
            if Instrument[iter]['SWC'] == 'NONE':
                Instrument[iter]['INSTMODE'] = 'GRISM_LWC'
            elif Instrument[iter]['LWC'] == 'NONE':
                Instrument[iter]['INSTMODE'] = 'GRISM_SWC'
            else:
                Instrument[iter]['INSTMODE'] = 'GRISM_DUAL'
            if Instrument[iter]['SWC'] == 'FOR_XG063':
                Instrument[iter]['INSTMODE'] = 'GRISM_XD'
            if Instrument[iter]['SWC'] == 'FOR_XG111':
                Instrument[iter]['INSTMODE'] = 'GRISM_XD'

        # some STYLE values need to be translated
        if Instrument[iter]['STYLE'] == 'NPC_CHOP_ALONG_SLIT':
            Instrument[iter]['STYLE'] = 'Nod_Perp_Chop_CAS'
        if Instrument[iter]['STYLE'] == 'NPC_NOD_ALONG_SLIT':
            Instrument[iter]['STYLE'] = 'Nod_Perp_Chop_NAS'

        # BORESITE value determination:
        # SLIT = NONE -> BORESITE = IMAGE
        # SLIT = FOR_LS* -> BORESITE = LSLIT
        # SLIT = FOR_SS24 -> BORESITE = SSLIT
        if Instrument[iter]['SLIT'] == 'NONE':
            Instrument[iter]['BORESITE'] = 'IMAGE'
        elif Instrument[iter]['SLIT'] == 'FOR_SS24':
            Instrument[iter]['BORESITE'] = 'SSLIT'
        else:
            Instrument[iter]['BORESITE'] = 'LSLIT'

        # IRSCRTYPE value determination (available in SSpot v2.4.1 and up)
        # only IMAGING_* have no IRSourceType defined in SSpot, so set value
        # to 'Unknown'
        if Instrument[iter]['IRSRCTYPE'] == 'IRSourceType':
            Instrument[iter]['IRSRCTYPE'] = 'Unknown'

        # Acquisition AORs: Cycles is always 1, DITHCOORD is always Array,
        # BORESITE is always LSLIT
        if item.find('instrument/data/ObsPlanConfig').text == 'ACQUISITION':
            Instrument[iter]['CYCLES'] = '1'
            Instrument[iter]['DITHCOORD'] = 'Array'
            Instrument[iter]['BORESITE'] = 'LSLIT'

    # Populate Target Block w/ values from AORs
    newkey,oldkey = [list(x) for x in zip(*Tmap)]

    for n, o in zip(newkey, oldkey):  # fill dictionary w/ dummy values
        Target[n] = o

    # Non-sidereal targets: fill in dummy coordinate-related values
    TargetType=vector.find('Request/target')
    if TargetType.get('class') == 'TargetMovingSingle':
        Target['RA'] = '0.0'
        Target['DEC'] = '0.0'
        Target['PMRA'] = 'UNKNOWN'
        Target['PMDEC'] = 'UNKNOWN'
        Target['COORDSYS'] = 'J2000'
        Target['EPOCH'] = '2000.0'

    for n in newkey:
        if requests[0].find('.//'+n) is not None:    # get info from first AOR
            Target[n]=item.find('.//'+n).text   # replace dummy w/ real

    # dictionary -> list
    targetinfo = deepcopy(Tmap)
    for pair in targetinfo:                # Target Block info
        pair[1] = Target[pair[0]]

    instinfo = [] #cleaned = []
    for iter in Instrument:          # Instrument Block info
        temp = deepcopy(Imap)
        for pair in temp:
        	# don't want deltaRaV or deltaDecW, so skip them
            if pair[0] != 'deltaRaV' or pair[0] != 'deltaDecW':
                pair[1] = deepcopy(Instrument[iter][pair[0]])
        # add DITHERN RA DEC here
        if int(Instrument[iter]['DITHER']) > 0:
            for i in range(int(Instrument[iter]['DITHER'])):
                something = ['DITHER' + str(i+1),
                            Instrument[iter]['DITHER' + str(i+1)]]
                # DITHERN goes between DITHCOORD and SWC
                temp.insert(-5, something)
        instinfo.append(temp)

    return targetinfo, instinfo


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

def proc_aorid(aorfile, PropID, combo):
    '''Process aorids in configs'''
    tree = ET.parse(aorfile)
    vector = tree.find('list/vector')
    requests = vector.findall('Request')

    # remove non-relevant AORs
    for elem in requests:
        name = elem.find('target/name').text
        inst = elem.find('instrument/data/InstrumentName').text
        if (name, inst) != combo:
            vector.remove(elem)

    # create root, PropID_Target_Inst
    root = PropID + '_' + combo[0] # + '_' + combo[1]  ### don't need _FORCAST
    # replace some reserved characters with '_'
    root = replace_badchar(root)
    log.info('Processing AORID %s' % root)

    # call the post-processing function for this instrument
    #   currently only FORCAST is supported
    inst = combo[1]
    if inst != 'FORCAST':
        #raise NotImplementedError('Only FORCAST is supported')
        return None,None,None,None

    data = [(int(elem.findtext('instrument/data/order')),elem) for elem in vector]
    data.sort(key=lambda x: x[0])  # sort by priority only
    # insert the last item (second item; element) from each tuple
    vector[:] = [item[-1] for item in data]

    # get the AORID of each AOR
    aorID_ObsConfig = [
        [(item.text) for item in
         vector.findall('Request/instrument/data/aorID')],
        [(item.text) for item in
         vector.findall('Request/instrument/data/ObsPlanConfig')]]

    vector = FO_rename_tags(vector)

    clean_aor = FO_clean_aor(vector)

    preamble = OrderedDict(clean_aor[0])
    log.info('Generated preamble for %s' % preamble['TARGET'])
    configs = [OrderedDict(cfg) for cfg in clean_aor[1]]
    log.info('Generated configs for %s' % ', '.join([cfg['AORID'] for cfg in configs]))

    # add config type to dict for post-processing
    for aorid,obstype in zip(*aorID_ObsConfig):
        for cfg in configs:
            if cfg['AORID'] == aorid:
                cfg['TYPE'] = obstype

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
    tree = ET.parse(aorfile)

    vector = tree.find('list/vector')

    # Extract Target from each AOR
    targets = [(item.text) for item in vector.findall('Request/target/name')]
    # Extract Instrument from each AOR
    instruments = [(item.text) for item in vector.findall('Request/instrument/data/InstrumentName')]

    # Unique combinations of target and instrument
    target_inst = list(set(zip(targets,instruments)))

    # get Proposal ID
    PropID = tree.find('list/ProposalInfo/ProposalID').text
    if PropID == None:
        PropID = "00_0000"    # indicates no PropID

    # Loop over Target-Instrument combo
    faors = (proc_aorid(aorfile, PropID, combo) for combo in target_inst)
    faors = [FAOR(root,preamble,config,run) for root,preamble,config,run in faors]
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

    args = parser.parse_args()

    #logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(message)s')

    for aorfile in args.aors:
        # Generate aors
        faors = AOR_to_FAORdicts(aorfile)

        # Write out basic FAORs
        for faor in faors:
            # Make output directory
            outdir = Path('_'.join(faor.root.split('_')[0:2])) # PropID
            outdir.mkdir(exist_ok=True)

            outfile = (outdir/faor.root).with_suffix('.basic.faor')
            faor.write(outfile)

if __name__ == '__main__':
    main()
    #f = FAOR.read('../AOR_translator/201906_FO_v6/DEBORAH/faors/Leg07__07_0053_HR_5171A.faor')
    #f.config[0] = FAOR.modify_dithers(f.config[0],10,10)
    #print(f.config[0])
    #f = FAOR.read('/home/gordon/Projects/AOR_translator/201906_FO_v0/DAPHNE/faors/Leg06__07_0155_NGC3603-f3-1.faor')
    #f = FAOR.read('/home/gordon/Projects/FAOR/OTTO/06_0162/06_0162_CYGXN46_FORCAST.basic.faor')
    #f = FAOR.read('/home/gordon/Projects/AOR_translator/201906_FO_v0/DAPHNE/faors/Leg09__07_0049_HD_100546_(1).faor')
    #print(f.get_comments_as_tables())
    #f.write('temp.faor')
