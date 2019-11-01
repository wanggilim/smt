#! /usr/bin/env python
import re
from configparser import ConfigParser
import astropy.units as u
import argparse
import numpy as np
from dcs import DCS
from itime3 import itime3
from scipy.optimize import minimize, basinhopping
from functools import partial
from astropy.table import Table,hstack,vstack,unique,Column
from astropy.utils.console import ProgressBar
from collections import deque
from CMap import CMap
from pathlib import Path
from itertools import cycle
from blessings import Terminal
from dcs import FAOR as dcsFAOR
from datetime import datetime,timedelta,timezone
from astropy.time import Time, TimeDelta
import json

# ROF rate
ROF_RE = re.compile('\[([\+\-\d\.]*)\,\s?([\+\-\d\.]*)\]')

#ALLOWED_INTERVALS = (15,30,45,60,75,90,120)

VALIDSTR = ('force','yes','on','true')


# simple Table wrapper for the sole purpose of bolding the first line
class FTable(Table):

    def pprint(self,fline=0,**kwargs):
        lines, outs = self.formatter._pformat_table(self,**kwargs)
        if outs['show_length']:
            lines.append('Length = {0} rows'.format(len(self)))

        n_header = outs['n_header']

        # init terminal
        term = Terminal()
        
        for i, line in enumerate(lines):
            if i < n_header:
                print('%s%s%s%s'%(term.bold,term.red,line,term.normal))
                #color_print(line, 'red')
            elif i == (n_header + fline):
                print('%s%s%s'%(term.standout,line,term.normal))
            else:
                print(line)
        

def get_ROF_rate(mistab, leg=None, key='rate'):
    '''Get ROF rate for each leg, or specific one.
         leg can be of the form 'LegX' or just 'X'
    '''

    if leg is None:
        # extract metadata for each leg
        num_legs = int(mistab.meta['Legs'])
        legs = [mistab.meta.get('Leg%i'%(i+1)) for i in range(0,num_legs)]
    else:
        # extract specified leg
        if isinstance(leg,str):
            try:
                legno = int(leg)
                leg = 'Leg%i'%legno
            except ValueError:
                # leg is LegX
                leg = leg
        else:
            leg = 'Leg%i'%leg
            
        legs = [mistab.meta.get(leg)]
    
    # get ROFrate string
    rofrate = [leg.get(key) for leg in legs]

    # parse ROF list e.g. '[+0.10, +0.09] deg/min'
    rofrate = [ROF_RE.findall(r)[0] if r else None for r in rofrate]

    # convert to quantity
    rofrate = [(float(r[0]),float(r[1]))*u.deg/u.min if r else None for r in rofrate]

    if leg:
        return rofrate[0]
    else:
        return rofrate


def compute_LOS_cadence(dur, rofrate, fastspan=2.5*u.deg, slowspan=2*u.deg, slowlim=0.25*u.deg/u.min):
    '''Calculate number of LOS rewinds required for leg.
    fastspan is the relaxed allowed rotation angle for fast rotators.
    slowspan is the maximum rotation angle for slow rotators.
    the cut-off between fast and slow is slowlim

    rofrate is the rotation rate for the fov
    dur can be the leg duration or aor duration
    '''
    rofrate = np.abs(rofrate)
    rofrate = u.Quantity(rofrate,slowlim.unit)
    
    # choose span based on rotation rate
    span = slowspan if rofrate < slowlim else fastspan

    if np.isclose(rofrate.value,0):
        try:
            return np.inf*u.s,[0 for x in list(dur)],span
        except TypeError:
            return np.inf*u.s,[0],span
    
    tlos = span/rofrate
    
    nlos = dur/tlos
    # round down

    try:
        nlos = [int(n) for n in list(nlos)]
    except TypeError:
        try:
            nlos = [int(nlos)]
        except TypeError:
            nlos = [0]
            
    return tlos,nlos,span

def round_time(time, resolution=15*u.minute):
    '''Round to nearest 15 minute increment'''
    resolution = resolution.to(u.minute).value

    if time.shape:
        dts = time.to_datetime(timezone.utc)
    else:
        dts = [time.to_datetime(timezone.utc)]

    #new = [(dt.minute // resolution + (1 if direction == 'up' else 0)) * resolution for dt in dts]
    up   = [((dt.minute // resolution) + 1) * resolution for dt in dts]
    down = [((dt.minute // resolution) + 0) * resolution for dt in dts]
    # time past hour
    tph = [timedelta(minutes=dt.minute,seconds=dt.second) for dt in dts]
    updeltas   = [timedelta(minutes=n)-tp for tp,n in zip(tph,up)]
    downdeltas = [tp-timedelta(minutes=n) for tp,n in zip(tph,down)]
    
    #downdeltas = [timedelta(minutes=dt.minute-n,seconds=dt.second-n) for dt,n in zip(dts,down)]

    # choose closest
    new = [dt+u if u < d else dt-d for dt,u,d in zip(dts,updeltas,downdeltas)]
    #new = [dt + timedelta(minutes=n-dt.minute) - timedelta(seconds=dt.second) for dt,n in zip(dts,new)]
    new = [dt.strftime('%Y-%b-%d %H:%M:%S') for dt in new]
    new = Time([Time.strptime(x,'%Y-%b-%d %H:%M:%S') for x in new])

    return new

def split_leg(utctab, interval, rofrate=None):
    '''Split leg into intervals, calculating rof rates for each interval along the leg'''
    if not isinstance(interval,u.Quantity):
        interval = u.Quantity(interval, u.minute)
    else:
        interval = interval.to(u.minute)

    # remove rows where ROFrt == 'N/A'--setup time
    idx = np.where(utctab['ROFrt'] == 'N/A')
    if idx:
        utctab.remove_rows(idx[0])

    date = utctab.meta['summary']['Takeoff'].split()[0]
    times = Time([Time.strptime(' '.join((t,date)),'%H:%M:%S %Y-%b-%d') for t in utctab['UTC']])

    # set start time to closest 15 min interval
    inittime = round_time(times[0])

    # increment by interval until end of leg
    intervals = deque()
    if interval == 0*u.minute:
        # no interval, so don't split
        start,stop = inittime, times[-1]
        rof = utctab['ROF'][0]
        rofr = np.max([np.abs(float(x)) for x in utctab['ROFrt']])
        #label = start.strftime('%H%M')[0]
        label = ''
        intervals.append((start,stop,rof,rofr,label))
    else:
        start,stop = inittime, inittime+interval

    while stop < times[-1]:
        # find times in interval
        tidx = np.where((times >= start) & (times <= stop))
        rof = utctab['ROF'][tidx][0]  # initial ROF in interval
        rofr = np.max([np.abs(float(x)) for x in utctab['ROFrt'][tidx]])  # max ROFrate in interval
        label = start.strftime('%H%M')[0]

        intervals.append((start[0],stop[0],rof,rofr,label))
        start = stop
        stop = start+interval
    else:
        intervals = list(intervals)
        # set first interval to beginning
        stop = intervals[0][1]
        label = intervals[0][-1]
        del intervals[0]
        tidx = np.where((times <= stop))
        rof = utctab['ROF'][tidx][0]
        rofr = np.max([np.abs(float(x)) for x in utctab['ROFrt'][tidx]])
        intervals.insert(0,(times[0],stop,rof,rofr,label))
        
        # set last interval stop to end of times
        start = intervals[-1][0]
        label = intervals[-1][-1]
        del intervals[-1]
        tidx = np.where((times >= start))
        rof = utctab['ROF'][tidx][0]
        rofr = np.max([np.abs(float(x)) for x in utctab['ROFrt'][tidx]])
        intervals.append((start,times[-1],rof,rofr,label))
    

    start,stop,ROF,ROFrt,labels = zip(*intervals)
    #start = Time([x[0] for x in start])
    #stop  = Time([x[0] for x in stop])

    if rofrate:
        # override rofrate
        ROFrt = [rofrate]*len(ROFrt)

    intervals = Table()
    intervals['start'] = Time(start)
    intervals['stop'] = Time(stop)
    intervals['ROF'] = ROF*u.deg
    intervals['ROFrt'] = ROFrt*u.deg/u.min
    intervals.meta['interval'] = interval

    # get los timing just for display purposes in --dry-run mode
    tlos,_,span = zip(*[compute_LOS_cadence(0,rt) for rt in intervals['ROFrt']])
    intervals['TLOS'] = u.Quantity(tlos,u.s)
    intervals['TLOS'].format = '%.1f'
    intervals['TLSPAN'] = u.Quantity(span,u.deg)
    intervals['label'] = labels
    
    return intervals

def optimize_aor(tup, maximize, tlos, taor, los_buffer=30, maxfrac=.04,**kwargs):
    '''Optimization function for itime3'''
    # maxfrac is used in inttime optimization mode. it reflects what fraction of the integration time
    #  over the request is allowed

    if any([np.isnan(x) for x in tup]):
        return np.inf

    noddwell, repeats, rewinds, dithers, loops = tup
    # repeats,rewinds,loops need to be integers
    repeats = int(np.rint(repeats))
    dithers = int(np.rint(dithers))
    rewinds = int(np.rint(rewinds))
    loops = int(np.rint(loops))
    ## ALSO NODDWELL??
    #noddwell = int(np.rint(noddwell))
    

    if 'c2nc2' in kwargs:
        if kwargs['c2nc2'] and dithers % 2 != 0:
            # c2nc must have even dithers
            return np.inf

    if repeats > 0:
        '''
        if repeats > 2 and (repeats % 2 == 1):
            # repeats must be even for c2n
            return np.inf
        '''
        if (rewinds % repeats != 0) and (rewinds != 0):
            # Repeats must be an integer multiple of rewinds
            return np.inf
    durperrew,fdur,inttime,eff = itime3(noddwell,repeats,rewinds,dithers,loops=loops,**kwargs)

    if durperrew > (tlos+los_buffer):
        # duration per rewind needs to be less than rewind cadence (with a buffer)
        return np.inf

    # minimize this function so return 1/maximize
    if maximize == 'eff':
        return 1/eff
    elif maximize == 'inttime':
        if inttime > (1+maxfrac)*taor:
            return np.inf
        return 1/inttime
    elif maximize == 'dur':
        return np.abs(taor-fdur)
    else:
        raise ValueError('Allowed values for maximize are (eff,inttime,dur)')

def do_basin(inp,func,x0,niter=1000,stepsize=4,interval=50,bounds=None,fixed=False):
    '''This actually runs the basinhopping algorithm'''
    seed,maximize = inp
    if fixed:
        return {'x':x0}
    res = basinhopping(func, x0, niter=niter,stepsize=stepsize,interval=interval,seed=seed,
                       minimizer_kwargs={'args':(maximize,),'bounds':bounds})
    return res
    

def sort_and_eval(table,eff=False):
    '''Perform final selection of configuration based on efficiency, inttime, and FDUR'''
    if eff:
        S = 1/table['EFF']
    else:
        S = 1/table['EFF'] + np.abs(table['FDUR']-table['TAOR'])
    idx = np.argsort(S)
    return table[idx]

def read_rofdict(rofstr):
    '''Process config 'rofdict' key'''
    #keys = ('noddwell','dithers','repeats','rewinds','loops')
    #for key in keys:
    #    rofstr = rofstr.replace(key,'"%s"'%key)
    #print(rofstr)
    #rofdict = json.loads(rofstr.replace('\n',''))
    #exit()
    lines = rofstr.split('\n')
    for idx,line in enumerate(lines):
        if line[0] == '{':
            line = line[1:]
        if line[-1] == ',':
            line = line[:-1]
        line = line.replace('}','')
        key,params = line.split(':{')
        params = [p.split(':') for p in params.split(',')]
        params = [(p[0],float(p[1])) if p[0] in ['noddwell'] else (p[0],int(p[1])) for p in params]
        lines[idx] = (float(key),dict(params))

    rofdict = dict(lines)
    return rofdict
    

def plan_obsblock(obsblock,mistab,
                  legno=None,
                  swc_opt = 30*u.s,
                  lwc_opt = 90*u.s,
                  dual_opt = 40*u.s,
                  #grism_opt = 30*u.s,
                  swc_bounds = [10,35],
                  lwc_bounds = [10,65],
                  dual_bounds = [20,55],
                  niter=8,
                  los_buffer=30,
                  multiprocess=True,
                  basic=False,
                  quick=False,
                  config=False,
                  ROF=None,
                  rofrate=None,
                  label=None,
                  dcs=None,refresh_cache=False,
                  fixed=False,
                  allaorids=False,
                  **kwargs):
    '''Plan observing cadence for each obsblock'''
    if not dcs:
        # initialize DCS link
        dcs = DCS.DCS(refresh_cache=refresh_cache)

    # get globals from cfg if present
    if config:
        cfg = config
    else:
        cfg = {}


    # ensure these are in secs
    swc_opt = u.Quantity(swc_opt,u.s)
    lwc_opt = u.Quantity(lwc_opt,u.s)
    dual_opt = u.Quantity(dual_opt,u.s)
    swc_bounds = u.Quantity(swc_bounds,u.s)
    lwc_bounds = u.Quantity(lwc_bounds,u.s)
    dual_bounds = u.Quantity(dual_bounds,u.s)
    
    #grism_opt = u.Quantity(grism_opt,u.s)

    # skip blank legs
    if obsblock in (None,'','--'):
        return None

    # pop rofdict if present
    rofdict = kwargs.pop('rofdict',None)


    # get leg
    mistab = MIS.as_table(mistab)
    
    leg = mistab[mistab['ObsBlk'] == obsblock]
    if len(leg) > 1:
        leg = mistab[mistab['Leg'] == legno]

    if not leg:
        return None

    # get leg duration
    dur = [float(x) for x in leg['ObsDur'][0].split(':')]
    dur = (dur[0]*u.hour+dur[1]*u.minute+dur[2]*u.second).to(u.minute)

    # get max rofrate
    if rofrate is None:
        rofrate = get_ROF_rate(mistab,leg=leg['Leg'])
        rofrate = np.max(np.abs(rofrate))
    else:
        rofrate = rofrate if hasattr(rofrate,'unit') else rofrate*u.deg/u.min

    # get AOR
    aor = dcs.getAORs(leg['ObsBlk'])

    if aor is None or len(aor) == 0:
        # weirdly, no AORIDs associated with obsblk
        aor = d.getAOR(leg['planID'])

        # likely a calibrator, in which case flag as calibrator
        ###cal = True
        ### THIS IS UNSAFE.  WARN INSTEAD.
        raise RuntimeWarning('ObsBlk %s is empty. Generating FAORs for every AORID in ObsPlan' % leg['ObsBlk'])
    #else:
    #    cal = False
    
    print('Leg%02i'%leg['Leg'])
    print('-----')
    aor = AOR.as_table(aor)
    aor.pprint()
    print()
    
    # requested duration
    req = aor['TotalTime']*u.s
    
    
    # compute part of leg allocated to each aorid
    legshare = ((req/req.sum())*dur).to(u.minute)

    # num rewind cadence in leg per aorid
    tlos, nlos, span = compute_LOS_cadence(legshare,rofrate)
    if tlos == 0:
        tlos = dur

    # fix dither column
    dithers = aor['DitherPattern']
    aor.add_column(Column([int(x[0]) if x not in ('','None',None) else None for x in dithers],name='NumDith'))

    # iterate over each aorid in leg
    rows = deque()
    keep = deque()
    for rewinds,taor,r,treq in zip(nlos,legshare,aor,req):
        if rewinds == 0:
            rewindbounds = [0,0]
        else:
            # all rewinds -1 to +2
            rewindbounds = [rewinds-1,rewinds+2]

        # get num dithers
        dithers = int(r['NumDith']) if r['NumDith'] else 0
        ditherbounds = [dithers,dithers]

        # set default nods
        if 'LWC' in r['InstrumentConfiguration']:
            noddwell = lwc_opt
            nodbounds = lwc_bounds
        elif 'DUAL' in r['InstrumentConfiguration']:
            noddwell = dual_opt
            nodbounds = dual_bounds
        else:
            noddwell = swc_opt
            nodbounds = swc_bounds
            
        if r['ObsPlanMode'] == 'C2NC2':
            repeats = 0
            repeatbounds = [0,0]
            c2nc2 = True
            loops = 1
            loopbounds = [1,10]

            if dithers in (0,1,3,5,7,9) or dithers > 10:
                dithers = 2
            rewindbounds = [0,0]
            ditherbounds = [2,10]
            rewinds = 0 # rewinds not allowed in c2nc2

        else:
            repeats = 1
            repeatbounds = [1,35]
            c2nc2 = False
            loops = 1
            loopbounds = [1,1]

        if r['ObsPlanMode'] == 'NXCAC':
            nxcac = True
        else:
            nxcac = False
            
        if 'ACQ' in r['ObsPlanConfig']:
            # acquisition mode set to 10 seconds
            repeats = 2
            repeatbounds = [2,2]
            c2nc2 = False
            loops = 1
            loopbounds = [1,1]
            noddwell = 10*u.s
            nodbounds = [10,10]*u.s
            dithers = 0
            ditherbounds = [0,0]
            nitrun = 1
        else:
            nitrun = niter

        '''
        if cal and not 'ACQ' in r['Mode']:
            # set calibrator mode
            repeats = 1
            repeatbounds = [1,1]
            noddwell = 5*u.s
            nodbounds = [5,5]*u.s
            c2nc2 = False
            loops = 1
            loopbounds = [1,1]

            # override niter
            nitrun = 1
        else:
            nitrun = nitrun
        '''


        flightid = mistab['FlightPlan'][0]
        # override values with config file if present
        # >>
        aorid = r['aorID']

        planid = '_'.join(aorid.split('_')[0:2])
        aoridmod = r['aorID']#'%s_%i'%(planid,int(aorid.split('_')[-1])) # _1
        mission = '_'.join((flightid,aorid))
        missionmod = '_'.join((flightid,aoridmod))
        missionplan = '_'.join((flightid,planid))
        legsec = 'Leg%02d_%s'%(leg['Leg'],mission)
        legsecmod = 'Leg%02d_%s'%(leg['Leg'],missionmod)
        legsecplan = 'Leg%02d_%s'%(leg['Leg'],missionplan)
        legmission = 'Leg%02d_%s'%(leg['Leg'],flightid)
        blk = str(leg['ObsBlk'][0])

        aorsec = (legsec,legsecmod,legsecplan,legmission,mission,missionmod,missionplan,blk,aorid,aoridmod,planid,'DEFAULT')

        # chainmap checks each section for values, starting with aorid
        chmap = CMap(cfg,aorsec,**kwargs)
        
        noddwell = chmap.getquantity('noddwell',u.s,noddwell)
        repeats = chmap.getint('repeats',repeats)
        try:
            rewinds = chmap.getint('rewinds',rewinds)
            if rewinds in (-1,-2):
                rewinds = 0
        except ValueError:
            rewinds = chmap.get('rewinds',rewinds)
            if rewinds.lower() in VALIDSTR:
                rewinds = 0
            else:
                raise ValueError('rewind value %s not understood'%rewinds)
            
        loops = chmap.getint('loops',loops)
        nodbounds = chmap.getquantity('nodbounds',u.s, nodbounds)
        repeatbounds = chmap.getints('repeatbounds', repeatbounds)
        loopbounds = chmap.getints('loopbounds', loopbounds)
        rewindbounds = chmap.getints('rewindbounds', rewindbounds)
        nitrun = chmap.getint('niter',nitrun)

        dithers = chmap.getint('dithers',dithers)
        ditherbounds = chmap.getints('ditherbounds', ditherbounds)

        # if plan id and rofdict set, grab values from that based on rofrate
        if rofdict and ((rofrate.value in rofdict) or (any([np.isclose(rofrate.value,key) for key in rofdict]))):
            params = rofdict[rofrate.value]
            for k,v in params.items():
                if k not in kwargs:
                    if k == 'noddwell':
                        noddwell = u.Quantity(v,u.s)
                    elif k == 'dithers':
                        dithers = int(v)
                    elif k == 'repeats':
                        repeats = int(v)
                    elif k == 'loops':
                        loops = int(v)
                    elif k == 'rewinds' and v.lower() not in VALIDSTR:
                        rewinds = int(v)
                    else:
                        continue
                    chmap[k] = str(v)
            fixed=True
            

        # re-set dithers
        r['NumDith'] = str(dithers)
        
        # see if aorid is in aorid list
        if allaorids is False:
            aoridlist = chmap.get('aorids')
        else:
            aoridlist = None

        if aoridlist is not None:
            aoridlist = aoridlist.split(',')
            if aorid in aoridlist or aoridmod in aoridlist:
                keep.append(aoridmod)

        # fix bounds
        if 'noddwell' in chmap:
            nodbounds = u.Quantity([noddwell,noddwell])
        if 'repeats' in chmap:
            repeatbounds = [repeats,repeats]
        if 'loops' in chmap:
            loopbounds = [loops,loops]
        if 'rewinds' in chmap and chmap.get('rewinds').lower() not in VALIDSTR:
            rewindbounds = [rewinds,rewinds]
        if 'dithers' in chmap:
            ditherbounds = [dithers,dithers]

        if chmap.getboolean('totaltime',False):
            # if totaltime is true, do not divide legshare
            #taor = r['Total Exp']*u.s
            cyc = cycle(('inttime',))
            aortime = r['TotalTime']*u.s
            nitrun = 4
        else:
            ## cycle through each of the three maximization modes
            ##cyc = cycle(('dur','eff','inttime'))
            ##aortime = taor
            ##nitrun = 2
            cyc = cycle(('inttime',)) ### UPDATE: only do inttime
            aortime = r['TotalTime']*u.s
            nitrun = 4
        # <<

        # init values
        x0 = [noddwell.to(u.s).value,repeats,rewinds,dithers,loops]

        if fixed:
            nitrun = 1
        
        if basic or 'ACQ' in r['ObsPlanConfig']:
            res = [{'x':x0}]

        else:
            # define optimizer function
            optimize = partial(optimize_aor,tlos=tlos.to(u.s).value,taor=aortime.to(u.s).value,
                               c2nc2=c2nc2,nxcac=nxcac,
                               los_buffer=los_buffer)
            #res = minimize(optimize,x0,bounds=(nodbounds.to(u.s).value,repeatbounds,rewindbounds,[1,None]))

            bounds = (nodbounds.to(u.s).value,repeatbounds,rewindbounds,ditherbounds,loopbounds)

            # define execution function
            #basinfunc = partial(do_basin,func=optimize,x0=x0,bounds=bounds,fixed=fixed or cal)
            if quick:
                basinfunc = partial(do_basin,func=optimize,x0=x0,bounds=bounds,fixed=fixed,
                                    niter=100,stepsize=2,interval=10)
            else:
                basinfunc = partial(do_basin,func=optimize,x0=x0,bounds=bounds,fixed=fixed)

            # define niter seeds
            seeds = np.random.randint(0,1000,size=nitrun)
            # attach maximization functions to seeds
            seeds = list(zip(seeds,cyc))

            # deploy basinhopper niter times
            if nitrun > 1:
                res = ProgressBar.map(basinfunc,seeds,
                                      multiprocess=multiprocess)
            else:
                res = [basinfunc(seeds[0])]
                
        print()
        print(', '.join((r['aorID'],r['InstrumentConfiguration'],r['ObsPlanMode'],'%i Dithers'%dithers,'%.2f @ %.1f deg LOS'%(rofrate.value,span.to(u.deg).value))))

        results = deque()
        '''
        #for run in ProgressBar(range(0,niter)):
            #res = basinhopping(optimize, x0, niter=2500,stepsize=2,interval=20,
            #                   minimizer_kwargs={'bounds':bounds})
        '''

        for re in res:
            inputs = {k:v for k,v in zip(('noddwell','repeats','rewinds','dithers','loops'),re['x'])}
            #inputs['dithers'] = dithers
            inputs['noddwell'] = np.around(inputs['noddwell'],1) # round to 1/10th of a second
            inputs['TAOR'] = taor.to(u.s).value

            for k in ('repeats','repeats','rewinds','dithers','loops'):
                # make sure these are integers
                inputs[k] = int(np.rint(inputs[k]))

            # evaluate solution
            solution = itime3(**inputs,c2nc2=c2nc2,nxcac=nxcac)

            outputs = {k:v for k,v in zip(('DURPERREW','FDUR','INTTIME','EFF'),solution)}
            outputs['TLOS'] = tlos.to(u.s).value

            inputs.update(outputs)
            results.append(inputs)

        results = FTable(rows=list(results))
        results = results[['EFF','INTTIME','FDUR','TAOR','DURPERREW','TLOS','noddwell','repeats','rewinds','dithers','loops']]

        # add requested time
        results.add_column(Column([treq.value for _ in results],name='TREQ'),index=2)

        for col in ('EFF','INTTIME','TREQ','FDUR','DURPERREW','TAOR','TLOS'):
            # set formatter to two decimal places
            results[col] = np.around(results[col],2)
            results[col].format = '%.2f'

        #results = vstack(results)
        results = unique(results,keys=['noddwell','repeats','rewinds','dithers','loops'])

        # Sort by closeness to TAOR
        #idx = np.argsort(np.abs(results['FDUR']-results['TAOR']))
        #results = results[idx]
        #results.sort(['EFF','INTTIME'])
        #results.reverse()

        results = sort_and_eval(results,eff=True)

        # print and highlight chosen line
        results.pprint(fline=0)
        results.add_column(Column([r['aorID'].strip()]*len(results),name='AORID'),index=0)
        results.add_column(Column([span.value]*len(results),name='TLSPN'),index=results.index_column('noddwell'))

        # override rewinds if -1
        try:
            rewinds = chmap.getint('rewinds')
            if rewinds in (-1,-2):
                results['rewinds'][0] = rewinds
        except ValueError:
            rewinds = chmap.get('rewinds')
            if rewinds.lower() in VALIDSTR:
                results['rewinds'][0] = -2

        # get stop
        if 'stop' in chmap and chmap.get('stop').lower() in VALIDSTR:
            results.add_column(Column([True]*len(results),name='STOP'))
        else:
            results.add_column(Column([False]*len(results),name='STOP'))
            
        
        rows.append(results[0])


    # combine results
    table = Table(rows=rows,names=rows[0].colnames)
    for col in ('INTTIME','TREQ','FDUR','TAOR','DURPERREW','TLOS','noddwell'):
        table[col].unit = u.s
    table['EFF'].unit = '%'
    table['TLSPN'].unit = u.deg

    # remove TAOR cuz it is confusing
    table.remove_column('TAOR')

    print()
    table.pprint()
    print('---------------------\n')

    # add metadata
    table.meta['CFILENAME'] = aor['FILENAME'][0]
    #table.meta['DitherPattern'] = {r['_AORID']:int(r['Num Dith']) for r in aor}
    table.meta['DitherPattern'] = {row['AORID']:int(row['dithers']) for row in table}
    
    #comments = [', '.join((r['Mode'],r['Obs Type'],'%i Dithers'%int(r['Num Dith']),'%.2f @ %.1f deg LOS'%(rof.value,span.to(u.deg).value))) for r in aor]
    comments = [', '.join((r['InstrumentConfiguration'],r['ObsPlanMode'],'%i Dithers'%int(row['dithers']),'%.2f @ %.1f deg LOS'%(rofrate.value,span.to(u.deg).value))) for r,row in zip(aor,table)]
    # add obs details to comment
    for i,comment in enumerate(comments):
        #comments = ['%s\n%s'%(comment,str(t[
        #det = str(Table(table[i])[['EFF','INTTIME','FDUR','TAOR','DURPERREW','TLOS']])
        det = str(Table(table[i])[['EFF','INTTIME','TREQ','FDUR','DURPERREW','TLOS']])
        comments[i] = ('%s\n%s\n' % (comment,det)).replace('\n','\n#   ')

    table.add_column(Column(comments,name='Comments'))
    table.meta['Flight'] = mistab['FlightPlan'][0]
    table.meta['Leg'] = int(leg['Leg'])
    table.meta['KEEPIN'] = list(keep)

    # add leg and obs dur info
    #lkey = 'Leg%i'%leg['Leg']
    for key,nkey in zip(('Duration','ObsDur'),('Leg Dur','Obs Dur')): ### names in FAOR table
        misleg = mistab[mistab['Leg'] == leg['Leg']]
        if key in mistab.colnames:
            tstr = misleg[key][0]              # orig string
            tobj = datetime.strptime(tstr,'%H:%M:%S')  # datetime obj
            tsec = timedelta(hours=tobj.hour,minutes=tobj.minute,seconds=tobj.second)
            tsec = '%i'%int(tsec.total_seconds())
            table.meta[nkey] = '%s (%s s)' % (tstr,tsec)

    totdur = u.Quantity(np.sum(table['FDUR']),u.s).value
    mm, ss = divmod(totdur, 60)
    hh, mm = divmod(mm, 60)
    totobj = "%02d:%02d:%02d" % (hh, mm, ss)

    table.meta['Tot Dur'] = '%s (%s s)' % (totobj,'%i'%totdur)

    if ROF is not None:
        table.meta['ROF'] = '%.1f deg' % u.Quantity(ROF,u.deg).value
    if label:
        table.meta['Start'] = label

    if 'rewinds' in table.colnames:
        table.rename_column('rewinds','rewind')

    return table

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
    parser = argparse.ArgumentParser(description='Plan FORCAST observations, accounting for LOS rewind cadence.')
    parser.add_argument('flightid',type=str,help='Flight ID (e.g. 201902_FO_OTTO)')
    parser.add_argument('-niter',type=int,default=6,help='Number of iterations per AORID (default=6)')
    parser.add_argument('--basic',action='store_true',help='If specified, generate basic FAORs without optimization')
    parser.add_argument('--quick',action='store_true',help='If specified, generate FAORs with quick optimization.  May not find optimal solution.')
    #parser.add_argument('-nbounds',nargs=2,type=float,default=[1,120],help='Min Max noddwell times in sec (default=1 120)')
    parser.add_argument('-losbuffer',type=float,default=30,help='Allowed time in sec to go over LOS rewind cadence (default=30)')
    parser.add_argument('--debug',action='store_false',help='If specified, disable multiprocessing')
    parser.add_argument('-leg',type=int,default=None,help='Only process this leg')
    parser.add_argument('-interval',type=int,default=None,
                        help='Split leg into intervals of this many minutes (requires -leg)')
    parser.add_argument('--allaors',
                        action='store_true',help='If specified, ignore aorids listed in config')
    parser.add_argument('-alias',action='append',nargs=2,help='Alias obsblocks (e.g. blk1 is blk2)')
    parser.add_argument('-rofrate',type=float,default=None,help='Fix rof rate (in deg/min) for all legs')
    parser.add_argument('--dry-run',dest='dry',action='store_true',help='If specified, print plan and exit')
    parser.add_argument('-o',type=str,help='Output directory (default=flight series)')
    parser.add_argument('-v','--version',
                        dest='version',type=str,
                        default='',help='Specify Rev number/letter')
    parser.add_argument('-local',type=str,default=None,
                        help='Specify local directory for .mis files, else query DCS')
    parser.add_argument('-cfg',type=str,default=None,
                        help='Specify .cfg file for additional options')
    parser.add_argument('-mcfg',type=str,default='dcs/DBmodels.cfg',help='Model config file (default=dcs/DBmodels.cfg)')
    parser.add_argument('--comment-out','-co',
                        dest='comment',action='store_true',help='If specified, leave unplanned run blocks in FAORs rather than delete them')
    parser.add_argument('--fixed',
                        dest='fixed',action='store_true',help='If specified, do not iterate, and use initial values.')
    parser.add_argument('-noddwell',type=float,default=None,help='Specify noddwell for all legs')
    parser.add_argument('-loops',type=int,default=None,help='Specify loops for all legs')
    parser.add_argument('-repeats',type=int,default=None,help='Specify repeats for all legs')
    parser.add_argument('-dithers',type=int,default=None,help='Specify dithers for all legs')
    parser.add_argument('-r','--refresh-cache',
                        dest='refresh_cache',action='store_true',
                        help='Force update from DCS')

    args = parser.parse_args()

    # initialize DCS link and database
    mcfg = ConfigParser()
    mcfg.read(args.mcfg)

    dcs = DCS.DCS(refresh_cache=args.refresh_cache,modelcfg=mcfg)
    register_models(dcs)

    args.flightid = args.flightid.upper()

    # check for alt status
    alt = 'ALT' in args.flightid
    if alt:
        args.flightid = args.flightid.replace('_ALT','')

    # flight name/series
    split = args.flightid.split('_')
    if len(split) == 3:
        # single flight
        flightids = [args.flightid]
        if alt:
            flightids = ['%s_ALT'%flightids[0]]
        names = [split[-1]]
        if alt:
            names = ['%s_ALT'%names[0]]
        series = '_'.join(split[:-1])
    
    elif len(split) == 2:
        # series, get flightids from dcs
        flightids = dcs.getFlightSeries(args.flightid, get_ids=True, local=args.local)
        names = [f.split('_') for f in flightids]
        names = ['_'.join(name[-2:]) if 'ALT' in name else name[-1] for name in names]
        series = args.flightid
        
    else:
        raise ValueError('Flight ID "%s" not understood'%args.flightid)

    if args.o:
        odir = Path(args.o)

        if len(names) == 1:
            odirs = [odir for name in names]
        else:
            odirs = [(odir/name) for name in names]
    else:
        odir = series
        if args.version:
            odir = '%s_v%s' % (series,args.version)
        odir = Path(odir)
        odirs = [(odir/name) for name in names]

    if args.allaors:
        odirs = [odir/'allfaors' for odir in odirs]
    else:
        odirs = [odir/'faors' for odir in odirs]

    # get specified command line args
    kwargs = {k:getattr(args,k) for k in ('noddwell','repeats','loops','dithers') if getattr(args,k) is not None}

    # alias mapping
    if args.alias:
        args.alias = {a[0]:a[1] for a in args.alias}

    # read config
    if args.cfg:
        cfg = ConfigParser()
        cfg.read(args.cfg)
    else:
        cfg = {}

    # process each flightid
    for fid,odir in zip(flightids,odirs):
        mis = dcs.getFlightPlan(fid, local=args.local)

        # get utctab for dividing legs into intervals
        utctabs = dcs.getFlightPlan(fid, local=args.local, utctab=True)
        utctabs = list(filter(lambda x:'Leg' in x.meta, utctabs))

        # generate chainmap to prioritize overrides from command line
        '''
        if args.leg is None:
            keys = ['_'.join(('Leg%02d'%leg,fid)) for leg in mis['Leg'] if leg]
        else:
            keys = ['_'.join(('Leg%02d'%args.leg,fid))]
        '''

        # get cmd args as dict
        '''
        cmdargs = {k:v for k,v in vars(args).items() \
                   if ((k in ('interval','rof',
                              'noddwell','loops','repeats','dithers')) and \
                       (k in vars(args)) and \
                       (vars(args)[k] is not None))}
        '''
        #chmap = CMap(cfg,keys,**cmdargs)
        #print(chmap)
        
        '''
        # if args.leg, only look for 'interval' in cfg
        if args.leg is not None:
            if args.interval is None:
                key = '_'.join(('Leg%02d'%args.leg,fid))
                if key in cfg and 'interval' in cfg[key]:
                    args.interval = cfg[key].getint('interval')
                    getUTC = True
                else:
                    getUTC = False
            else:
                getUTC = False
        else:
            # look for legs/intervals in cfg
            keys = ['_'.join(('Leg%02d'%leg,fid)) for leg in mis['Leg'] if leg]
            if any([key in cfg for key in keys]):
                getUTC = True

        # get utctab for dividing legs into intervals
        if args.interval:
            if args.leg is None:
                raise RuntimeError('leg must be specified for interval leg splitting.')
            utctabs = d.getFlightPlan(fid, local=args.local, utctab=True)
            utctabs = list(filter(lambda x:'Leg' in x.meta, utctabs))
        '''

        # process each leg/obsblock
        for leg in mis:
            if args.leg is not None:
                if leg['Leg'] != args.leg:
                    continue

            # if no obsblk, skip
            if not leg['ObsBlk']:
                continue
            

            # apply aliases, if any
            if args.alias and leg['ObsBlk'] in args.alias:
                orig = leg['ObsBlk']
                leg['ObsBlk'] = args.alias[orig]
                leg['AOR'] = '_'.join(args.alias[orig].split('_')[1:3])

            # cfg key for leg/mission
            key = '_'.join(('Leg%02d'%leg['Leg'],fid))

            # IGNORE ALT IN KEY
            ##key = key.replace('_ALT','')
            

            # look in cfg for leg rofrate, but override with cmdline
            if args.rofrate:
                rofrate = args.rofrate
            elif ((key in cfg) and cfg[key].getfloat('rofrate',None)):
                interval = cfg[key].getfloat('rofrate')
            else:
                rofrate = None


            # process intervals, if any
            if args.interval:
                interval = args.interval
            elif ((key in cfg) and cfg[key].getint('interval',None)):
                interval = cfg[key].getint('interval')
                #if interval not in ALLOWED_INTERVALS:
                #    raise ValueError('Intervals must be in %s'%str(ALLOWED_INTERVALS))
            else:
                interval = 0

            # calculate interval start times
            utctab = list(filter(lambda x:x.meta['Leg'] == leg['Leg'], utctabs))[0]
            utctab.pprint()
            iTab = split_leg(utctab, interval, rofrate=rofrate)

            # just print out table if dry-run is set
            if args.dry:
                print('Leg%02i'%leg['Leg'])
                print('-----')
                if all([not l for l in iTab['label']]):
                    iTab.remove_column('label')
                iTab.add_column(Column([leg['ObsBlk']]*len(iTab),name='ObsBlk'),index=0)
                dur = (((iTab['stop']-iTab['start']).sec)/60)*u.min
                iTab.add_column(Column(dur,name='duration',format='%.1f'),index=3)
                iTab.pprint()
                print()
                print()
                continue

            # if rofdict in Plan ID cfg, override all values in plan_obsblock
            planID = '_'.join(leg['ObsBlk'].split('_')[1:3])
            if planID in cfg and 'rofdict' in cfg[planID]:
                # proc rofdict
                rofstr = cfg[planID]['rofdict']
                rofdict = read_rofdict(rofstr)
                kwargs['rofdict'] = rofdict

            for row in iTab:
                # generate paramtable for each interval (just one if interval == 0)
                paramtable = plan_obsblock(leg['ObsBlk'],mis,
                                           legno=leg['Leg'],
                                           niter=args.niter,
                                           multiprocess=args.debug,
                                           los_buffer=args.losbuffer,
                                           basic=args.basic,
                                           quick=args.quick,
                                           config=cfg,
                                           ROF=row['ROF'],
                                           rofrate=rofrate if rofrate else row['ROFrt'],
                                           label=row['label'],
                                           dcs=dcs,refresh_cache=args.refresh_cache,
                                           fixed=args.fixed,
                                           allaorids=args.allaors,
                                           **kwargs)
                
                if paramtable is None:
                    # invalid leg
                    continue
                
                # generate basic FAOR and restrict to AORIDs in this obsblk
                faors = dcsFAOR.AOR_to_FAORdicts(paramtable.meta['CFILENAME'],
                                                 aorids=list(paramtable['AORID']),
                                                 comment=args.comment)

                # create output directory
                odir.mkdir(parents=True,exist_ok=True)

                print('Leg%02i   %s'%(leg['Leg'],row['label']))
                print('-----')

                # update faors with optimal params
                for faor in faors:
                    faor.modify_run(paramtable)

                    ### THIS IS DANGEROUS with commment=True--empty faors wont load
                    if paramtable.meta['KEEPIN']:
                        faor.keep_aorids(paramtable.meta['KEEPIN'],comment=args.comment)
                        
                    if faor.is_empty():
                        # all run blocks removed, so don't write out
                        continue


                    if row['label']:
                        faor.root = 'Leg%02i__%s_%s'%(leg['Leg'],faor.root,row['label'])
                    else:
                        faor.root = 'Leg%02i__%s'%(leg['Leg'],faor.root)
                    outfile = (odir/faor.root).with_suffix('.faor')

                    faor.write(outfile)
                    print(outfile)

                    rows = dcs.ACTIVE_MODELS['FAOR'].to_rows(outfile, mcfg['FAOR'],faor=faor)
                    if rows:
                        dcs.ACTIVE_MODELS['FAOR'].replace_rows(dcs.db, rows)

                print()

            # pop rofdict from command line args for next leg
            kwargs.pop('rofdict',None)
        print()
        

if __name__ == '__main__':
    main()
