#! /usr/bin/env python
import argparse
from dcs import DCS
from pathlib import Path
from . import SCT as dcsSCT    
from configparser import ConfigParser
import warnings
from astropy.utils.console import ProgressBar
from collections import deque
from functools import partial

def register_models(dcs):
    """Register active DB models in DCS instance globally"""
    active_models = dcs.get_active_models()
    if active_models:
        for name,model in active_models.items():
            globals()[name] = model

    # finally, register dcs
    globals()['dcs'] = dcs

    return active_models


def process_aor(tup,ocfg,mcfg):
    aor,odir = tup
    # generate sct dicts
    scts = dcsSCT.AOR_to_SCTdicts(aor,ocfg,mcfg,mp=False)

    # create output directory
    odir.mkdir(parents=True,exist_ok=True)

    # write files
    write_func = partial(dcsSCT.write_files,odir=odir)
    files = list(map(write_func,scts))
    for f,s in zip(files,scts):
        s['FILENAME'] = Path(f[0]).resolve()
        # save timestamp
        stats = Path(f[0]).stat()
        ts = stats.st_mtime if stats.st_mtime > stats.st_ctime else stats.st_ctime
        s['TIMESTAMP'] = ts
        
    return files,scts



def _argparse():
    """Separate argparse callable for sphinx documentation"""

    # get defaults
    mcfg_DEFAULT = str(Path(DCS.__file__).parent.resolve()/'DBmodels.cfg')
    sct_DEFAULT = str(Path(__file__).parent.resolve()/'obsmaker.cfg')

    parser = argparse.ArgumentParser(description='Generate FIFI sct and map files.')
    parser.add_argument('flightid',type=str,help='Flight ID (e.g. 201902_FO_OTTO)')
    parser.add_argument('-leg',type=int,default=None,help='Only process this leg')
    parser.add_argument('-alias',action='append',nargs=2,help='Alias obsblocks (e.g. blk1 is blk2)')
    parser.add_argument('-o',type=str,help='Output directory (default=flight series)')
    parser.add_argument('-v','--version',
                        dest='version',type=str,
                        default='',help='Specify Rev number/letter')
    parser.add_argument('-local',type=str,default=None,
                        help='Specify local directory for .mis files, else query DCS')
    parser.add_argument('-ocfg',type=str,default=sct_DEFAULT,help='Obsmaker config file (default=/library/path/obsmaker.cfg)')
    parser.add_argument('-mcfg',type=str,default=mcfg_DEFAULT,help='Model config file (default=/library/path/DBmodels.cfg)')
    parser.add_argument('-r','--refresh-cache',
                        dest='refresh_cache',action='store_true',
                        help='Force update from DCS')
    return parser

def main():
    """Main CLI entry-point."""
    args = _argparse().parse_args()

    # initialize DCS link and database
    mcfg = ConfigParser()
    mcfg.read(args.mcfg)
    ocfg = ConfigParser()
    ocfg.read(args.ocfg)

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

    odirs = [odir/'sct' for odir in odirs]

    # alias mapping
    if args.alias:
        args.alias = {a[0]:a[1] for a in args.alias}

    # process each flightid
    aors = deque()
    print('Processing flights')
    for fid,odir in ProgressBar(list(zip(flightids,odirs))):
        mis = dcs.getFlightPlan(fid, local=args.local)

        # process each leg/obsblock
        for leg in mis:
            if args.leg is not None:
                if leg['Leg'] != args.leg:
                    continue

            # if no obsblk, skip
            if not leg['ObsBlkID']:
                continue

            # apply aliases, if any
            if args.alias and leg['ObsBlkID'] in args.alias:
                orig = leg['ObsBlkID']
                leg['ObsBlkID'] = args.alias[orig]
                leg['AOR'] = '_'.join(args.alias[orig].split('_')[1:3])

            # get aors
            aor = dcs.getAORs(leg['ObsBlkID'])
            if aor is None or len(aor) == 0:
                # weirdly, no AORIDs associated with obsblk
                aor = dcs.getAORs(leg['planID'])

                # likely a calibrator, in which case flag as calibrator
                ### THIS IS UNSAFE.  WARN INSTEAD.
                #raise RuntimeWarning('ObsBlk %s is empty. Generating FAORs for every AORID in ObsPlan' % leg['ObsBlkID'][0])
                warnings.warn('ObsBlk %s is empty. Generating FAORs for every AORID in ObsPlan' % leg['ObsBlkID'],category=RuntimeWarning)
                for r in aor:
                    r['ObsBlkID'] = leg['ObsBlkID']

            aors.append((aor,odir))

    # process aor, odir pairs
    proc_func = partial(process_aor,ocfg=ocfg,mcfg=mcfg)
    print('Generating sct and maps')
    res = ProgressBar.map(proc_func,aors,multiprocess=True)
    files,scts = zip(*res)
    scts = (s for sct in scts for s in sct)

    dcs.ACTIVE_MODELS['SCT'].replace_rows(dcs.db, list(scts))
    


if __name__ == '__main__':
    main()
