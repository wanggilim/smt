#! /usr/bin/env python
from dcs import DCS
import argparse
from pathlib import Path
from dossier.dossier import parse_flightid, get_leg
from astropy.utils.console import ProgressBar
from functools import partial

# get default path for model config
mcfg_DEFAULT = str(Path(DCS.__file__).parent.resolve()/'DBmodels.cfg')
plankey = 'ObsPlanID'
obsblkkey = 'ObsBlkID'

def main():
    parser = argparse.ArgumentParser(description='Refresh cache for specified flight/series/obsplan.')
    parser.add_argument('ids',nargs='+',help='ObsPlan IDs, flight plans, or series')
    parser.add_argument('-local',type=str,default=None,
                        help='Specify local directory for .mis files, else query DCS')
    parser.add_argument('-alias',action='append',nargs=2,help='Alias obsblocks (e.g. blk1 is blk2)')
    parser.add_argument('-mcfg',type=str,default=mcfg_DEFAULT,help='Model config file (default=/library/path/DBmodels.cfg)')

    args = parser.parse_args()

    d = DCS.DCS(modelcfg=args.mcfg,refresh_cache=True)

    # alias mapping for obsblks
    if args.alias:
        args.alias = {a[0]:a[1] for a in args.alias}

    for name in args.ids:
        # check if name is ObsPlan
        split = name.split('_')
        if len(split[0]) == 2:
            d.getAORs(name)

        else:
            # parse input flightid query
            flightids, names, series = parse_flightid(name, d, local=args.local)
            for fid in flightids:
                mis = d.getFlightPlan(fid,local=args.local)
                if args.alias:
                    for leg in mis:
                        # apply aliases, if any
                        if leg[obsblkkey] in args.alias:
                            newblk = args.alias[leg[obsblkkey]]
                            leg[obsblkkey] = newblk
                            leg[plankey] = '_'.join(newblk.split('_')[1:3])

                # for each leg, download AOR and proposals
                get_func = partial(get_leg, dcs=d)
                print(fid)
                print('Processing legs...')
                legs = ProgressBar.map(get_func,mis,
                                       multiprocess=False)

if __name__ == '__main__':
    main()
