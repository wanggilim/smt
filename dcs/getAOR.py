#! /usr/bin/env python
import argparse
from dcs import DCS
from shutil import copy
from pathlib import Path
from astropy.utils.console import ProgressBar

def main():
    parser = argparse.ArgumentParser(description='Download AORs locally')
    parser.add_argument('planids',nargs='+',help='Input planID(s)')
    parser.add_argument('-o',type=str,default='.',
                        help='Output directory (default=.)')
    parser.add_argument('-pdf',action='store_true',
                        help='Download proposal pdfs as well')
    parser.add_argument('-r','--refresh-cache',
                        dest='refresh_cache',action='store_true',
                        help='Force update from DCS')

    args = parser.parse_args()

    d = DCS(refresh_cache=args.refresh_cache)

    odir = Path(args.o)
    odir.mkdir(parents=True,exist_ok=True)

    for planID in ProgressBar(args.planids):
        cfile = dcs.getAORs(planID, raw=True)
        outfile = '%s.aor' % planID
        copy(cfile,odir/outfile)

        if args.pdf:
        # download PDF file
            cfile = dcs.getProposal(planID)
            outfile = '%s.pdf'%planID
            copy(cfile,odir/outfile)


if __name__ == '__main__':
    main()
