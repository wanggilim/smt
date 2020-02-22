#! /usr/bin/env python
from dcs import DCS
import argparse
from pathlib import Path

# get default path for model config
mcfg_DEFAULT = str(Path(DCS.__file__).parent.resolve()/'DBmodels.cfg')

def main():
    parser = argparse.ArgumentParser(description='Force refresh of smt.db file from local cache')
    parser.add_argument('-mcfg',type=str,default=mcfg_DEFAULT,help='Model config file (default=/library/path/DBmodels.cfg)')

    args = parser.parse_args()

    d = DCS.DCS(modelcfg=args.mcfg)

    d._force_db_sync()

if __name__ == '__main__':
    main()
