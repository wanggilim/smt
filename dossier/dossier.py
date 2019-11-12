#! /usr/bin/env python
#####################################################
##
## dossier.py
## SOFIA Dossier Tool
## Author: Michael S. Gordon (mgordon@sofia.usra.edu)
##
## version: 1.0.1
##
######################################################
from dcs import DCS
import argparse
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy.table import join,vstack,Column,Table
from astropy.table.row import Row
from collections import deque
from pathlib import Path
import subprocess
from . import TEX
from shutil import copy
from astropy.utils.console import ProgressBar
from functools import partial
from configparser import ConfigParser

def get_raw_leg(leg,dcs,aordir,propdir=None,proposal=False,plankey='ObsPlanID'):
    """Get raw aor from DCS

    Download raw AOR and (optionally) proposal pdf files from DCS.

    Args:
        leg (dict): MIS dict for flight leg.
        dcs (DCS): DCS instance.
        odir (str): str or Path to output directory.
        proposal (bool): If True, download proposal pdf.
        plankey (str): Key in leg dict for ObsPlan.

    Returns:
        leg (dict): Return input MIS leg.
    
    """
    planID = leg[plankey]
    if not planID:
        return None

    #msg = '  Retrieving AOR %s for %s' % (leg['AOR'],leg['LegName'])
    odir = Path(aordir)
    cfile = dcs.getAORs(planID,raw=True)
    outfile = '%s.aor'%planID
    copy(cfile,odir/outfile)

    if proposal:
        #msg = '\n  Retrieving proposal for AOR %s' % (leg['AOR'])
        cfile = dcs.getProposal(planID)
        outfile = '%s.pdf'%planID
        if propdir:
            copy(cfile,propdir/outfile)
        else:
            copy(cfile,odir/outfile)

    return leg
    

def get_raw_AORs(flightid, odir,
                 dcs=None, refresh_cache=False,
                 proposal=False,
                 local=None):
    """Get raw AORs from DCS by mission ID

    Download raw AORs for all legs of a flight.

    Args:
        flightid (str): FlightID e.g. 201909_HA_FEMKE
        odir (str): str or Path to output directory.
        dcs (DCS): DCS instance.  If None, a new one will be generated (possibly requiring login).
        refresh_cache (bool): If True, refresh data from DCS.  If False, use cache first.
        proposal (bool): If True, download proposal pdf.
        local (str): str or Path to local .mis file.

    Returns:
        mis (list): MIS list of dict legs.
    
    """
    aordir = Path(odir)
    aordir.mkdir(parents=True,exist_ok=True)

    if not dcs:
        # initialize DCS link
        dcs = DCS.DCS(refresh_cache=refresh_cache)
    
    # get MIS table for flightid
    mis = dcs.getFlightPlan(flightid,local=local)

    if mis is None:
        raise ValueError('No Matching MIS found for Flight Plan %s' % flightid)

    
    # for each leg, download AOR, proposals
    get_func = partial(get_raw_leg,dcs=dcs,aordir=aordir,proposal=proposal)

    print('Downloading AORs...')
    ProgressBar.map(get_func,mis,
                    multiprocess=False)
    print()
    return mis


def get_leg(leg, dcs, plankey='ObsPlanID', obsblkkey='ObsBlkID'):
    """Get AOR for input leg

    Download AOR for flight leg.

    Args:
        leg (dict): MIS dict for flight leg.
        dcs (DCS): DCS instance.
        plankey (str): Key in leg dict for ObsPlanID.
        obsblkkey (str): Key in leg dict for ObsBlkID.

    Returns:
        aor (list): AOR list of dicts.
    
    """
    if not leg[plankey]:
        return None
    
    #print('  Retrieving AOR %s for %s' % (leg['AOR'],leg['LegName']))
    #aor = dcs.getAORs(['07_0225','07_0225_8'],guide=True,pos=True,as_table=True)
    #print(aor)
    aor = dcs.getAORs(leg[obsblkkey])
    if isinstance(aor,dict):
        aor = [aor]

    elif aor is None:
        # likely a calibrator
        return leg

        
    aor = [dict(leg,**a) if a else None for a in aor]

    # Update meta data for easier retrieval of PI info later by LEG
    #aor.meta['LEG%iPI'%leg['Leg']] = aor.meta['PI']
    #print('    ObsBlk: %s'%leg['ObsBlk'])

    return aor


def generate_dossier(flightid, odir,
                     dcs=None, refresh_cache=False,
                     version='',
                     template='template.tex',
                     local=None,
                     config=None,
                     mconfig=None,
                     alias=None,
                     faor=False,
                     posfiles=False,
                     reg=False,
                     sio=False,
                     savefits=False,
                     irsurvey=None,
                     writetex=True,
                     plankey='ObsPlanID',obsblkkey='ObsBlkID'):
    """Generate dossier from flightid

    Generate dossier from template for given flightid.

    Args:
        flightid (str): Mission/flight ID (e.g. 201909_HA_FABIO).
        odir (str): Output directory.
        dcs (DCS): DCS instance.  If none, a new one will be generated (possibly requiring login).
        refresh_cache (bool): If True, refresh data from DCS.  If False, use cache first.
        version (str): Version or revision number/letter for dossier.
        template (str): Template .tex file with jinja2 render syntax.
        local (str): Specify local directory with .mis files rather than pull them from DCS.
        config (str): Config file with options for ObsBlks/Legs in dossier.
        mconfig (str): Config file with database schema.
        alias (dict): Dictionary of ObsBlks to replace or alias (e.g. for dummy FORCAST cal files).
        faor (bool): If True, process FORCAST configurations from .faor files.  These files must have been generated from companion 'planner.py' script and must exist in an 'faors/' folder within the output directory structure.
        posfiles (bool): If True, pull down .pos files from DCS.
        reg (bool): If True, generate DS9 region files of overlays.
        sio (bool): If True, overwrite ObsBlk comments with SIO instructions (HAWC+ only).
        savefits (bool): If True, copy SkyView .fits files to local directory structure.
        irsurvey (str): Survey in SkyView imaging catalogs to download additional .fits images. This can be any imaging archive known to SkyView.  Requires savefits=True.
        writetex (bool): If False, suppress .tex file generation.  Mostly for testing.
        plankey (str): Key in leg dict for ObsPlanID.
        obsblkkey (str): Key in leg dict for ObsBlkID.

    Returns:
        output (Path): Path to output .tex file.
    """

    # create output directory
    odir = Path(odir)
    odir.mkdir(parents=True,exist_ok=True)

    if not dcs:
        # initialize DCS link
        dcs = DCS.DCS(refresh_cache=refresh_cache,modelcfg=mconfig)

    if 'ALT' in flightid:
        name = '_'.join(flightid.split('_')[-2:])
    else:
        name = flightid.split('_')[-1]
    title = flightid.replace('_','\_')
    if version != '':
        title = '%s Rev %s' % (title,version)
    
    # get MIS table for flightid
    mis = dcs.getFlightPlan(flightid, local=local)

    if mis is None:
        raise ValueError('No Matching MIS found for Flight Plan %s' % flightid)

    if alias:
        for leg in mis:
            # apply aliases, if any
            if leg[obsblkkey] in alias:
                newblk = alias[leg[obsblkkey]]
                leg[obsblkkey] = newblk
                leg[plankey] = '_'.join(newblk.split('_')[1:3])
    
    #print('Retrieving MIS file for %s' % flightid)
    #mis.pprint()
    print()

    # for each leg, download AOR, proposals
    get_func = partial(get_leg,dcs=dcs)

    print(flightid)
    print('Processing legs...')
    legs = ProgressBar.map(get_func,mis,
                           multiprocess=False)

    legs = [leg for leg in legs if leg is not None]

    # download POS file
    '''
    try:
        pos = dcs.getPOS(flightid,guide=guide)
        guide = pos.guide

        # just keep target and AORID
        pos = pos[['AORID','Target']]
        pos.rename_column('Target','POSname')

        mis = join(mis,pos,join_type='left',keys=['AORID'])
    except IndexError:
        # failed to read pos file
        print('Skipping .pos for %s' % flightid)
        mis.add_column(Column(mis['Target'],name='POSname'))
        guide = None

    exit()
    mis = MIS.as_table(legs)
    print(mis)
    exit()

    mis.sort(['Leg','AORID'])
    #mis.pprint()

    # group by legs
    legs = mis.group_by('Leg')
    for leg in legs:
        leg.meta['Flight Plan ID'] = flightid

    # obsblk comments
    #legs.meta['ObsBlkComments'] = comments

    '''

    # tableify
    '''
    tables = list(map(MIS.as_table,legs))

    if guide:
        for table in tables:
            aorids = list(table['aorID'])
            guide = dcs.getPOS(aorids,guide=True,as_table=True)
            table.guide = guide
    '''
    tables = legs

    # output tex file
    output = odir/('%s.tex'%flightid)

    print()
    #print('Writing to %s...'%output)
    print('Generating %s...' % output)
    TEX.write_tex_dossier(tables, name, title, output,
                          template=template,
                          config=config,
                          mconfig=mconfig,
                          refresh_cache=refresh_cache,
                          faor=faor,
                          posfiles=posfiles,
                          dcs=dcs,
                          local=local,
                          reg=reg,
                          sio=sio,
                          savefits=savefits,
                          irsurvey=irsurvey,
                          writetex=writetex)
    return output


def _register_models(dcs):
    """Register active DB models in DCS instance globally

    Make AOR, MIS, GUIDE (etc.) database models accessible as global classes.  Used for debugging.

    Args:
        dcs (DCS): DCS instance with active models

    Returns:
        active_models (dict): Dictionary of modelname:cls for active models in the database.
    """
    active_models = dcs.get_active_models()
    if active_models:
        for name,model in active_models.items():
            globals()[name] = model

    # finally, register dcs
    globals()['dcs'] = dcs

    return active_models


def _argparse():
    """Separate argparse callable for sphinx documentation"""

    # get defaults
    mcfg_DEFAULT = str(Path(DCS.__file__).parent.resolve()/'DBmodels.cfg')
    tex_DEFAULT = Path(TEX.__file__).parent.resolve()/'template.tex'
    
    parser = argparse.ArgumentParser(description='Generate dossiers from flight ID')
    parser.add_argument('flightid',type=str,help='Flight ID. Can be a single flight (201809_HA_KEANU), or a flight series (e.g. 201809_HA).')
    parser.add_argument('-v','--version',
                        dest='version',type=str,
                        default='',help='Specify Rev number/letter')
    parser.add_argument('-o',type=str,help='Output directory (default=flight series)')
    parser.add_argument('-r','--refresh-cache',
                        dest='refresh_cache',action='store_true',
                        help='Force update from DCS')
    parser.add_argument('-c','--compile',
                        dest='compile',action='store_true',
                        help='Compile tex to pdf')
    parser.add_argument('--faor',action='store_true',help='Populate values from FAORs (FORCAST only)')
    parser.add_argument('--pos',dest='posfiles',action='store_true',help='Download .pos files')
    #parser.add_argument('--guide',action='store_true',help='Highlight guide stars in images.')
    parser.add_argument('-sio','--ioinstructions',
                        dest='sio',action='store_true',
                        help='Override ObsBlk comments with instrument operator instructions (HAWC+ only)')
    #parser.add_argument('-texcmd',type=str,default='pdflatex',
    #                    help='tex compiler in $PATH (default="pdflatex")')
    parser.add_argument('-local',type=str,default=None,
                        help='Specify local directory for .mis files, else query DCS')
    parser.add_argument('-cfg',type=str,default=None,
                        help='Specify .cfg file for additional options')
    #parser.add_argument('-scfg',type=str,default='server.cfg',help='Server config file (default=server.cfg)')
    parser.add_argument('-alias',action='append',nargs=2,help='Alias obsblocks (e.g. blk1 is blk2)')
    parser.add_argument('-s','--save-aors',
                        dest='save',action='store_true',
                        help='Save AORs and proposals to mission folder')
    parser.add_argument('-se','--save-aors-exit',
                        dest='sexit',action='store_true',
                        help='Save AORs and proposals to mission folder and exit without generating dossier')
    parser.add_argument('-reg','--save-regs',
                        dest='reg',action='store_true',
                        help='Save region files to mission folder')
    parser.add_argument('-nowrite','--no-write',
                        dest='writetex',action='store_false',
                        help='Do not [over]write .tex file.')
    parser.add_argument('-sf','--save-fits',
                        dest='savefits',action='store_true',
                        help='Save .fits files of fields.')
    parser.add_argument('-irs','--ir-survey',
                        dest='irsurvey',
                        default=None,
                        help='Save .fits files of specified IR survey. Requires --save-fits options.')
    parser.add_argument('-mcfg',type=str,default=mcfg_DEFAULT,help='Model config file (default=dcs/DBmodels.cfg)')
    parser.add_argument('-t','--template',
                        dest='template',type=str,default=tex_DEFAULT,
                        help='Template file (default=dossier/template.tex)')


    return parser
    

def main():
    """Main CLI entry-point."""
    args = _argparse().parse_args()

    # initialize DCS link and database
    mcfg = ConfigParser()
    mcfg.read(args.mcfg)

    dcs = DCS.DCS(refresh_cache=args.refresh_cache,modelcfg=mcfg)
    _register_models(dcs)

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

    if alt:
        flightids = list(filter(lambda x: '_ALT' in x,flightids))
        names = list(filter(lambda x: '_ALT' in x,names))


    if args.o:
        odir = Path(args.o)
    else:
        odir = series
        if args.version:
            odir = '%s_v%s' % (series,args.version)
        odir = Path(odir)

    odirs = [odir/name for name in names]

    # alias mapping for obsblks
    if args.alias:
        args.alias = {a[0]:a[1] for a in args.alias}

    # process each flightid
    outfiles = deque()
    for fid,odir in zip(flightids,odirs):
        dcs.getFlightPlan(fid)

        if args.save or args.sexit:
            get_raw_AORs(fid,odir/'aor',
                         dcs=dcs, refresh_cache=args.refresh_cache,
                         proposal=True,
                         local=args.local)
            if args.sexit:
                continue
        
        output = generate_dossier(fid, odir,
                                  dcs=dcs, refresh_cache=args.refresh_cache,
                                  version=args.version,
                                  template=args.template,
                                  local=args.local,
                                  config=args.cfg,
                                  mconfig=args.mcfg,
                                  alias=args.alias,
                                  faor=args.faor,
                                  posfiles=args.posfiles,
                                  reg=args.reg,
                                  sio=args.sio,
                                  irsurvey=args.irsurvey,
                                  savefits=args.savefits,
                                  writetex=args.writetex)
        outfiles.append(output)
    
    if args.compile and len(outfiles) and args.writetex:
        from pdflatex import PDFLaTeX
        pdfl = (PDFLaTeX.from_texfile(o) for o in outfiles)
        print('Compiling to pdf...')
        for p,odir in ProgressBar(list(zip(pdfl,odirs))):
            p.set_output_directory(odir)
            _,_, res = p.create_pdf(keep_pdf_file=True, cwd=odir)
            
    

    

if __name__ == '__main__':
    main()

    #dcs = DCS.DCS()
    

    #_register_models(dcs)
    #dcs._force_db_sync()
    #generate_dossier('201909_HA_FEMKE',odir='test2',dcs=dcs,guide=True)
