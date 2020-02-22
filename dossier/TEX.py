#! /usr/bin/env python
from jinja2 import Environment, Template, FileSystemLoader
import os.path
from astropy.table import Table,unique,MaskedColumn
from io import StringIO
#from aladin import Aladin
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import numpy as np
import re
import matplotlib.pyplot as plt
from astroquery.skyview import SkyView
from matplotlib.cbook.deprecation import MatplotlibDeprecationWarning
from astropy.utils.exceptions import AstropyDeprecationWarning,AstropyUserWarning,AstropyWarning
from matplotlib.lines import Line2D
from aplpy import FITSFigure
from pathlib import Path
from astropy import log
from astropy.table import Table,Column,vstack,hstack,MaskedColumn
from collections import deque,defaultdict
from astropy.utils.console import ProgressBar
from regions import LineSkyRegion
from configparser import ConfigParser
from functools import partial
#from FAOR import FAOR
from .CMap import CMap
from . import MSX
#from planner import split_leg
from requests.exceptions import ChunkedEncodingError, SSLError, ConnectionError
from regions import RectangleSkyRegion,PointSkyRegion,RegionMeta,RegionVisual,write_ds9,ds9_objects_to_string,DS9Parser
from pylatexenc.latexencode import utf8tolatex
from shapely.geometry import MultiPoint, Polygon
import matplotlib.patheffects as path_effects
import tarfile
import shutil
import warnings
warnings.filterwarnings('ignore',category=MatplotlibDeprecationWarning)
warnings.filterwarnings('ignore',category=AstropyDeprecationWarning)
warnings.filterwarnings('ignore',category=AstropyUserWarning)
warnings.filterwarnings('ignore',category=AstropyWarning)
np.warnings.filterwarnings('ignore')

DEBUG = True
MP = False if DEBUG else True

if DEBUG is False:
    log.disable_warnings_logging()
    log.setLevel('ERROR')

# formatters for astropy table columns
COL_FORMATTER = lambda x: x.replace('_','\_')
INT_FORMATTER = lambda x: '%i'%int(np.float(x)) if x not in (None,'None',-9999,'--') else ''
ZERO_INT_FORMATTER = lambda x: '%i'%int(np.float(x)) if x not in (None,'None',-9999,'--',0,'0') else ''

def INT_CONVERTER(row,cols):
    """Convert row to ints from strings or floats"""
    for col in cols:
        if col not in row:
            continue
        val = row[col]
        if isinstance(val,str):
            try:
                val = np.float(val)
            except TypeError:
                continue
        if val is None or np.isnan(val):
            row[col] = ''
        elif np.float(val).is_integer():
            row[col] = int(val)
        else:
            continue
    return row
        

ROLL_RE = re.compile('\[([\d\.]*)\,\s?([\d\.]*)\]')
COMMENT_RE = re.compile('\%\s\<COMMENT\slegname=Leg\s(?P<legnum>\d\d?).*\>\n(?P<comment>[\s\S]*?(?=\n\%\s\<COMMENT\/\>))')

#COLORS = ['#ff0000','#00ff00','#0000ff']
COLORS = ['#d62728','#1f77b4','#2ca02c','#ff7f0e','#9467bd','#17becf','#e377c2']
GCOLOR = '#FFD700'
FCOLOR = '#9467bd'

FOV = {'A_TOT':(2.8,1.7),'B_TOT':(4.2,2.7),'C_TOT':(4.2,2.7),'D_TOT':(7.4,4.6),'E_TOT':(10.0,6.3),
       'A_POL':(1.4,1.7),'B_POL':(2.1,2.7),'C_POL':(2.1,2.7),'D_POL':(3.7,4.6),'E_POL':(5.0,6.3),
       'A_C2N':(1.4,1.7),'B_C2N':(2.1,2.7),'C_C2N':(2.1,2.7),'D_C2N':(3.7,4.6),'E_C2N':(5.0,6.3),
       'FORCAST_IMG':(3.4,3.2),'FORCAST_GSM':(.04,3.18),'FIF_BLUE':(.5,.5),'FIF_RED':(1,1)}

PIXSIZE = {'A_TOT':2.34,'B_TOT':4.00,'C_TOT':4.02,'D_TOT':6.90,'E_TOT':8.62,
           'A_POL':2.34,'B_POL':4.00,'C_POL':4.02,'D_POL':6.90,'E_POL':8.62,
           'A_C2N':2.34,'B_C2N':4.00,'C_C2N':4.02,'D_C2N':6.90,'E_C2N':8.62,
           'FORCAST_IMG':0.768,'FORCAST_GSM':0.768,'FIF_BLUE':12,'FIF_RED':12}


IMGOPTIONS = {'width':0.4*u.deg, 'height':0.4*u.deg,
              'survey':'DSS2 Red',
              'vmin':None, 'vmax':None,
              'recenter':None,'roll':True,
              'invert':True,'irsurvey':None,
              'compass':True,'nofigure':False,
              'observed':False}

TARFOFFSET = {'HAWC_PLUS':3*u.deg,
              'FORCAST':40.6*u.deg,
              'FIFI-LS':0*u.deg}

INST_REPL = {'University':'Univ','Universitaet':'Univ',
             'Department':'Dept',' and':' \&',' und':' \&',
             'Institute':'Inst','Institut':'Inst',
             'Observatory':'Obs',
             'fuer ':'f.\ ','der ':'d.\ ',
             'Astrophysics':'Ast','Astrophysik':'Ast',
             'Dr. ':'','Mr. ':'','Ms. ':'','Mrs. ':'',
             'Prof ':'','Prof. ':'',
             '. ':'.\ '}

HAWC_SIO = {
    'Lissajous':r"""
    - Select boresight/pupil \\
    - Select instrument configuration \\
    - Install Low/HighAlt bias \\
    - Perform Lissajous observations at each band until SNR $>5$ \\
    """,
    'Polarimetry':r"""
    - Select boresight/pupil \\
    - Select instrument configuration. \\
    - Install Low/HighAlt bias \\
    - Perform chop-nod polarimetric observations \\
    - Ensure chop-angle and chop-throw are correct \\
    """,
    'LISPOL':r"""
    - Select boresight \\
    - Select instrument configuration \\
    - Install Low/HighAlt bias \\
    - Perform scan-pol observations using Lissajous scans \\
    - Do NOT change LOS within a set of 4 scans \\
    """
    }

STRIKETHROUGH_REPL = r'\1[-1.7ex]\n\\hline\\noalign{\\vspace{\\dimexpr 1.7ex-\\doublerulesep}}'


def get_latex_env(directory):
    """Create jinja environment to load template"""
    # modify jinja environment to work with latex files
    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string = '\BLOCK{',
        block_end_string = '}',
        variable_start_string = '\VAR{',
        variable_end_string = '}',
        comment_start_string = '\#{',
        comment_end_string = '}',
        line_statement_prefix = '%%',
        line_comment_prefix = '%#',
        trim_blocks = True,
        autoescape = False,
        #loader = FileSystemLoader(os.path.abspath('.'))
        loader = FileSystemLoader(str(Path(directory).resolve()))
    )
    return latex_jinja_env

def get_latex_template(filename):
    """Return jinja template"""
    filename = Path(filename)
    env = get_latex_env(filename.parent.resolve())
    return env.get_template(filename.name)


def make_box(center, width, height, angle=0*u.deg, TARFoffset=0*u.deg,label=None,
             linewidth=2, color='#ff0000',name=None,split=False,**kwargs):
    '''Generate box overlay'''
    if split:
        # split is for the HAWC+ R0/R1 gap
        try:
            splitgap = PIXSIZE[label]*u.arcsec
        except KeyError:
            if 'G' in label:
                splitgap = PIXSIZE['FORCAST_GSM']*u.arcsec
            else:
                splitgap = PIXSIZE['FORCAST_IMG']*u.arcsec
        r_width = (width-splitgap)/2
        boxangle = (np.arctan2(height,width).to(u.deg)+90*u.deg)/2
        offset = angle + TARFoffset
        
        r0_center = center
        r1_center = center.directional_offset_by(-90*u.deg+offset,r_width+splitgap)

        recenter = center.directional_offset_by(-90*u.deg+offset,r_width/2+splitgap/2)
        
        r0 = make_box(r0_center,r_width,height,angle,TARFoffset,label,
                      linewidth,color,name,reglabel='_R0')
        r1 = make_box(r1_center,r_width,height,angle,TARFoffset,label,
                      linewidth,color,name,reglabel='_R1')

        regs = [r0.get('reg'),r1.get('reg')]

        boxdict = {'box':[r0['box'][0],r1['box'][0]],'linewidth':linewidth,
                   'color':color,'center':center,'label':label,'name':name,
                   'recenter':recenter,'reg':regs}

    elif kwargs.get('scan'):
        try:
            splitgap = PIXSIZE[label]*u.arcsec
        except KeyError:
            if 'G' in label:
                splitgap = PIXSIZE['FORCAST_GSM']*u.arcsec
            else:
                splitgap = PIXSIZE['FORCAST_IMG']*u.arcsec
        r_width = width/2
        boxangle = (np.arctan2(height,width).to(u.deg)+90*u.deg)/2
        offset = angle + TARFoffset
        # keep boresite at R0 center for HAWC+
        if 'TOT' in label:
            center = center.directional_offset_by(-90*u.deg+offset,r_width/2-splitgap*2)
            
        return make_box(center,width,height,angle,TARFoffset,label,
                        linewidth,color,name,split=False,scan=False,reglabel='_scan')

    else:
        diag = np.hypot(width/2,height/2)
        boxangle = np.arctan2(height,width).to(u.deg)+90*u.deg

        offset = angle + TARFoffset

        tl = center.directional_offset_by( boxangle+offset,diag)
        tr = center.directional_offset_by(-boxangle+offset,diag)
        bl = center.directional_offset_by(-(180*u.deg+boxangle)+offset,diag)
        br = center.directional_offset_by( 180*u.deg+boxangle+offset,diag)

        box = (tl,tr,br,bl,tl)
        box = np.array([[coord.ra.value,coord.dec.value] for coord in box])
        boxdict = {'box':[box],'linewidth':linewidth,'color':color,
                   'center':center,'label':label,'name':name}

        # make region
        if name is None:
            boxdict['reg'] = None
        else:
            # split name---last two entries are ra/dec
            rlabel = ' '.join((' '.join(name.split()[0:-2]),label,kwargs.get('reglabel','')))
            meta = RegionMeta({'label':rlabel})
            vmeta = {'color':color}
            if ('_scan' in rlabel) or ('_d' in rlabel):
                vmeta['dash'] = 1
            vmeta = RegionVisual(vmeta)
            reg = RectangleSkyRegion(center,width,height,angle=offset,
                                     meta=meta,visual=vmeta)
            boxdict['reg'] = ds9_objects_to_string([reg])
    return boxdict


def make_dithers(center,scale,angle=0*u.deg):
    """Generate box for dithers"""
    diag = np.hypot(scale,scale)
    posangle = angle+45*u.deg
    tl = center.directional_offset_by(posangle,diag)
    tr = center.directional_offset_by(posangle-90*u.deg,diag)
    bl = center.directional_offset_by(posangle-90*u.deg,-diag)
    br = center.directional_offset_by(posangle,-diag)
    return (tl,tr,bl,br)

def make_NMC(center,chopthrow=300*u.arcsec,chopangle=0*u.arcsec,label=None):
    '''Given a center and chop nod parameters, calculate the nod throws'''
    if label == 'FIF_RED':
        pass
    
    nodA = center.directional_offset_by(chopangle,chopthrow)
    nodB = center.directional_offset_by(180*u.deg+chopangle,chopthrow)

    return (nodA,nodB)

def make_C2NC2(center,chopthrow,chopangle,nodthrow,nodangle):
    '''Generate asymmetric chopping centers'''
    nodAchopB = center.directional_offset_by(chopangle,chopthrow)
    nodBchopA = center.directional_offset_by(nodangle,nodthrow)
    nodBchopB = nodBchopA.directional_offset_by(chopangle,chopthrow)
    return (nodAchopB,nodBchopA,nodBchopB)


def get_cfgoptions(cfg, blk):
    '''Get obsblk options from cfg'''
    options = {}
    for k,v in IMGOPTIONS.items():
        o = cfg.get(blk, k, fallback=v)
        try:
            o = float(o)
        except (ValueError,TypeError):
            pass
        options[k] = o

    options['height'] = u.Quantity(options['height'], u.deg)
    options['width']  = u.Quantity(options['width'],  u.deg)
    options = {k:v for k,v in options.items() if v is not None}
    if cfg.has_option(blk,'nofigure') and cfg.getboolean(blk,'nofigure'):
        options['nofigure'] = True
    if cfg.has_option(blk,'roll'):
        roll = cfg.get(blk,'roll')
        if roll.lower() in ('true','on','yes'):
            options['roll'] = cfg.getboolean(blk,'roll')
        else:
            options['roll'] = float(roll)

    return options

def get_img_options(cmap):
    '''Get img options from cmap'''
    options = {}
    for k,v in IMGOPTIONS.items():
        o = cmap.get(blk, k, fallback=v)
        try:
            o = float(o)
        except (ValueError,TypeError):
            pass
        options[k] = o

    options['height'] = u.Quantity(options['height'], u.deg)
    options['width']  = u.Quantity(options['width'],  u.deg)
    options = {k:v for k,v in options.items() if v is not None}
    if cfg.has_option(blk,'nofigure') and cfg.getboolean(blk,'nofigure'):
        options['nofigure'] = True
    if cfg.has_option(blk,'roll'):
        roll = cfg.get(blk,'roll')
        if roll.lower() in ('true','on','yes'):
            options['roll'] = cfg.getboolean(blk,'roll')
        else:
            options['roll'] = float(roll)

    return options

def get_image(overlays,survey='DSS2 Red',width=0.2*u.deg,height=0.2*u.deg,
              reticle=False,reticle_style_kwargs=None,compass=True,
              vmin=None,vmax=None,recenter=None,invert=True,fpi=False,
              irsurvey=None,**kwargs):
    '''Get image from skyview'''
    
    if overlays is None or not overlays:
        return None

    # flatten overlays---FIFI has two
    if recenter:
        if isinstance(recenter,SkyCoord):
            center = recenter
        else:
            center = SkyCoord(recenter,unit=(u.hourangle,u.deg))

    else:
        center = overlays[0]['center']

    try:
        im = SkyView.get_images(center,survey=survey,
                                width=width,height=height,
                                show_progress=DEBUG)
    except (SSLError,ChunkedEncodingError,ConnectionError):
        warnings.warn('Cannot query SkyView service. Skipping image.')
        print('Cannot query SkyView service. Skipping image.')
        return None

    try:
        hdu = im[0][0]
    except (IndexError,TypeError):
        warnings.warn('Cannot process SkyView response. Skipping image.')
        print('Cannot process SkyView response. Skipping image.')
        return None
        
    fig = FITSFigure(hdu)
    fig.show_grayscale(vmin=vmin,vmax=vmax,invert=invert)

    #regs = deque()

    for idx,overlay in enumerate(overlays):
        fig.show_polygons(overlay['box'],edgecolor=overlay['color'],lw=overlay['linewidth'])
        fig.show_markers(overlay['center'].ra.value,overlay['center'].dec.value,marker='*',edgecolor=overlay['color'])

        if 'overlay2' in overlay:
            fig.show_polygons(overlay['overlay2']['box'],edgecolor=overlay['color'],lw=overlay['linewidth'])

        '''
        if overlay.get('reg'):
            rs = [r for r in overlay['reg'] if r]
            regs.extend(rs)
        '''

        if overlay['label']:
            if 'overlay2' in overlay:
                #change label to aorid for FIFI
                overlay['label'] = overlay['aorID']
            if len(overlay['label']) < 6:
                # note, ignore this used to be (0.87, 0.95)
                fig.add_label(0.75, 0.95, '%s%s'%('\n\n'*idx,overlay['label']),
                              horizontalalignment='left', weight='bold',
                              relative=True, size='large',color=overlay['color'])
            else:
                fig.add_label(0.75, 0.95, '%s%s'%('\n\n'*idx,overlay['label']),
                              horizontalalignment='left', weight='bold',
                              relative=True, size='large',color=overlay['color'])


        if overlay['name']:
            fig.add_label(0.02, 0.95, '%s%s'%('\n\n'*idx,overlay['name']),
                          horizontalalignment='left', weight='bold',
                          relative=True, size='large',color=overlay['color'])

        if fpi:
            fpiradius = (4.5*u.arcmin).to(u.deg).value
            fig.show_circles(overlay['center'].ra.value,
                             overlay['center'].dec.value,
                             fpiradius,
                             #edgecolor=FCOLOR,
                             edgecolor=overlay['color'],
                             linestyle='dashed',lw=1,
                             alpha=0.5)
            if 'nods' in overlay and len(overlay['nods']) > 2:
                # c2nc2 mode
                fig.show_circles(overlay['nods'][1]['center'].ra.value,
                                 overlay['nods'][1]['center'].dec.value,
                                 fpiradius,edgecolor=overlay['color'],
                                 linestyle='dashed',lw=1,
                                 alpha=0.4)


    if reticle:
        pixel_width = hdu.data.shape[0]
        inner,outer = 0.03,0.08
        if reticle_style_kwargs is None:
            reticle_style_kwargs = {}
        reticle_style_kwargs.setdefault('linewidth', 2)
        reticle_style_kwargs.setdefault('color', 'm')
        ax = fig.ax

        ax.axvline(x=0.5*pixel_width, ymin=0.5+inner, ymax=0.5+outer,
                   **reticle_style_kwargs)
        ax.axvline(x=0.5*pixel_width, ymin=0.5-inner, ymax=0.5-outer,
                   **reticle_style_kwargs)
        ax.axhline(y=0.5*pixel_width, xmin=0.5+inner, xmax=0.5+outer,
                   **reticle_style_kwargs)
        ax.axhline(y=0.5*pixel_width, xmin=0.5-inner, xmax=0.5-outer,
                   **reticle_style_kwargs)

    if compass:
        ax = fig.ax
        x,y = 0.95, 0.05
        dispTrans = ax.transData.inverted()
        dispOrigin = dispTrans.transform(ax.transAxes.transform((x,y)))
        origin = fig.pixel2world(*dispOrigin)
        w = WCS(hdu.header)

        coo_origin = SkyCoord(ra=origin[0],dec=origin[1],unit=(u.deg,u.deg))
        delta = hdu.header['CDELT2']*u.deg*20
                
        coo_e = coo_origin.directional_offset_by(90*u.deg, delta)
        coo_n = coo_origin.directional_offset_by( 0*u.deg, delta)

        for c in (coo_e,coo_n):
            line_sky = LineSkyRegion(start=coo_origin,end=c)
            line_pix = line_sky.to_pixel(w)
            line_pix.meta['line'] = 1
            line_pix.visual['line'] = 1
            line_pix.visual['linewidth'] = 2
            line_pix.plot(ax=ax)
        ax.text(0.881,0.056,'E',color='g',transform=ax.transAxes,weight='bold')
        ax.text(0.942,0.122,'N',color='g',transform=ax.transAxes,weight='bold')

        roll_s = coo_origin.directional_offset_by(overlay['roll'][0],delta*0.75)
        roll_e = coo_origin.directional_offset_by(overlay['roll'][1],delta*0.75)

        for c,color in zip((roll_s,roll_e),('#1f77b4','#d62728')):
            line_sky = LineSkyRegion(start=coo_origin,end=c)
            line_pix = line_sky.to_pixel(w)
            line_pix.meta['line'] = 1
            line_pix.visual['line'] = 1
            line_pix.visual['linewidth'] = 2
            line_pix.visual['color'] = color
            line_pix.plot(ax=ax)


    #regs = list(filter(lambda x: x is not None,regs))
    # remove duplicates
    #regs = set(regs)   # NOTE: this effectively does nothing if aorids are given to make_box, since each line will be unique
    
    # convert back into objects--this is all because regions don't pickle
    #regs = [DS9Parser(reg).shapes.to_regions()[0] for reg in regs]

    # add ir image
    if irsurvey is not None:
        if 'msx' in irsurvey.lower():
            band = irsurvey.split()[-1]
            dfile = MSX.query_region(center,band=band,show_progress=True)
            if dfile is not None:
                irhdu = fits.open(dfile)
                irhdu[0].header['SURVEY'] = 'MSX Band %s'%band
                hdu = fits.HDUList([hdu,irhdu[0]])
        else:
            try:
                im = SkyView.get_images(center,survey=irsurvey,
                                        width=width,height=height,
                                        show_progress=DEBUG)
                irhdu = im[0][0]
                hdu = fits.HDUList([hdu,irhdu])
            except (SSLError,ChunkedEncodingError,ConnectionError,IndexError,TypeError):
                pass
    
    return fig, hdu
    

def make_overview(leg, tex=True):
    '''Make table of overview stats for leg'''

    # overview cols to extract
    ocols = ('Start','ObsDur','Target','ObsBlkID','Priority','RA','DEC')

    # metacols
    hcols = ('Leg','Name','PI')
    mcols = ('Elev','ROF','ROFRT','MoonAngle','MoonIllum','THdg','THdgRT')

    l = leg[0]
    overview = {k:l.get(k,'') for k in ocols}

    # make metadata
    if tex:
        if 'PI' not in l:
            l['PI'] = ''
        overview['header'] = '\\captionline{Leg %i (%s)}{%s}' % (l['Leg'],l['Name'].replace('_','\_'),l['PI'])
        # shorten header
        for k,v in INST_REPL.items():
            overview['header'] = overview['header'].replace(k,v)
    else:
        overview['header'] ={k:l.get(k,'') for k in hcols}

    # footer holds mis file info
    footer = {}
    for k in mcols:
        if '%s_start'%k in l:
            if 'RT' in k:
                footer[k] = '[%+.2f, %+.2f]' % (l['%s_start'%k],l['%s_end'%k])
            else:
                footer[k] = '[%.1f, %.1f]' % (l['%s_start'%k],l['%s_end'%k])
        elif k == 'MoonAngle':
            footer[k] = '%i$^{\circ}$'%int(l[k])
        else:
            footer[k] = l[k]

    if tex:
        footer1 = '\quad '.join(['%s: %s'%(k,footer[k]) for k in mcols[0:3]])
        footer1 = ' '.join((footer1,'deg/min'))
        # add priority
        prior = overview.pop('Priority')
        if prior:
            footer1 = '\quad '.join((footer1,'Priority: %s'%prior))
        footer2 = '\quad '.join(['%s: %s'%(k,footer[k]) for k in mcols[-4:]])
        footer2 = ' '.join((footer2,'deg/min'))
        footer = '%s\\\\\n%s'%(footer1,footer2)
        footer = footer.replace('ROFRT','rate') \
                       .replace('THdgRT','rate') \
                       .replace('Moon','Moon ') \
                       .replace('%','\%')
        footer = ''.join((r'\\[0.5em]','\n',footer))

    overview['footer'] = footer

    if tex:
        overview = generate_overview_tex(overview)
    
    return overview


def generate_overview_tex(overview, metakeys=('header','footer')):
    '''Generate latex string of overview table'''

    # col align param must be special for boldface header line
    col_align = ['c']*len(overview)
    col_align = '|^'.join(col_align)
    col_align = '|$%s|'%col_align

    if len(overview['Target']) > 15:
        # target is too long, so make cells smaller
        preamble = r'\setlength{\tabcolsep}{0.5em}'
        overview['footer'] += '\n'+r'\setlength{\tabcolsep}{1em}'
    else:
        preamble = ''

    # make meta dict
    meta = {mkey:overview.pop(mkey) for mkey in metakeys}

    # convert to table
    overview = Table(data=[overview],names=overview.keys(),meta=meta)

    # remove ra/dec cols if non-sidereal object
    if overview['RA'][0] is None:
        overview.remove_columns(('RA','DEC'))
        #overview['RA'][0] = ''
        #overview['DEC'][0] = ''

    else:
        # reformat
        coord = SkyCoord(ra=overview['RA'][0],dec=overview['DEC'][0],unit=(u.hourangle,u.deg))
        ra,dec = coord.to_string('hmsdms',sep=':',precision=2).split()
        overview['RA'][0] = ra
        overview['DEC'][0] = dec

    # safe convert obsblkid,target
    for col in ('ObsBlkID','Target'):
        try:
            overview[col] = [k.replace('_','\_') for k in overview[col]]
        except AttributeError:
            overview[col] == ''

        
    
    # rename cols to have headercolor
    for col in overview.colnames:
        newcol = '\\cellcolor{headercolor}%s' % col
        overview.rename_column(col,newcol)

    with StringIO() as f:
        
        overview.write(f,format='latex',
                       latexdict={'header_start':r'\hline\rowstyle{\bfseries}',
                                  'tablealign':'h!',
                                  'caption':overview.meta['header'],
                                  'col_align':col_align,
                                  'data_end':'\hline',
                                  'preamble':preamble,
                                  'tablefoot':overview.meta['footer']})
        texcode = f.getvalue()
        if 'RA' not in str(overview.colnames):
            # push right
            texcode = texcode.replace(r'\begin{tabular}','\\hspace*{2cm}\n\\begin{tabular}')
    return texcode

def make_details(tab, tex=True, faor=False):
    '''Make observation details'''

    if tab[0]['aorID'] in ('99_9999_99','--'):
        return ''

    instrument = tab[0]['InstrumentName']
    if instrument == 'HAWC_PLUS':
        #keys = ('ObsPlanConfig','aorID','Name','InstrumentSpectralElement1','Repeat','NodTime','ChopThrow','ChopAngle','ScanTime','ScanAmplitudeEL','ScanAmplitudeXEL','ScanRate','ChopAngleCoordinate')
        keys = ('ObsPlanConfig','aorID','target','InstrumentSpectralElement1','Repeat','NodTime','ChopThrow','ChopAngle','ScanTime','ScanAmplitudeEL','ScanAmplitudeXEL','ScanRate','ChopAngleCoordinate')
        key_map = {'ObsPlanConfig':'Mode','aorID':'AORID','ChopAngleCoordinate':'Sys','InstrumentSpectralElement1':'Band/Bore','ScanAmplitudeEL':'ScanAmp','target':'Name'}

        # replace some values
        for t in tab:
            # change coordsys
            sys = t['ChopAngleCoordinate']
            t['ChopAngleCoordinate'] = 'ERF' if sys == 'Sky' else 'SIRF'

            # store filter for boresite
            filt = t['InstrumentSpectralElement1'][-1]
            spec2 = t['InstrumentSpectralElement2']

            # scan mode
            if t['ObsPlanMode'] == 'OTFMAP':
                # combine scan amps
                el,xel = t['ScanAmplitudeEL'], t['ScanAmplitudeXEL']
                if el == xel:
                    t['ScanAmplitudeEL'] = str(el)
                else:
                    t['ScanAmplitudeEL'] = '%s/%s'%(el,xel)

                # change scan modes
                if t['ObsPlanConfig'] == 'TOTAL_INTENSITY':
                    if t['ScanType'] == 'Box':
                        t['ObsPlanConfig'] = 'BOX'
                    elif t['ScanType'] == 'Lissajous':
                        t['ObsPlanConfig'] = 'LIS'
                    else:
                        t['ObsPlanConfig'] = '?'

                    t['InstrumentSpectralElement1'] = '/'.join((filt,'Open'))
                        
                elif t['ObsPlanConfig'] == 'POLARIZATION':
                    t['ObsPlanConfig'] = 'LISPOL'
                    t['InstrumentSpectralElement1'] = '/'.join((filt,spec2[-1]))
                else:
                    t['ObsPlanConfig'] = '?'

                # set chops to none
                t['ChopAngle'] = -9999
                t['ChopThrow'] = -9999
                t['NodAngle'] = -9999
                t['NodThrow'] = -9999
                t['NodTime'] = None

            # C2N
            else:
                if t['ObsPlanConfig'] == 'POLARIZATION':
                    t['ObsPlanConfig'] = 'POL'
                    t['InstrumentSpectralElement1'] = '/'.join((filt,spec2[-1]))
                else:
                    # C2N total_intensity
                    t['ObsPlanConfig'] = 'C2N'
                    spec2 = 'Open' if spec2 == 'OPEN' else spec2[-1]
                    t['InstrumentSpectralElement1'] = '/'.join((filt,spec2))
                    
                
        # keep certain keys
        detail = [{key:t[key] for key in keys} for t in tab]

        # make table
        detail = Table(detail,names=keys)

        # rename columns
        detail.rename_columns(tuple(key_map.keys()),tuple(key_map.values()))

        # remove extra scanamp col
        detail.remove_column('ScanAmplitudeXEL')

        # if all modes are scanning, drop chop/nod params
        if all((mode in ('LIS','LISPOL','BOX') for mode in detail['Mode'])):
            detail.remove_columns(('NodTime','ChopThrow','ChopAngle','Sys'))
        # if all modes are pol, drop scan params
        if all(mode in ('POL','C2N') for mode in detail['Mode']):
            detail.remove_columns(('ScanTime','ScanAmp','ScanRate'))

        # if any dithering, make dither footer
        unit_map = {'Sky':'arcsec','Array':'pix'}
        if any((t.get('DitherPattern') for t in tab)):
            dithscale = (t.get('DitherScale') for t in tab)
            dithunit = (unit_map.get(t.get('DitherCoord')) for t in tab)
            dithband = (t['InstrumentSpectralElement1'][0] for t in tab)
            dithscale = (str(int(scale)).rjust(2).replace(' ','~') if scale else '' for scale in dithscale)
            footer = ['\t%s: %s %s' % (band,scale,unit) for band,scale,unit \
                      in zip(dithband,dithscale,dithunit) if scale]
            #footer = set(footer)  # REMOVES ORDER
            footer = list(dict.fromkeys(footer).keys())  # basically a set operation that preserves order
            if len(footer) == 1:
                footer = footer.pop()
            else:
                #footer = '\\\\\n'.join(sorted(footer))
                footer = '\\\\\n'.join(footer)
            footer = 'dither\quad\quad %s\\\\'%footer
            detail.meta['footer'] = footer

    

    elif instrument == 'FORCAST':
        keys = ['ObsPlanConfig','ObsPlanMode','aorID','Name',
                'InstrumentSpectralElement1','InstrumentSpectralElement2',
                'Repeat','NodTime',
                'ChopThrow','ChopAngle','ChopAngleCoordinate',
                'NodThrow','NodAngle','TotalTime']
        key_map = {'ObsPlanConfig':'Mode','aorID':'AORID','ObsPlanMode':'Type','ChopAngleCoordinate':'Sys','InstrumentSpectralElement1':'SWC','InstrumentSpectralElement2':'LWC'}
        faor_keys = ['Nod','Loop','Dithers','Scale','FDUR','TREW','TLOS','TLSPN','DitherCoord',
                     'Rewind','IntTime']
        mode_map = {'ACQUISITION':'ACQ','GRISM':'GSM','IMAGING':'IMG'}

        for t in tab:
            sys = t['ChopAngleCoordinate']
            t['ChopAngleCoordinate'] = 'ERF' if sys == 'Sky' else 'SIRF'

            # shorten filter config
            t['InstrumentSpectralElement1'] = t['InstrumentSpectralElement1'].replace('FOR_','').replace('OPEN','')
            t['InstrumentSpectralElement2'] = t['InstrumentSpectralElement2'].replace('FOR_','')

            # combine filters if dual mode
            if 'DUAL' in t.get('ObsPlanConfig',''):
                t['InstrumentSpectralElement1'] = '/'.join((t['InstrumentSpectralElement1'],
                                                            t['InstrumentSpectralElement2']))
                t['InstrumentSpectralElement2'] = ''

            # drop second element if OPEN
            elif t['InstrumentSpectralElement2'] == 'OPEN':
                t['InstrumentSpectralElement2'] = ''

            # shorten chop mode
            if t['NodType'] == 'Nod_Match_Chop':
                t['ObsPlanMode'] = 'NMC'

            # shorten config
            t['ObsPlanConfig'] = mode_map.get(t['ObsPlanConfig'],t['ObsPlanConfig'])


        # keep certain keys
        if 'FAORfile' in tab[0]:
            keys += faor_keys
        detail = [{key:t[key] for key in keys} for t in tab]

        # make table
        detail = Table(detail,names=keys)
        
        # rename columns
        detail.rename_columns(tuple(key_map.keys()),tuple(key_map.values()))

        # if 'Nod' present from FAOR, then replace NodTime
        if 'Nod' in detail.colnames:
            detail['NodTime'] = detail['Nod']
            detail.remove_column('Nod')

        # if all modes are NMC, drop nod/c2nc2 params
        if all((mode == 'NMC' for mode in detail['Type'])):
            detail.remove_columns(('NodAngle','NodThrow'))
            if 'Loop' in detail.colnames:
                detail.remove_column('Loop')

        # if all modes are c2nc2, drop repeats params
        if all((mode == 'C2NC2' for mode in detail['Type'])):
            if 'Repeat' in detail.colnames:
                detail.remove_column('Repeat')

        # if there are no dithers, remove dither cols
        try:
            if not any(detail['Dithers']):
                detail.remove_columns(('Dithers','Scale'))
        except KeyError:
            pass
            

        # remove repeats, rewinds, or loops if all are None
        for col in ('Repeat','Rewind','Loop'):
            if col in detail.colnames and (not any(detail[col])):
                detail.remove_column(col)

        # add TLOS info to footer
        if any([mode == 'C2NC2' for mode in detail['Type']]):
            # leave TLOS and TLSPN in there
            detail.meta['footer'] = ''
        else:
            try:
                tl = detail['TLOS'][0]
                span = detail['TLSPN'][0].split()[0]
                if tl == 'inf' or np.isinf(tl):
                    tl = '--'
                losdet = 'TLOS = %s s @ %s deg'%(tl,span)
                detail.meta['footer'] = '%s\t\\hspace{2in}' % losdet
                detail.remove_columns(('TLSPN','TLOS'))
            except KeyError:
                pass

        if 'TLSPN' in detail.colnames:
            detail.replace_column('TLSPN',
                                  Column([float(x.split()[0]) for x in detail['TLSPN']],name='TLSPN'))
            

    elif instrument == 'FIFI-LS':
        keys = ('PrimeArray','aorID','Name','TimePerPoint','Repeat','ChopType','ChopThrow','ChopAngle',
                'ChopAngleCoordinate','MapRotationAngle','TotalTime')
        sysmap = {'J2000':'ERF','HORIZON':'SIRF'}
        key_map = {'MapRotationAngle':'FOVAngle','aorID':'AORID','ChopAngleCoordinate':'Sys','PrimeArray':'Prime',
                   'TimePerPoint':'NodTime'}

        for t in tab:
            # change coordsys
            sys = t['ChopAngleCoordinate']
            t['ChopAngleCoordinate'] = sysmap.get(sys,sys)
            t['TotalTime'] = t['TimePerPoint'] * t['Repeat']
            
        # keep certain keys
        detail = [{key:t[key] for key in keys} for t in tab]

        # make table
        detail = Table(detail,names=keys)

        # rename columns
        detail.rename_columns(tuple(key_map.keys()),tuple(key_map.values()))
        
    else:
        #raise NotImplementedError('Instrument %s not implemented. %s' % (instrument, tab[0]['ObsBlkID']))
        warnings.warn('WARNING: Instrument %s not implemented. %s' % (instrument, tab[0]['ObsBlkID']))
        print('WARNING: Instrument %s not implemented. %s' % (instrument, tab[0]['ObsBlkID']))
        return ''

    if tex:
        # set formatter to replace '_' with '\_'
        for col in ('AORID','Name'):
            detail[col].format = COL_FORMATTER

        # set int formatter
        blank_zero_cols = ('ScanDur','ChopThrow','ScanTime','NodTime',
                           'ScanAmp','ScanRate','NodThrow')   # make a zero blank
        for col in ('NodTime','Repeat','ScanDur','ChopThrow','ChopAngle','ScanTime',
                    'ScanAmp','ScanRate','NodThrow','NodAngle','TotalTime','IntTime',
                    'Rewind','Loop','Dithers','FDUR','TREW','TLOS','Scale'):
            try:
                try:
                    floats = detail[col].filled(0)
                except AttributeError:
                    floats = MaskedColumn(detail[col]).filled(0)
                for i,f in enumerate(floats):
                    try:
                        if f in (None,'None','NONE','--') or np.isnan(f):
                            floats[i] = 0
                    except TypeError:
                        continue
                floats = (np.float(x).is_integer() for x in floats)
                if all(floats):
                    if col in blank_zero_cols:
                        detail[col].format = ZERO_INT_FORMATTER
                    else:
                        detail[col].format = INT_FORMATTER
            except (KeyError,ValueError):
                continue

        # force some cols to round
        for col in ('IntTime','FDUR','TLOS'):
            if col in detail.colnames:
                detail.replace_column(col,Column(np.rint(detail[col]),name=col))
                detail[col].format = INT_FORMATTER

        # set units
        detail.meta['units'] = {'NodTime':'s','ChopThrow':r'$^{\prime\prime}$','ChopAngle':r'$^\circ$','ScanDur':'s','ScanAmp':r'$^{\prime\prime}$','ScanRate':'$^{\prime\prime}$/s','TotalTime':'s','NodDwell':'s','NodAngle':r'$^\circ$','NodThrow':r'$^{\prime\prime}$','IntTime':'s','FDUR':'s','TLOS':'s','TREW':'s','TLSPN':r'$^\circ$'}

        caption = '\\captionline{Observation details %s:}{}' % tab[0]['ObsBlkID']
        detail.meta['caption'] = caption.replace('_','\_')

        # add strikethrough
        observed_dict = {t['aorID']:t.get('observed',False) for t in tab}
        detail.meta['observed'] = observed_dict


        # if FORCAST, split into two tables
        if instrument == 'FORCAST' and 'FAORfile' in tab[0]:
            detail2 = detail.copy()

            # fix metadata
            if 'footer' in detail.meta:
                del detail.meta['footer']
            if 'caption' in detail2.meta:
                del detail2.meta['caption']

            d_keep = filter(lambda x: x in detail.colnames, ('Mode','Type','AORID','Name','SWC','LWC',
                                                             'ChopThrow','ChopAngle','NodThrow','NodAngle',
                                                             'Sys','TotalTime'))
            d2_keep = filter(lambda x: x in detail2.colnames, ('AORID','Repeat','NodTime','Dithers','Scale','Loop',
                                                               'FDUR','TREW','TLOS','TLSPN','IntTime'))

            detail.keep_columns(list(d_keep))
            detail.rename_column('TotalTime','ReqTime')
            detail.meta['units']['ReqTime'] = 's'

            detail2.keep_columns(list(d2_keep))
            if 'Scale' in detail2.colnames:
                detail2.meta['units']['Scale'] = r'$^{\prime\prime}$'

            detail = [detail,detail2]
            
        # return tex string
        detail = generate_details_tex(detail)

    else:
        # convert to list of dicts
        detail = detail.to_pandas().to_dict('records')
                
        # set int formatter
        intcols = ('NodTime','Repeat','ScanDur','ChopThrow','ChopAngle','ScanTime',
                   'Rewind','Loop','Dithers',
                   'ScanAmp','ScanRate','NodThrow','NodAngle','TotalTime','IntTime')
        intformat_func = partial(INT_CONVERTER,cols=intcols)
        detail = list(map(intformat_func,detail))
    return detail


def generate_details_tex(detail):
    '''Generate latex string of details table'''
    if isinstance(detail,(Table,dict)):
        # single detail table
        detail = [detail]

    texcodes = deque()
    
    for d in detail:
        # col align param must be special for boldface header line
        col_align = ['c']*len(d.colnames)
        col_align = '|^'.join(col_align)
        col_align = '|$%s|'%col_align

        preamble = r'\setlength{\tabcolsep}{0.25em}'
        #preamble += '\n\\centering\n\\captionsetup{justification=centering}'
        tablefoot = r'\\[0.5em]' + '\n' + r'\raggedleft{' + \
                    d.meta.get('footer','')+'}\n'+r'\setlength{\tabcolsep}{1em}'
        colnames = d.colnames.copy()

        # rename colnames to have headercolor
        newcols = {col:'\\cellcolor{headercolor}%s'%col for col in d.colnames}
        for col,newcol in newcols.items():
            #newcol = '\\cellcolor{headercolor}%s'%col
            d.rename_column(col,newcol)
            d.meta['units'][newcol] = d.meta['units'].get(col,'')

        # add gray color to unit cells
        units = {col:'\cellcolor{headercolor}%s'%d.meta['units'].get(col,'') for col in d.colnames}
        
        with StringIO() as f:
            d.write(f,format='latex',
                         latexdict={'header_start':r'\hline\rowstyle{\bfseries}',
                                    'tablealign':'h!',
                                    'caption':d.meta.get('caption',''),
                                    'col_align':col_align,
                                    'units':units,
                                    'data_end':'\hline',
                                    'preamble':preamble,
                                    'tablefoot':tablefoot})
            texcode = f.getvalue()

        # pull left if too long
        try:
            if (((max([len(name) for name in d[newcols['Name']]]) > 13) and ('ChopThrow' in colnames)) or (('ChopThrow' in colnames) and ('ScanAmp' in colnames))):
                texcode = texcode.replace(r'\begin{tabular}','\\hspace*{-1cm}\n\\begin{tabular}')
            elif 'NodThrow' in colnames or 'NodThw' in colnames:
                texcode = texcode.replace(r'\begin{tabular}','\\hspace*{-1cm}\n\\begin{tabular}')
            else:
                pass
        except KeyError:
            if 'Scale' in colnames:
                texcode = texcode.replace(r'\begin{tabular}','\\hspace*{-1cm}\n\\begin{tabular}')
            #texcode = texcode.replace(r'\begin{tabular}','\\hspace*{-1cm}\n\\begin{tabular}')

        #make small
        texcode = texcode.replace(r'\begin{tabular}','\\footnotesize\n\\begin{tabular}')
            
        # shrink whatever we can
        texcode = texcode.replace('Polarimetry','POL')
        texcode = texcode.replace('Lissajous','LIS')
        texcode = texcode.replace('Nod Time','Nod')
        texcode = texcode.replace('nan','')
        texcode = texcode.replace('None','')
        texcode = texcode.replace('Chop Throw','ChpT')
        texcode = texcode.replace('ChopThrow','ChpT')
        texcode = texcode.replace('Chop Ang','ChpA')
        texcode = texcode.replace('ChopAngle','ChpA')
        texcode = texcode.replace('Nod Throw','NodT')
        texcode = texcode.replace('NodThrow','NodT')
        texcode = texcode.replace('Nod Ang','NodA')
        texcode = texcode.replace('NodAngle','NodA')
        texcode = texcode.replace('Nod Time','NodT')
        if 'Nod Dwell' in colnames:
            texcode = texcode.replace('Nod Dwell','Nod')


        # strikethrough 'observed' aorIDs
        texcode = strikethrough(d.meta['observed'],texcode)
            
        texcodes.append(texcode)

    texcodes = '\n\\vspace*{-3em}\n'.join(list(texcodes))
    texcodes = '\\vspace*{-3em}\n%s' % texcodes
    return texcodes


def make_positions(tab, tex=True):
    '''Make position tables'''

    if tab[0]['aorID'] in ('99_9999_99','--'):
        return ''
    
    rows = deque()
    for t in tab:
        if t.get('RA') is None:
            continue
        aorid = t['aorID']
        #name = t['POSName'] if t.get('POSName') else t['Name']
        name = t['POSName'] if t.get('POSName') else t['target']
        coord = SkyCoord(ra=t['RA'],dec=t['DEC'],unit=(u.hourangle,u.deg))
        ra,dec = coord.to_string('hmsdms',sep=':',precision=2).split()
        order = t['order']
        rows.append((aorid,name,ra,dec,order))

    if not rows:
        return ''
    position = Table(rows=list(rows),names=('AORID','Name','RA','DEC','Order'))
    position.meta['caption'] = '\\captionline{Positions}{}'

    # if all positions are the same, remove duplicates
    origlen = len(position)
    position = unique(position,keys=('RA','DEC'))
    position.sort(['Order','AORID'])
    if len(position) == 1 and origlen != 1:
        position['AORID'][0] = '_'.join(position['AORID'][0].split('_')[0:2]) + '_*'

    position.remove_column('Order')
    if tex:
        for col in ['AORID','Name']:
            position[col].format = COL_FORMATTER

        # add strikethrough
        observed_dict = {t['aorID']:t.get('observed',False) for t in tab}
        position.meta['observed'] = observed_dict
        position = generate_pos_tex(position)
    else:
        position = position.to_pandas().to_dict('records')

    return position

def generate_pos_tex(position):
    '''Generate latex string of positions table'''

    if position is None:
        return ''

    # col align param must be special for boldface header line
    col_align = ['c']*len(position.colnames)
    col_align = '|^'.join(col_align)
    col_align = '|$%s|'%col_align

    preamble = '\\centering\n\\captionsetup{justification=centering}'

    # rename colnames to have headercolor
    for col in position.colnames:
        newcol = '\\cellcolor{headercolor}%s'%col
        position.rename_column(col,newcol)

    with StringIO() as f:
        position.write(f,format='latex',
                       latexdict={'header_start':r'\hline\rowstyle{\bfseries}',
                                  'tablealign':'h!',
                                  'caption':position.meta['caption'],
                                  'preamble':preamble,
                                  'col_align':col_align,
                                  'data_end':'\hline'})
        texcode = f.getvalue()

    # strikethrough 'observed' aorIDs
    texcode = strikethrough(position.meta['observed'],texcode)

    return texcode


def strikethrough(aordict,texcode):
    for aorid,obs in aordict.items():
        if obs:
            a = aorid.replace('_','\\\_')
            pattern = r'(%s\s\&.*\s\\\\)'%a
            texcode = re.sub(pattern,STRIKETHROUGH_REPL,texcode)
    return texcode


def get_pos_bundle(tab, dcs, odir):
    '''Download pos bundle and filter by included AORs'''
    posnames = [t['POSName'] if t.get('POSName') else t['Name'] for t in tab]
    fpid = tab[0]['FlightPlan']
    posfiles = ['./%s/%s_%s.pos'%(fpid,fpid,x) for x in posnames]

    postarfile = dcs.getPOSBundle(fpid)
    if postarfile is None:
        return None
    with tarfile.open(postarfile) as t:
        members = t.getmembers()
        members = list(filter(lambda x: x.name in posfiles, members))
        t.extractall(odir,members=members)
        # move to directory above
        for fname in Path(odir/fpid).glob('*.pos'):
            dest = fname.parent.parent/(fname.name.split(fpid)[1][1:])
            shutil.move(fname,dest)
        try:
            Path(odir/fpid).rmdir()
        except FileNotFoundError:
            pass
    return postarfile

def generate_overlays(table):
    # add row index for color cycle
    for idx,row in enumerate(table):
        row['cidx'] = idx
    overlays = list(map(generate_overlay,table))

    for overlay,row in zip(overlays,table):
        row['overlay'] = overlay
        '''
        if isinstance(overlay,dict):
            row['overlay'] = overlay
        else:
            # FIFI has two overlays per aorid
            row['overlayB'] = overlay[0]
            row['overlayR'] = overlay[1]
        '''
    return table

def generate_overlay(row,nod=True,dithers=True,FIFI_label=None):
    #tab = unique(tab,keys=['RA_aor','DEC_aor'])

    if row['aorID'] in ('99_9999_99','--'):
        return None
    
    try:
        coord = SkyCoord(ra=row['RA'],dec=row['DEC'],unit=(u.hourangle,u.deg))
    except ValueError:
        # likely a solar system object
        return None
    
    # get roll angle
    rolls = (row['ROF_start'],row['ROF_end'])

    # override displayed roll if specified in config
    roll = row.get('roll')
    if isinstance(roll,bool) and roll:
        if row['InstrumentName'] == 'FIFI-LS':
            roll = float(row['MapRotationAngle'])*u.deg
        else:
            roll = float(rolls[0])*u.deg
    if isinstance(roll,(float,int,np.float,np.int)):
        roll = roll*u.deg
    elif isinstance(roll,u.Quantity):
        roll = roll
    elif isinstance(roll,str):
        roll = u.Quantity(np.float(roll),u.deg)
    else:
        if row['InstrumentName'] == 'FIFI-LS':
            roll = float(row['MapRotationAngle'])*u.deg
        else:
            roll = float(rolls[0])*u.deg

    try:
        TARFoffset = TARFOFFSET[row['InstrumentName']]
    except KeyError:
        # likely wrong instrument?
        TARFoffset = 0*u.deg

    # get band and mode for FOV
    band = row['InstrumentSpectralElement1'].split('_')[-1] if 'FIF' not in row['InstrumentSpectralElement1'] else row['InstrumentSpectralElement1']
    mode = row['ObsPlanConfig']

    if row['ObsPlanMode'] == 'C2N' and mode == 'TOTAL_INTENSITY':
        mode = 'C2N'
    else:
        mode = 'TOT' if mode == 'TOTAL_INTENSITY' else 'POL'

    label = '%s_%s'%(band,mode)
    try:
        fov = FOV[label]
    except KeyError:
        # assume FORCAST
        if 'G' in label:
            fov = FOV['FORCAST_GSM']
        else:
            fov = FOV['FORCAST_IMG']

    if row['InstrumentName'] == 'FIFI-LS':
        if FIFI_label:
            labels = [FIFI_label]
        else:
            # make both FIFI cameras
            labels = ['FIF_BLUE','FIF_RED']
        overlays = deque()
        for label in labels:
            width,height = [f*u.arcmin for f in FOV[label]]
            name = '%s %s' % (row['Name'],coord.to_string('hmsdms',precision=2,sep=':'))
            if label == 'FIF_RED' and not FIFI_label:
                # need to tweak FIFI red boresite
                coord = coord.directional_offset_by(roll+TARFoffset+90*u.deg,-0.162*u.arcmin)

            
            o = make_box(coord,width,height,roll,TARFoffset,label=label,name=name,
                         color=COLORS[row['cidx']%len(COLORS)],split=False,aorid=row['aorID'])
            overlays.append(o)
        overlays = list(overlays)
    else:
        # HAWC and FORCAST
        width,height = [f*u.arcmin for f in fov]
        name = '%s %s' % (row['Name'],coord.to_string('hmsdms',precision=2,sep=':'))

        split = True if mode == 'TOT' and label in FOV else False

        if label not in FOV:
            # FORCAST
            label = label.replace('_TOT','')
            label = label.replace('_POL','')

        overlays = [make_box(coord,width,height,roll,TARFoffset,label=label,name=name,
                             color=COLORS[row['cidx']%len(COLORS)],split=split,aorid=row['aorID'])]

    for overlay in overlays:
        overlay['roll'] = (float(rolls[0])*u.deg,float(rolls[1])*u.deg)


        if row['NodType'] != 'OTFMAP':
            # we are chop/nod dithering
            if dithers and row['DitherPattern'] not in (None,'None'):
                if row['ChopAngleCoordinate'] == 'Sky':
                    scale = row['DitherScale']*u.arcsec
                else:
                    try:
                        scale = row['DitherScale']*PIXSIZE[label]
                    except KeyError:
                        if 'G' in label:
                            scale = row['DitherScale']*PIXSIZE['FORCAST_GSM']
                        else:
                            scale = row['DitherScale']*PIXSIZE['FORCAST_IMG']

                diths = make_dithers(overlay['center'],scale=scale,angle=roll)

                overlay['dithers'] = [make_box(dith, width, height, angle=roll, TARFoffset=TARFoffset, label=label, split=split, color=overlay['color'], reglabel='_d',name=name) for dith in diths]

            if nod and row['ObsPlanMode'] == 'C2NC2':
                chopthrow = row['ChopThrow']*u.arcsec
                chopangle = row['ChopAngle']*u.deg
                nodthrow = row['NodThrow']*u.arcsec
                nodangle = row['NodAngle']*u.deg

                if row['ChopAngleCoordinate'] == 'Array':
                    chopangle += roll

                nodAchopB,nodBchopA,nodBchopB = make_C2NC2(overlay['center'],
                                                           chopthrow=chopthrow,chopangle=chopangle,
                                                           nodthrow=nodthrow,nodangle=nodangle)
                nodAchopBdict = row.copy()
                nodBchopAdict = row.copy()
                nodBchopBdict = row.copy()

                for ntab,n in zip((nodAchopBdict,nodBchopAdict,nodBchopBdict),(nodAchopB,nodBchopA,nodBchopB)):
                    ra,dec = n.to_string('hmsdms').split()
                    ntab['RA'] = ra
                    ntab['DEC'] = dec

                overlay['nods'] = [generate_overlay(n, nod=False,dithers=False) \
                                   for n in (nodAchopBdict,nodBchopAdict,nodBchopBdict)]

            elif nod:
                chopthrow = row['ChopThrow']*u.arcsec
                chopangle = row['ChopAngle']*u.deg

                if row['ChopAngleCoordinate'] == 'Array':
                    chopangle += roll

                nodA,nodB = make_NMC(overlay['center'],
                                     chopthrow=chopthrow,
                                     chopangle=chopangle,label=overlay['label'])

                nodAdict = row.copy()
                nodBdict = row.copy()
                
                ra,dec = nodA.to_string('hmsdms').split()
                nodAdict['RA'] = ra
                nodAdict['DEC'] = dec
                ra,dec = nodB.to_string('hmsdms').split()
                nodBdict['RA'] = ra
                nodBdict['DEC'] = dec

                if row['InstrumentName'] == 'FIFI-LS':
                    # only make current label
                    fiflabel = overlay['label']
                else:
                    fiflabel = None
                overlay['nods'] = [generate_overlay(n, nod=False,dithers=False, FIFI_label=fiflabel) \
                                   for n in (nodAdict,nodBdict)]


        else:
            ampx,ampy = u.Quantity(row['ScanAmplitudeEL'],u.arcsec), \
                        u.Quantity(row['ScanAmplitudeXEL'],u.arcsec)
            ampx += width
            ampy += height
            overlay['dithers'] = [make_box(coord,ampx,ampy,roll,TARFoffset,label=label,name=name,
                                           color=COLORS[row['cidx']%len(COLORS)],scan=True,reglabel='scan',
                                           aorid=row['aorID'])]

        overlay['aorID'] = row['aorID']
        overlay['InstrumentName'] = row['InstrumentName']
    return overlays

def get_overlay_params(tab):
    #####UNUSED
    overlays = deque()
    for idx,row in enumerate(tab):
        if 'IMGOVERRIDES' in tab.meta:
            key = 'Leg%i__%s'%(row['Leg'],row['aorID'])
            if key in tab['IMGOVERRIDES']:
                roll = u.Quantity(tab['IMGOVERRIDES'][key],u.deg)
            elif tab['IMGOVERRIDES'].get('roll',False):
                roll = u.Quantity(row['Angle'],u.deg)
            else:
                roll = None
        else:
            roll = None
        overlay = generate_overlay(tab,idx=idx,roll=roll,TARFoffset=TARFoffset)
        overlays.append(overlay)
        
    overlays = [overlay for overlay in overlays if overlay is not None]
    return overlays

def get_recenter_image(overlaylist):
    '''Get longest wavelength "recenter"'''
    recenter = None
    for overlay in overlaylist:
        if 'recenter' not in overlay:
            continue
        recenter = overlay['recenter']

    return recenter

def make_figures(table,fdir,reg=False,guidestars=None,irsurvey=None,savefits=False,fpi=False,**kwargs):
    '''Generate figure from each overlay'''

    # make defaults from table
    vrows = list(filter(lambda row:row.get('overlay'),table))
    imgoptions = IMGOPTIONS.copy()
    if vrows and vrows[0]['InstrumentName'] == 'FIFI-LS':
        imgoptions['width'] = 0.2*u.deg
        imgoptions['height'] = 0.2*u.deg
    overlays = [{k:row.get(k,v) for k,v in imgoptions.items()} for row in vrows]

    for o,row in zip(overlays,vrows):
        if row.get('overlay'):
            if isinstance(row['overlay'],dict):
                o.update(row['overlay'])
            else:
                # if there is a second overlay (FIFI), add it
                if len(row['overlay']) == 2:
                    vals = o.copy()
                    o['overlay2'] = vals
                    o.update(row['overlay'][0])
                    o['overlay2'].update(row['overlay'][1])
                else:
                    o.update(row['overlay'][0])


    #overlays = [row['overlay'] for row in table if row.get('overlay')]
    
    # remove any with 'nofigure' flag
    overlays = list(filter(lambda o:not o.get('nofigure',False),overlays))
    if not overlays:
        return None

    # image options are grabbed from first row in blk
    options = overlays[0].copy()
    if irsurvey is not None:
        options['irsurvey'] = irsurvey

    if fpi:
        options['fpi'] = fpi

    if 'width' in options and isinstance(options['width'],str):
        options['width'] = u.Quantity(float(options['width']),u.deg)
    if 'height' in options and isinstance(options['height'],str):
        options['height'] = u.Quantity(float(options['height']),u.deg)

    try:
        fig,hdu = get_image(overlays,**options)
    except TypeError:
        warnings.warn('Issue querying SkyView. Skipping figure.',RuntimeWarning)
        return None

    if fig is None:
        return None

    if guidestars is not None:
        # get hdu footprint
        try:
            footprint = WCS(hdu).calc_footprint()
        except AttributeError:
            # likely hdu is an hdu list from irhdu
            footprint = WCS(hdu[0]).calc_footprint()
        box = Polygon(footprint)
        
        guides = guidestars[table[0]['ObsBlkID']]
        try:
            guidecoord = SkyCoord([g['COORD'] for g in guides])
        except IndexError:
            guidecoord = None

        if guidecoord:
            # do this check in case the above try/except fails
            fig.show_markers(guidecoord.ra,guidecoord.dec,
                             marker='o',s=80,
                             linewidths=2,edgecolor=GCOLOR)

            points = MultiPoint([(g.ra.value,g.dec.value) for g in guidecoord])
            inside = (box.contains(p) for p in points)

            gregs = deque()

            for g,i in zip(guides,inside):
                if not i:
                    continue
                fig.add_label(g['COORD'].ra.value, g['COORD'].dec.value-0.01, g['Name'],layer=g['Name'],size=8)
                txt = fig._layers[g['Name']]
                txt.set_path_effects([path_effects.Stroke(linewidth=1, foreground='white'),
                                      path_effects.Normal()])
                greg = PointSkyRegion(g['COORD'],
                                      meta=RegionMeta({'label':' '.join((g['Name'],g['Imager'],g['Catalog']))}),
                                      visual=RegionVisual({'color':GCOLOR}))
                gregs.append(greg)
            
        
    '''
        guides = guidestars[table[0]['ObsBlkID']]
        guides = list(filter(lambda g: g['Radius'] < np.hypot(table[0]['width'],table[0]['height']),guides))
        guidecoord = SkyCoord([g['COORD'] for g in guides])
        fig.show_markers(guidecoord.ra,guidecoord.dec,
                         marker='o',s=80,
                         linewidths=2)#edgecolor='#FFD700')
        for g in guides:
            fig.add_label(g['COORD'].ra.value, g['COORD'].dec.value, g['Name'])
    '''

    # show chop/nod/dithers
    for o in overlays:
        if 'dithers' in o:
            for dith in o['dithers']:
                fig.show_polygons(dith['box'],edgecolor=dith['color'],lw=1,
                                  linestyle='dashed',alpha=0.7)
        if 'nods' in o:
            for nod in o['nods']:
                try:
                    fig.show_polygons(nod['box'],edgecolor=o['color'],lw=1.5,
                                      linestyle='dotted',alpha=0.7)
                except TypeError:
                    for n in nod:
                        fig.show_polygons(n['box'],edgecolor=o['color'],lw=1.5,
                                          linestyle='dotted',alpha=0.7)

        if 'overlay2' in o:
            o = o.get('overlay2','')
            if 'dithers' in o:
                for dith in o['dithers']:
                    fig.show_polygons(dith['box'],edgecolor=dith['color'],lw=1,
                                      linestyle='dashed',alpha=0.7)
            if 'nods' in o:
                for nod in o['nods']:
                    try:
                        fig.show_polygons(nod['box'],edgecolor=o['color'],lw=1.5,
                                          linestyle='dotted',alpha=0.7)
                    except TypeError:
                        for n in nod:
                            fig.show_polygons(n['box'],edgecolor=o['color'],lw=1.5,
                                              linestyle='dotted',alpha=0.7)
            


    #outfile = fdir/('Leg%02d.png'%tab['Leg'][0])
    fdir.mkdir(exist_ok=True)
    outfile = fdir/('Leg%02d.pdf'%table[0]['Leg'])
    
    fig.savefig(outfile,dpi=300)
    fig.close()

    if savefits:
        fitsdir = fdir.parent/'images'
        fitsdir.mkdir(exist_ok=True)
        aorid = table[0]['planID']
        if not isinstance(hdu,fits.HDUList):
            hdu = [hdu]
        for h in hdu:
            surv = h.header['SURVEY'].strip().replace(' ','_')
            fname = '%s_%s_%s.fits'%(outfile.stem,aorid,surv)
            h.writeto(fitsdir/fname,overwrite=True,output_verify='silentfix+ignore')

    if reg:
        regs = deque()
        for row in table:
            reg = row['overlay'].get('reg')
            if isinstance(reg,str):
                reg = [reg]
            if reg:
                regs.extend(reg)
            if row['overlay'].get('dithers'):
                dithreg = [d.get('reg') for d in row['overlay']['dithers'] if d.get('reg')]
                if dithreg:
                    regs.extend(dithreg)
                
        regs = [DS9Parser(r).shapes.to_regions()[0] for r in set(list(regs))]
        if guidestars is not None:
            regs += gregs
        
        rdir = fdir.parent/'reg'
        rdir.mkdir(exist_ok=True)
        aorid = table[0]['planID']
        regfile = '%s_%s.reg'%(outfile.stem,aorid)
        regfile = rdir/regfile
        #with open(regfile,'w') as f:
        #    f.write('\n'.join(regs))
        write_ds9(regs,regfile)
        
    
    return outfile.relative_to(fdir.parent)

def make_comments(table):
    '''Generate tex for obsblk comments'''

    if table[0]['aorID'] in ('99_9999_99','--'):
        comment = table[0]['ObsBlkComment']
        if comment in (None,'None'):
            comment = ''
        comment = utf8tolatex(comment)
        comment = comment.replace('\n',r'\\')
        comment = comment.replace(r'{\textbackslash}{\textbackslash}',r'\\')
        return comment

    comments = table[0].get('ObsBlkComment','')
    if not comments:
        return ''

    #safe convert to latex
    comments = utf8tolatex(comments)
    comments = comments.replace(r'{\textbackslash}{\textbackslash}',r'\\')
    return comments

def hawc_sio_comments(table):
    '''Assign comment block based on mode'''

    if table[0]['aorID'] in ('99_9999_99','--'):
        comment = table[0]['ObsBlkComment']
        comment = utf8tolatex(comment)
        comment = comment.replace('\n',r'\\')
        comment = comment.replace(r'{\textbackslash}{\textbackslash}',r'\\')
        return comment
    
    head = r'Procedure for Instrument Operator: \\'
    
    if any([((row['ObsPlanConfig'] == 'POLARIZATION') and (row['ObsPlanMode'] == 'OTFMAP')) for row in table]):
        mode = 'LISPOL'
    elif any([((row['ObsPlanConfig'] == 'POLARIZATION') and (row['ObsPlanMode'] == 'C2N')) for row in table]):
        mode = 'Polarimetry'
    else:
        mode = 'Lissajous'

    comment = '%s%s' % (head, HAWC_SIO[mode])
    
    return comment

def copy_comments(filename):
    """Return comments from .tex file, if they exist"""
    # first find existing .tex file
    filename = Path(filename)
    if filename.exists():
        with open(filename,'r') as f:
            text = f.read()
    else:
        # try going back to one older version
        parents = [str(p) for p in filename.parents]
        top = parents[-2]
        if top.split('_')[-1][0] == 'v':
            # this is versioned
            version = top[-1]
            name = parents[-3].split('/')[-1]
            try:
                version = int(version)
                prev_version = version - 1
            except ValueError:
                prev_version = chr(ord(version) - 1)
            newparent = '%s%s'%(top[:-1],prev_version)
            newpath = Path(newparent).joinpath(name)/filename.name

            if Path(newpath).exists():
                filename = Path(newpath)
                with open(filename,'r') as f:
                    text = f.read()
            else:
                return None
            
        else:
            return None
    
    if r'% <COMMENT legname' not in text:
        return None

    print('Copying comments from %s' % filename)
    comments = COMMENT_RE.findall(text)
    comments = [x[1] for x in comments]
    return comments


def match_FAORs(tables,dcs):
    aorIDs = set((row['aorID'] for table in tables for row in table))
    faors = dcs.getFAORs(aorIDs,match=True)

    for table in ProgressBar(tables):
        for row in table:
            match = faors.get(row['aorID'])
            if match:
                row.update(match)
    return tables

def match_SCTs(tables,dcs):
    aorIDs = set((row['aorID'] for table in tables for row in table))
    scts = dcs.getSCTs(aorIDs,match=True)

    for table in ProgressBar(tables):
        for row in table:
            match = scts.get(row['aorID'])
            if match:
                row.update(match)
    return tables

def get_FAOR_map(faordir,keys=None,label=''):
    '''Locate faor files and return file mapping with aorid'''
    #print('Locating .faor files in %s...'%faordir)
    faorfiles = faordir.glob('*%s.faor'%label)
    faordict = {}
    for ffile in faorfiles:
        faor = FAOR.read(ffile)
        aorids = (cfg['AORID'] for cfg in faor.config)
        leg = int(ffile.name.split('__')[0][-2:])
        if keys is None:
            fdict = {'Leg%02d__%s'%(leg,aorid):ffile for aorid in aorids}
        elif isinstance(keys,str):
            fdict = {'Leg%02d__%s'%(leg,aorid):ffile for aorid in aorids if 'Leg%02d__%s'%(leg,aorid) == keys}
        else:
            fdict = {'Leg%02d__%s'%(leg,aorid):ffile for aorid in aorids if 'Leg%02d__%s'%(leg,aorid) in keys}
            
        faordict.update(**fdict)

    return faordict


def update_with_cfg(table, cfg):
    """Update table rows with cfg"""

    try:
        key = table[0]['fkey']
        blk = table[0]['ObsBlkID']
    except KeyError:
        return table

    # first, check if 'keep' keyword is present in fkey or obsblk
    if cfg.has_section(key) and 'aorids' in cfg[key]:
        keep = [x.strip() for x in cfg[key].get('aorids').split(',')]
    elif cfg.has_section(blk) and cfg.has_option(blk,'aorids'):
        keep = [x.strip() for x in cfg[blk].get('aorids').split(',')]
    else:
        keep = None

    # filter rows with 'keep' aorids
    if keep:
        tab = list(filter(lambda row: row['aorID'] in keep,table))
        if not tab:
            warnings.warn('All aorids would be removed from this AOR. Check cfg file. Cfg file ignored.', RuntimeWarning)
        else:
            table = tab

    cmaps = map(lambda aor: CMap.from_AOR_dict(cfg,aor),table)

    for row,cmap in zip(table,cmaps):
        # add other options
        imgoptions = IMGOPTIONS.copy()
        cmap.add_map(imgoptions)
        row.update(cmap.as_dict())

    return table


def write_tex_dossier(tables, name,title,filename,
                      template='template.tex',
                      config=None,
                      mconfig=None,
                      refresh_cache=False,
                      faor=False,
                      sct=False,
                      posfiles=False,
                      reg=False,
                      dcs=None,local=None,
                      sio=False,
                      fpi=False,
                      savefits=False,
                      irsurvey=None,
                      preserve_comments=False,
                      no_figure=False,
                      tex=True,
                      writetex=True):
    '''Write dossier pages for each table in tables'''

    # get template
    #template = latex_jinja_env.get_template(template)
    template = get_latex_template(template)

    if not dcs:
        # initialize DCS link
        dcs = DCS(refresh_cache=refresh_cache,modelcfg=mconfig)

    # process config file
    if config:
        cfg = ConfigParser()
        cfg.read(config)
        # merge config with rows
        tables = list(map(lambda t:update_with_cfg(t,cfg),tables))

    if faor:
        # merge faor information
        print('Matching faors to AORIDs...')
        tables = match_FAORs(tables,dcs)

    if sct:
        # merge sct information
        print('Matching scts to AORIDs...')
        tables = match_SCTs(tables,dcs)
        
    #tables = [[table] if isinstance(table,dict) else table for table in tables]
    #for table in tables:
    #    if 'aorID' not in table[0]:
    #        print(table)
    #exit()
        
    # make overview tables
    print('Generating overview...')
    over_func = partial(make_overview,tex=tex)
    overviews = ProgressBar.map(over_func,tables,multiprocess=MP)
    legnames = ['Leg %i (%s)'%(table[0]['Leg'],table[0]['Name'].replace('_','\_')) for table in tables]

    # make details tables
    print('Writing detail table...')
    detail_func = partial(make_details,tex=tex)
    details = ProgressBar.map(detail_func,tables,multiprocess=MP)

    # make pos tables
    print('Writing position table...')
    pos_func = partial(make_positions,tex=tex)
    positions = ProgressBar.map(pos_func,tables,multiprocess=MP)

    if posfiles:
        # copy pos files from cache to local dir
        print('Copying posfiles...')
        pdir = Path(filename).parent/'pos'
        pdir.mkdir(exist_ok=True)
        for tab in ProgressBar(tables):
            get_pos_bundle(tab,dcs,pdir)

    if no_figure:
        imagefiles = [None]*len(tables)
    else:
        # make figures
        print('Generating overlays...')
        tables = ProgressBar.map(generate_overlays,tables,multiprocess=MP)

        # get guidestars
        blkids = {row['ObsBlkID'] for table in tables for row in table \
                  if row['ObsBlkID'] not in ('--','None',None)}
        guides = dcs.getGuideStars(list(blkids))
        guidestars = defaultdict(list)
        for guide in guides:
            coord = SkyCoord(ra=guide['RA'],dec=guide['Dec'],unit=(u.deg,u.deg))
            guide['COORD'] = coord
            guidestars[guide['ObsBlkID']].append(guide)

    
        # output directory for figures
        fdir = Path(filename).parent/'figs'
        # define partial function for multiprocessing
        figfunc = partial(make_figures,fdir=fdir,reg=reg,guidestars=guidestars,
                      irsurvey=irsurvey,savefits=savefits,fpi=fpi)
        print('Generating figures...')
        imagefiles = ProgressBar.map(figfunc,tables,multiprocess=MP)

    if preserve_comments:
        oldcomments = copy_comments(filename)
        if oldcomments is None:
            oldcomments = ['']*len(tables)
    else:
        oldcomments = ['']*len(tables)

    if sio:
        print('Generating comments...')
        comments = ProgressBar.map(hawc_sio_comments,tables,multiprocess=MP)
    else:
        print('Gathering comments...')
        comments = ProgressBar.map(make_comments,tables,multiprocess=MP)

    # merge old and new
    comments = [old if old else new for old,new in zip(oldcomments,comments)]

    # finally, render latex template
    render = {'flightName':name.replace('_','\_'),
              'title':title,
              'tables':zip(legnames,overviews,
                           details,positions,imagefiles,comments)}

    pages = template.render(**render)

    if writetex:
        with open(filename,'w') as f:
            f.write(pages)

    return filename
