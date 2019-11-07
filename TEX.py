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
#from splinter import Browser
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
from CMap import CMap
#import MSX
import cloudpickle
import pickle
#from planner import split_leg
from requests.exceptions import ChunkedEncodingError, SSLError, ConnectionError
from regions import RectangleSkyRegion,RegionMeta,RegionVisual,write_ds9,ds9_objects_to_string,DS9Parser
from pylatexenc.latexencode import utf8tolatex
import tarfile
import shutil
import warnings
warnings.filterwarnings('ignore',category=MatplotlibDeprecationWarning)
warnings.filterwarnings('ignore',category=AstropyDeprecationWarning)
warnings.filterwarnings('ignore',category=AstropyUserWarning)
warnings.filterwarnings('ignore',category=AstropyWarning)
np.warnings.filterwarnings('ignore')

DEBUG = False

if DEBUG is False:
    log.disable_warnings_logging()
    log.setLevel('ERROR')

COL_FORMATTER = lambda x: x.replace('_','\_')
INT_FORMATTER = lambda x: '%i'%int(np.float(x)) if x not in (None,'None',-9999,'--') else ''

def INT_CONVERTER(row,cols):
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

#COLORS = ['#ff0000','#00ff00','#0000ff']
COLORS = ['#d62728','#1f77b4','#2ca02c','#ff7f0e','#9467bd','#17becf','#e377c2']

FOV = {'A_TOT':(2.8,1.7),'C_TOT':(4.2,2.7),'D_TOT':(7.4,4.6),'E_TOT':(10.0,6.3),
       'A_POL':(1.4,1.7),'C_POL':(2.1,2.7),'D_POL':(3.7,4.6),'E_POL':(5.0,6.3),
       'FORCAST_IMG':(3.4,3.2),'FORCAST_GSM':(.04,3.18)}

PIXSIZE = {'A_TOT':2.34, 'C_TOT':4.02,'D_TOT':6.90,'E_TOT':8.62,'FORCAST_IMG':0.768,'FORCAST_GSM':0.768}

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
    loader = FileSystemLoader(os.path.abspath('.'))
)


IMGOPTIONS = {'width':0.3, 'height':0.3,
              'survey':'DSS2 Red',
              'vmin':None, 'vmax':None,
              'recenter':None,'roll':True,
              'invert':True,'irsurvey':None,
              'compass':True,'nofigure':False}

TARFOFFSET = {'HAWC_PLUS':3*u.deg,
              'FORCAST':40.6*u.deg}


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

        #recenter = center.directional_offset_by(90*u.deg+offset, r_width/2+splitgap/2)
        recenter = center.directional_offset_by(-90*u.deg+offset,r_width/2+splitgap/2)
        
        r0 = make_box(r0_center,r_width,height,angle,TARFoffset,label,
                      linewidth,color,name)
        r1 = make_box(r1_center,r_width,height,angle,TARFoffset,label,
                      linewidth,color,name)

        boxdict = {'box':[r0['box'][0],r1['box'][0]],'linewidth':linewidth,
                   'color':color,'center':center,'label':label,'name':name,
                   'recenter':recenter}

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
                        linewidth,color,name,split=False,scan=False)

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
    try:
        center = center[0]
    except TypeError:
        pass
    if np.isnan(center.ra) or np.isnan(center.dec):
        boxdict['reg'] = None
    elif name is None:
        boxdict['reg'] = None
    else:
        # split name---last two entries are ra/dec
        name = ' '.join(name.split()[0:-2])
        if 'aorid' in kwargs:
            # IGNORE 'aorid'---messes up the comment format
            #reg = RectangleSkyRegion(center,width,height,angle=offset,meta=RegionMeta({'label':name,
            #                                                                           'comment':'(AORID:%s)'%kwargs['aorid']}))
            reg = RectangleSkyRegion(center,width,height,angle=offset,meta=RegionMeta({'label':name}))
        else:
            reg = RectangleSkyRegion(center,width,height,angle=offset,meta=RegionMeta({'label':name}))
        boxdict['reg'] = ds9_objects_to_string([reg])
    return boxdict


def make_dithers(center,scale,angle=0*u.deg):
    '''Generate box for dithers'''
    diag = np.hypot(scale,scale)
    posangle = angle+45*u.deg
    tl = center.directional_offset_by(posangle,diag)
    tr = center.directional_offset_by(posangle-90*u.deg,diag)
    bl = center.directional_offset_by(posangle-90*u.deg,-diag)
    br = center.directional_offset_by(posangle,-diag)
    return (tl,tr,bl,br)

def make_NMC(center,chopthrow=300*u.arcsec,chopangle=0*u.arcsec):
    '''Given a center and chop nod parameters, calculate the nod throws'''
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

    if recenter:
        if isinstance(recenter,SkyCoord):
            center = recenter
        else:
            center = SkyCoord(recenter,unit=(u.hourangle,u.deg))

    else:
        center = overlays[0]['center']

    try:
        im = SkyView.get_images(center,survey=survey,width=width,height=height,
                                show_progress=DEBUG)
    except (SSLError,ChunkedEncodingError,ConnectionError):
        warnings.warn('Cannot query SkyView service. Skipping image.')
        return None

    try:
        hdu = im[0][0]
    except (IndexError,TypeError):
        warnings.warn('Cannot process SkyView response. Skipping image.')
        return None
        
    fig = FITSFigure(hdu)
    fig.show_grayscale(vmin=vmin,vmax=vmax,invert=invert)

    regs = deque()

    for idx,overlay in enumerate(overlays):
        fig.show_polygons(overlay['box'],edgecolor=overlay['color'],lw=overlay['linewidth'])
        fig.show_markers(overlay['center'].ra.value,overlay['center'].dec.value,marker='*',edgecolor=overlay['color'])

        if 'reg' in overlay and overlay['reg'] and overlay['reg'].center is not None:
            regs.append(overlay['reg'])

        if overlay['label']:
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
                             fpiradius,edgecolor=overlay['color'],
                             linestyle='dashed',lw=1,
                             alpha=0.4)
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


    regs = list(filter(lambda x: x is not None,regs))
    # remove duplicates
    regs = set(regs)   # NOTE: this effectively does nothing if aorids are given to make_box, since each line will be unique
    
    # convert back into objects--this is all because regions don't pickle
    regs = [DS9Parser(reg).shapes.to_regions()[0] for reg in regs]

    # add ir image
    if irsurvey is not None:
        if irsurvey == 'msx':
            dfile = MSX.query_region(center,show_progress=True)
            if dfile is not None:
                irhdu = fits.open(dfile)
                irhdu[0].header['SURVEY'] = 'MSX Band E'
                hdu = fits.HDUList([hdu,irhdu[0]])
        else:
            try:
                im = SkyView.get_images(center,survey=irsurvey,width=width,height=height,
                                        show_progress=DEBUG)
                irhdu = im[0][0]
                hdu = fits.HDUList([hdu,irhdu])
            except (SSLError,ChunkedEncodingError,ConnectionError,IndexError,TypeError):
                pass
    
    return fig, regs, hdu
    

def make_overview(leg, tex=True):
    '''Make table of overview stats for leg'''

    # overview cols to extract
    ocols = ('Start','ObsDur','Target','ObsBlk','Priority','RA','DEC')

    # metacols
    hcols = ('Leg','Name','PI')
    mcols = ('Elev','ROF','ROFRT','MoonAngle','MoonIllum','THdg','THdgRT')

    l = leg[0]
    overview = {k:l.get(k,'') for k in ocols}

    # make metadata
    if tex:
        overview['header'] = '\\captionline{Leg %i (%s)}{%s}' % (l['Leg'],l['Name'],l['PI'])
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
            footer1 = '\quad '.join((footer1,'Priority %s'%prior))
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

    # safe convert obsblkid
    try:
        overview['ObsBlk'] = [blk.replace('_','\_') for blk in overview['ObsBlk']]
    except AttributeError:
        overview['ObsBlk'] == ''
    
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

    instrument = tab[0]['InstrumentName']
    if instrument == 'HAWC_PLUS':
        keys = ('ObsPlanConfig','aorID','Name','InstrumentSpectralElement1','Repeat','NodTime','ChopThrow','ChopAngle','ScanTime','ScanAmplitudeEL','ScanAmplitudeXEL','ScanRate','ChopAngleCoordinate')
        key_map = {'ObsPlanConfig':'Mode','aorID':'AORID','ChopAngleCoordinate':'Sys','InstrumentSpectralElement1':'Band/Bore','ScanAmplitudeEL':'ScanAmp'}

        # replace some values
        for t in tab:
            # change coordsys
            sys = t['ChopAngleCoordinate']
            t['ChopAngleCoordinate'] = 'ERF' if sys == 'Sky' else 'SIRF'

            # store filter for boresite
            filt = t['InstrumentSpectralElement1'][-1]

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
                    t['InstrumentSpectralElement1'] = '/'.join((filt,filt))
                else:
                    t['ObsPlanConfig'] = '?'

            # polarimetry
            else:
                t['ObsPlanConfig'] = 'POL'
                t['InstrumentSpectralElement1'] = '/'.join((filt,filt))
                
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
        if all(mode == 'POL' for mode in detail['Mode']):
            detail.remove_columns(('ScanTime','ScanAmp','ScanRate'))

        # if any dithering, make dither footer
        unit_map = {'Sky':'arcsec','Array':'pix'}
        if any((t.get('DitherPattern') for t in tab)):
            dithscale = (t.get('DitherScale') for t in tab)
            dithunit = (unit_map.get(t.get('DitherCoord')) for t in tab)
            dithband = (t['InstrumentSpectralElement1'][-1] for t in tab)
            footer = ['\t%s: %i %s' % (band,scale,unit) for band,scale,unit \
                      in zip(dithband,dithscale,dithunit) if scale]
            footer = set(footer)
            if len(footer) == 1:
                footer = footer.pop()
            else:
                footer = '\\\\\n'.join(sorted(footer))
            footer = 'dither\quad\quad %s\\\\'%footer
            detail.meta['footer'] = footer

    

    elif instrument == 'FORCAST':
        keys = ['ObsPlanConfig','ObsPlanMode','aorID','Name',
                'InstrumentSpectralElement1','InstrumentSpectralElement2',
                'Repeat','NodTime','TotalTime',
                'ChopThrow','ChopAngle','ChopAngleCoordinate',
                'NodThrow','NodAngle']
        key_map = {'ObsPlanConfig':'Mode','aorID':'AORID','ObsPlanMode':'Type','ChopAngleCoordinate':'Sys','InstrumentSpectralElement1':'Band','InstrumentSpectralElement2':'Slit'}
        faor_keys = ['Nod','Dithers','Scale','IntTime','FDUR','TLOS','TLSPN',
                     'Rewind','Loop']

        for t in tab:
            sys = t['ChopAngleCoordinate']
            t['ChopAngleCoordinate'] = 'ERF' if sys == 'Sky' else 'SIRF'

            # shorten filter config
            t['InstrumentSpectralElement1'] = t['InstrumentSpectralElement1'].replace('FOR_','')
            t['InstrumentSpectralElement2'] = t['InstrumentSpectralElement2'].replace('FOR_','')

            # combine filters if dual mode
            if 'DUAL' in t.get('ObsPlanConfig',''):
                t['InstrumentSpectralElement1'] = '/'.join((t['InstrumentSpectralElement1'],
                                                            t['InstrumentSpectralElement2']))
                t['InstrumentSpectralElement2'] = None

            # drop second element if OPEN
            elif t['InstrumentSpectralElement2'] == 'OPEN':
                t['InstrumentSpectralElement2'] = None

            # shorten chop mode
            if t['NodType'] == 'Nod_Match_Chop':
                t['ObsPlanMode'] = 'NMC'

        # keep certain keys
        if 'FAORfile' in tab[0]:
            keys += faor_keys
        detail = [{key:t[key] for key in keys} for t in tab]

        # make table
        detail = Table(detail,names=keys)

        # rename columns
        detail.rename_columns(tuple(key_map.keys()),tuple(key_map.values()))

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
        if not any(detail['Dithers']):
            detail.remove_columns(('Dithers','Scale'))

        # remove repeats, rewinds, or loops if all are None
        for col in ('Repeat','Rewind','Loop'):
            if col in detail.colnames and (not any(detail[col])):
                detail.remove_column(col)

        # add TLOS info to footer
        if any([mode == 'C2NC2' for mode in detail['Type']]):
            # leave TLOS and TLSPN in there
            detail.meta['footer'] = ''
        else:
            tl = detail['TLOS'][0]
            span = detail['TLSPN'][0].split()[0]
            if tl == 'inf' or np.isinf(tl):
                tl = '--'
            losdet = 'TLOS = %s s @ %s deg'%(tl,span)
            detail.meta['footer'] = '%s\t\\hspace{2in}' % losdet
            detail.remove_columns(('TLSPN','TLOS'))
            

    else:
        raise NotImplementedError('Instrument %s not implemented' % instrument)

    if tex:
        # set formatter to replace '_' with '\_'
        for col in ('AORID','Name'):
            detail[col].format = COL_FORMATTER

        # set int formatter
        for col in ('NodTime','Repeat','ScanDur','ChopThrow','ChopAngle','ScanTime',
                    'ScanAmp','ScanRate','NodThrow','NodAngle','TotalTime','IntTime'):
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
                    detail[col].format = INT_FORMATTER
            except (KeyError,ValueError):
                continue


        # set units
        detail.meta['units'] = {'NodTime':'s','ChopThrow':r'$^{\prime\prime}$','ChopAngle':r'$^\circ$','ScanDur':'s','ScanAmp':r'$^{\prime\prime}$','ScanRate':'$^{\prime\prime}$/s','TotalTime':'s','NodDwell':'s','NodAngle':r'$^\circ$','NodThrow':r'$^{\prime\prime}$','IntTime':'s'}

        caption = '\\captionline{Observation details %s:}{}' % tab[0]['ObsBlk']
        detail.meta['caption'] = caption.replace('_','\_')

        # if FORCAST, split into two tables
        #### NOT IMPLEMENTED YET

        # return tex string
        detail = generate_details_tex(detail)

    else:
        # convert to list of dicts
        detail = detail.to_pandas().to_dict('records')
                
        # set int formatter
        intcols = ('NodTime','Repeat','ScanDur','ChopThrow','ChopAngle','ScanTime',
                   'ScanAmp','ScanRate','NodThrow','NodAngle','TotalTime','IntTime')
        intformat_func = partial(INT_CONVERTER,cols=intcols)
        detail = list(map(intformat_func,detail))
    return detail


def generate_details_tex(detail):
    '''Generate latex string of details table'''
    if isinstance(detail,Table):
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
        for col in d.colnames:
            newcol = '\\cellcolor{headercolor}%s'%col
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
            if ((max([len(name) for name in d['Name']]) > 13) and ('Chop Throw' in colnames)) or (('Chop Throw' in colnames) and ('Scan Amp' in colnames)):
                texcode = texcode.replace(r'\begin{tabular}','\\hspace*{-1cm}\n\\begin{tabular}')
            elif 'Nod Throw' in colnames or 'Nod Thw' in colnames:
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
        texcode = texcode.replace('Chop Throw','ChpT')
        texcode = texcode.replace('Chop Ang','ChpA')
        texcode = texcode.replace('Nod Throw','NodT')
        texcode = texcode.replace('Nod Ang','NodA')
        if 'Nod Dwell' in colnames:
            texcode = texcode.replace('Nod Dwell','Nod')
        texcodes.append(texcode)

    texcodes = '\n\n'.join(list(texcodes))
    return texcodes


def make_positions(tab, tex=True):
    '''Make position tables'''
    rows = deque()
    for t in tab:
        if t.get('RA') is None:
            continue
        aorid = t['aorID']
        name = t['POSName'] if t.get('POSName') else t['Name']
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

    return texcode

def get_pos_bundle(tab, d, odir):
    '''Download pos bundle and filter by included AORs'''
    posnames = [t['POSName'] if t.get('POSName') else t['Name'] for t in tab]
    fpid = tab.meta['Flight Plan ID']
    posfiles = ['./%s/%s_%s.pos'%(fpid,fpid,x) for x in posnames]

    postarfile = d.getPOSBundle(fpid)
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
    return table

def generate_overlay(row,nod=True,dithers=True):
    #tab = unique(tab,keys=['RA_aor','DEC_aor'])
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
        roll = float(rolls[0])*u.deg
    if isinstance(roll,(float,int,np.float,np.int)):
        roll = roll*u.deg
    elif isinstance(roll,u.Quantity):
        roll = roll
    elif isinstance(roll,str):
        roll = u.Quantity(np.float(roll),u.deg)
    else:
        roll = float(rolls[0])*u.deg

    TARFoffset = TARFOFFSET[row['InstrumentName']]

    # get band and mode for FOV
    band = row['InstrumentSpectralElement1'].split('_')[-1]
    mode = row['ObsPlanConfig']

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
                    
    width,height = [f*u.arcmin for f in fov]
    name = '%s %s' % (row['Name'],coord.to_string('hmsdms',precision=2,sep=':'))

    split = True if mode == 'TOT' and label in FOV else False

    if label not in FOV:
        # FORCAST
        label = label.replace('_TOT','')

    
    overlay = make_box(coord,width,height,roll,TARFoffset,label=label,name=name,
                       color=COLORS[row['cidx']%len(COLORS)],split=split,aorid=row['aorID'])
    
    overlay['roll'] = (float(rolls[0])*u.deg,float(rolls[1])*u.deg)

        
    if row['NodType'] != 'OTFMAP':
        # we are chop/nod dithering
        if dithers and row['DitherPattern']:
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
        
            overlay['dithers'] = [make_box(dith, width, height, angle=roll, TARFoffset=TARFoffset, label=label, split=split, color=overlay['color']) for dith in diths]

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
                                 chopangle=chopangle)

            nodAdict = row.copy()
            nodBdict = row.copy()
            ra,dec = nodA.to_string('hmsdms').split()
            nodAdict['RA'] = ra
            nodAdict['DEC'] = dec
            ra,dec = nodB.to_string('hmsdms').split()
            nodBdict['RA'] = ra
            nodBdict['DEC'] = dec
            
            overlay['nods'] = [generate_overlay(n, nod=False,dithers=False) \
                               for n in (nodAdict,nodBdict)]


    else:
        ampx,ampy = u.Quantity(row['ScanAmplitudeEL'],u.arcsec), \
                    u.Quantity(row['ScanAmplitudeXEL'],u.arcsec)
        ampx += width
        ampy += height
        overlay['dithers'] = [make_box(coord,ampx,ampy,roll,TARFoffset,label=label,name=name,
                                       color=COLORS[row['cidx']%len(COLORS)],scan=True,
                                       aorid=row['aorID'])]
        # we are scanning
        '''
        if 'Scan Amp' in row.colnames and row['Scan Amp'] is not None:
            amp = row['Scan Amp']
            if '/' in amp:
                ampx,ampy = [float(x) for x in amp.split('/')]
            else:
                ampx = float(amp)
                ampy = ampx
            ampx *= u.arcsec
            ampy *= u.arcsec

            ampx += width
            ampy += height

            ####overlay['dithers'] = [make_box(overlay['center'],width=ampx,height=ampy,angle=roll,TARFoffset=TARFoffset,label=label,color=overlay['color'])]
            ## THIS IS STILL BROKEN
            overlay['dithers'] = make_box(overlay['center'],width=ampx,height=ampy,angle=roll,TARFoffset=TARFoffset,label=label,color=overlay['color'])
        '''

    return overlay

def get_overlay_params(tab):
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
        #overlays = [generate_overlay(tab,TARFoffset,idx=idx,roll=roll) for idx,row in enumerate(tab)]
    overlays = [overlay for overlay in overlays if overlay is not None]
    return overlays

def get_recenter_image(overlaylist):
    '''Get longest wavelength "recenter"'''
    recenter = None
    #overlaylist = sorted(overlaylist,key=lambda x: x['label'])
    for overlay in overlaylist:
        if 'recenter' not in overlay:
            continue
        recenter = overlay['recenter']
        #print(overlay['label'])
        #print(recenter.to_string('hmsdms'))
    #print()

    return recenter

def make_figures(table,fdir,reg=False,guidestars=None,irsurvey=None,savefits=False,**kwargs):
    '''Generate figure from each overlay'''

    overlays = [row['overlay'] for row in table if row.get('overlay')]
    # remove any with 'nofigure' flag
    overlays = list(filter(lambda o:not o.get('nofigure',False),overlays))
    if not overlays:
        return None

    # image options are grabbed from first row in blk
    options = overlays[0].copy()

    try:
        fig,regs,hdu = get_image(overlays,**options)
    except TypeError:
        warnings.warn('Issue querying SkyView. Skipping figure.',RuntimeWarning)
        return None

    if fig is None:
        return None

    if guidestars is not None:
        pass
    '''
        guides = guidestars[table[0]['ObsBlk']]
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
                fig.show_polygons(nod['box'],edgecolor=o['color'],lw=1.5,
                                  linestyle='dotted',alpha=0.7)


    #outfile = fdir/('Leg%02d.png'%tab['Leg'][0])
    outfile = fdir/('Leg%02d.pdf'%table[0]['Leg'])
    
    fig.savefig(outfile,dpi=300)
    fig.close()

    if savefits:
        fitsdir = fdir/'images'
        fitsdir.mkdir(exist_ok=True)
        aorid = table[0]['planID']
        if not isinstance(hdu,fits.HDUList):
            hdu = [hdu]
        for h in hdu:
            surv = h.header['SURVEY'].strip().replace(' ','_')
            print(surv)
            fname = '%s_%s.fits'%(aorid,surv)
            h.writeto(fitsdir/fname,overwrite=True,output_verify='silentfix+ignore')

    if reg and regs:
        rdir = fdir/'reg'
        rdir.mkdir(exist_ok=True)
        regfile = (rdir/outfile.name).with_suffix('.reg')
        #with open(regfile,'w') as f:
        #    f.write('\n'.join(regs))
        write_ds9(regs,regfile)
        
    
    return outfile.name

def make_comments(table):
    '''Generate tex for obsblk comments'''
    comments = table[0].get('ObsBlkComment','')
    if not comments:
        return ''

    #safe convert to latex
    comments = utf8tolatex(comments)
    comments = comments.replace(r'{\textbackslash}{\textbackslash}',r'\\')
    return comments

def hawc_sio_comments(table):
    '''Assign comment block based on mode'''

    head = r'Procedure for Instrument Operator: \\'
    
    if any([row['Mode'] == 'Polarimetry' and np.isnan(row['Nod Time']) for row in table]):
        mode = 'LISPOL'
    else:
        mode = table['Mode'][0]

    comment = '%s%s' % (head, HAWC_SIO[mode])
    
    return comment
    

def match_FAORs(tables,dcs):
    aorIDs = set((row['aorID'] for table in tables for row in table))
    faors = dcs.getFAORs(aorIDs,match=True)

    for table in ProgressBar(tables):
        for row in table:
            row.update(faors.get(row['aorID']))
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

    key = table[0]['fkey']
    blk = table[0]['ObsBlkID']

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
        cmap.add_map(IMGOPTIONS.copy())
        row.update(cmap.as_dict())

    return table


def write_tex_dossier(tables, name,title,filename,
                      template='template.tex',
                      config=None,
                      mconfig=None,
                      refresh_cache=False,
                      faor=False,
                      posfiles=False,
                      reg=False,
                      dcs=None,local=None,
                      sio=False,
                      savefits=False,
                      irsurvey=None,
                      tex=True,
                      writetex=True):
    '''Write dossier pages for each table in tables'''

    # get template
    template = latex_jinja_env.get_template(template)

    # extract obsblk comments
    #comments = tables.parent.meta['ObsBlkComments']
    
    # remove tables with non obs legs
    #tables = list(filter(lambda tab:tab['ObsBlk'][0] not in ['',None], tables))

    # remove duplicate aorids
    #tables = [unique(table,keys=['Mode','AORID','Name','Band']) for table in tables]

    # get utctab just for intervals
    #fid = tables[0].meta['Flight Plan ID']
    if not dcs:
        # initialize DCS link
        dcs = DCS(refresh_cache=refresh_cache,modelcfg=mconfig)
    #utctabs = d.getFlightPlan(fid, local=local, utctab=True)
    #utctabs = list(filter(lambda x:'Leg' in x.meta, utctabs))

    # sort by order
    #orders = [table['Order'] for table in tables]


    #for table in tables:
    #    table.sort('order')

    # process config file
    if config:
        cfg = ConfigParser()
        cfg.read(config)
        # merge config with rows
        tables = list(map(lambda t:update_with_cfg(t,cfg),tables))

    if faor:
        # merge faor information
        #  locate faors
        print('Matching faors to AORIDs...')
        tables = match_FAORs(tables,dcs)

    # make overview tables
    #overviews = [make_overview(tab) for tab in tables]
    print('Generating overview...')
    over_func = partial(make_overview,tex=tex)
    overviews = ProgressBar.map(over_func,tables,multiprocess=True)

    legnames = ['Leg %i (%s)'%(table[0]['Leg'],table[0]['Name'].replace('_','\_')) for table in tables]
    #legnames = [o.meta['legName'].replace('_','\_') for o in overviews]
    #overviews = ProgressBar.map(generate_overview_tex, overviews, multiprocess=True)


    # make details tables
    #details = [make_details(tab) for tab in tables]
    print('Writing detail table...')
    # strip meta data for multiprocessing
    #metas = [table.meta.copy() for table in tables]
    #for table in tables:
    #    table.meta = None
    detail_func = partial(make_details,tex=tex)
    details = ProgressBar.map(detail_func,tables,multiprocess=True)

    # make pos tables
    print('Writing position table...')
    pos_func = partial(make_positions,tex=tex)
    positions = ProgressBar.map(pos_func,tables,multiprocess=True)

    if posfiles:
        print('Copying posfiles...')
        pdir = Path(filename).parent/'pos'
        pdir.mkdir(exist_ok=True)
        for tab in tables:
            get_pos_bundle(tab,dcs,pdir)

    # get im
    #overlays = [get_overlay_params(tab) for tab in tables]
    print('Generating overlays...')
    tables = ProgressBar.map(generate_overlays,tables,multiprocess=True)
    '''
    if guide is not None:
        guidestars = SkyCoord(ra=guide['RA'],dec=guide['DEC'],unit=(u.hourangle,u.deg))
    else:
        guidestars = None
    '''

    # get guidestars
    blkids = {row['ObsBlk'] for table in tables for row in table}
    guides = dcs.getGuideStars(list(blkids))
    guidestars = defaultdict(list)
    for guide in guides:
        coord = SkyCoord(ra=guide['RA'],dec=guide['Dec'],unit=(u.deg,u.deg))
        guide['COORD'] = coord
        guidestars[guide['ObsBlk']].append(guide)
    
    fdir = Path(filename).parent

    # define partial function for multiprocessing
    figfunc = partial(make_figures,fdir=fdir,reg=reg,guidestars=guidestars,
                      irsurvey=irsurvey,savefits=savefits)

    print('Generating figures...')
    imagefiles = ProgressBar.map(figfunc,tables,multiprocess=False)

    if sio:
        print('Generating comments...')
        comments = ProgressBar.map(hawc_sio_comments,tables,multiprocess=True)
        #comments0 = ProgressBar.map(hawc_sio_comments,tables,multiprocess=True)
        #comments = ProgressBar.map(make_comments,tables,multiprocess=True)
        #comments = ['%s\n\n%s'%(c0,c1) for c0,c1 in zip(comments0,comments)]
    else:
        print('Gathering comments...')
        comments = ProgressBar.map(make_comments,tables,multiprocess=False)

    '''
    imagefiles = deque()
    with ProgressBar(len(overlays)) as bar:
        for overlay,tab in zip(overlays,tables):
            if overlay is None or not overlay:
                imagefiles.append(None)
                continue

            if isinstance(overlay,dict):
                options = overlay.copy()
                blk = overlay['obsblk']
            else:
                options = overlay[-1].copy()
                blk = overlay[-1]['obsblk']
                options['recenter'] = get_recenter_image(overlay)

            if config and cfg.has_section(blk):
                # override options with cfg
                cfgoptions = get_cfgoptions(cfg,blk)
                options.update(**cfgoptions)

            fig = get_image(overlay,**options)

            if guide is not None:
                # show ALL guide stars from POS
                # TODO: filter out guide stars if not within footprint
                fig.show_markers(guidestars.ra,guidestars.dec,
                                 marker='o',s=80,
                                 linewidths=2)#edgecolor='#FFD700')


            # show chop/nod/dithers
            for o in overlay:
                if 'dithers' in o:
                    for dith in o['dithers']:
                        fig.show_polygons(dith['box'],edgecolor=dith['color'],lw=1,
                                          linestyle='dashed',alpha=0.7)
                if 'nods' in o:
                    for nod in o['nods']:
                        fig.show_polygons(nod['box'],edgecolor=o['color'],lw=1,
                                          linestyle='dotted',alpha=0.7)
                

            #outfile = fdir/('Leg%02d.png'%tab['Leg'][0])
            outfile = fdir/('Leg%02d.pdf'%tab['Leg'][0])
            imagefiles.append(outfile.name)

            fig.savefig(outfile,dpi=300)
            fig.close()
            bar.update()
    '''

    
    render = {'flightName':name.replace('_','\_'),
              'title':title,
              'tables':zip(legnames,overviews,
                           details,positions,imagefiles,comments)}

    pages = template.render(**render)

    if writetex:
        with open(filename,'w') as f:
            f.write(pages)

    return filename

    '''
    exit()
    with Browser('chrome',headless=True) as browser:
        a = Aladin(refresh_cache=refresh_cache,browser=browser)
        imgfiles = [a(overlay['center'],overlays=overlay,fov=0.25) if overlay else None for overlay in overlays]
    print(imgfiles)
    imgfiles = [add_text(imgfile,overlay) for imgfile,overlay in zip(imgfiles,overlays)]
    
    exit()
    '''
