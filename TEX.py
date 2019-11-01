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
from astropy.table import Table,Column,vstack,hstack
from collections import deque
from astropy.utils.console import ProgressBar
from regions import LineSkyRegion
from configparser import ConfigParser
from functools import partial
#from FAOR import FAOR
#from CMap import CMap
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

DEBUG = True

if DEBUG is False:
    log.disable_warnings_logging()
    log.setLevel('ERROR')

COL_FORMATTER = lambda x: x.replace('_','\_')
INT_FORMATTER = lambda x: '%i'%int(np.float(x)) if x not in (None,'None',-9999) else ''


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


IMGDEFAULTS =  {'width':0.2, 'height':0.2,
                'survey':'DSS2 Red',
                'vmin':None, 'vmax':None,
                'recenter':None}


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
    for k,v in IMGDEFAULTS.items():
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

def get_image(overlays,survey='DSS2 Red',width=0.2*u.deg,height=0.2*u.deg,
              reticle=False,reticle_style_kwargs=None,compass=True,
              vmin=None,vmax=None,recenter=None,invert=True,fpi=False,
              irsurvey=None,
              **kwargs):
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
        else:
            footer[k] = l[k]

    if tex:
        footer1 = '\quad '.join(['%s: %s'%(k,footer[k]) for k in mcols[0:3]])
        footer1 = ' '.join((footer1,'deg/min'))
        footer2 = '\quad '.join(['%s: %s'%(k,footer[k]) for k in mcols[-4:]])
        footer2 = ' '.join((footer2,'deg/min'))
        footer = '%s\\\\\n%s'%(footer1,footer2)
        footer = footer.replace('ROFRT','rate') \
                       .replace('THdgRT','rate') \
                       .replace('Moon','Moon ')

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
    return texcode

def make_details(tab, tex=False, faor=False):
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
                t['Band'] = '/'.join((filt,filt))
                
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
                footer = '\\\\\n'.join(footer)
            footer = 'dither\quad\quad %s\\\\'%footer
            detail.meta['footer'] = footer

    

    elif instrument == 'FORCAST':
        keys = ['ObsPlanConfig','ObsPlanMode','aorID','Name',
                'InstrumentSpectralElement1','InstrumentSpectralElement2',
                'Repeat','NodTime','TotalTime',
                'ChopThrow','ChopAngle','ChopAngleCoordinate',
                'NodThrow','NodAngle']
        key_map = {'ObsPlanConfig':'Mode','aorID':'AORID','ObsPlanMode':'Type','ChopAngleCoordinate':'Sys','InstrumentSpectralElement1':'Band','InstrumentSpectralElement2':'Slit'}
        faor_keys = ['Nod','Dithers','Scale','IntTime','FDUR','TLOS']

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

        # if all modes are NMC, drop nod params
        if all((mode == 'NMC' for mode in detail['Type'])):
            detail.remove_columns(('NodAngle','NodThrow'))

        detail.pprint()
        print(detail.colnames)
        exit()
            
                        

    else:
        raise NotImplementedError('Instrument %s not implemented' % instrument)
    print(tab)
        
    exit()
    
    dithoffset = tab['Dith Offset']
    dithunit = tab['Dith Unit'] if 'Dith Unit' in tab.colnames else tab['Num Dith']

    if sort and sort in tab.colnames:
        tab.sort(sort)

    try:
        detail = tab[['Mode','AORID','Name','POSname','Band','#Iter','Nod Time','Chop Throw','Chop Angle','Scan Dur','Scan Amp','Scan Rate','Sys']]
    except KeyError:
        # all FORCAST
        if 'Slit' in tab.colnames:
            detail = tab[['Mode','AORID','Name','POSname','Band','Slit','Nod Dwell','Chop Throw','Chop Angle','Total Exp','Obs Type','Sys']]
        else:
            detail = tab[['Mode','AORID','Name','POSname','Band','Nod Dwell','Chop Throw','Chop Angle','Total Exp','Obs Type','Sys']]

    if 'Nod Angle' in tab.colnames:
        detail.add_columns((tab['Nod Throw'],tab['Nod Angle']))
        detail.rename_column('Nod Angle','Nod Ang')


    # revert to POSname if available
    names = [row['POSname'] if row['POSname'] else row['Name'] for row in detail]
    detail['Name'] = names
    detail.remove_column('POSname')

    detail.rename_column('Chop Angle','Chop Ang')
    #detail.rename_column('RA_aor','RA')
    #detail.rename_column('DEC_aor','DEC')

    # set units
    detail.meta['units'] = {'Nod Time':'s','Chop Throw':r'$^{\prime\prime}$','Chop Ang':r'$^\circ$','Scan Dur':'s','Scan Amp':r'$^{\prime\prime}$','Scan Rate':'$^{\prime\prime}$/s','Total Exp':'s','Nod Dwell':'s','Nod Ang':r'$^\circ$','Nod Throw':r'$^{\prime\prime}$'}

    caption = '\\captionline{Observation details %s:}{}' % tab['ObsBlk'][0]
    #caption = 'Observation details %s:' % tab['ObsBlk'][0]
    detail.meta['caption'] = caption.replace('_','\_')

    # set formatter to replace '_' with '\_'
    if latex_format:
        for col in ['AORID','Name']:
            detail[col].format = COL_FORMATTER

    # Scan pol
    for row in detail:
        try:
            if row['Mode'] == 'Polarimetry' and np.isnan(row['Nod Time']):
                row['Mode'] = 'LISPOL'
                row['Nod Time'] = -9999
                row['Chop Throw'] = -9999
                row['Chop Ang'] = -9999
            elif row['Mode'] == 'Polarimetry':
                row['Scan Dur'] = -9999
                row['Scan Rate'] = -9999
        except KeyError:
            continue


    # set int formatter
    for col in ['Nod Time','#Iter','Scan Dur','Chop Throw','Chop Ang','Scan Amp','Scan Rate','Total Exp','Nod Throw','Nod Ang']:
        try:
            try:
                floats = detail[col].filled(0)
            except AttributeError:
                floats = MaskedColumn(detail[col]).filled(0)
            for i,f in enumerate(floats):
                try:
                    if np.isnan(f):
                        floats[i] = 0
                except TypeError:
                    continue
            floats = [np.float(x).is_integer() for x in floats]
            
            #if all([np.float(x).is_integer() for x in detail[col].filled(0)]):
            if all(floats):
                detail[col].format = INT_FORMATTER
        except (KeyError,ValueError):
            continue

    # if all modes are lissajous, drop nod/chop
    if all([mode in ('Lissajous','Box','LIS') for mode in detail['Mode']]):
        detail.remove_columns(['Nod Time','Chop Throw','Chop Ang','Sys'])
        if 'Nod Throw' in detail.colnames:
            detail.remove_column('Nod Throw')
        elif 'NodT' in detail.colnames:
            detail.remove_column('NodT')

        if 'Nod Ang' in detail.colnames:
            detail.remove_column('Nod Ang')
        elif 'NodA' in detail.colnames:
            detail.remove_column('NodA')

    # if all modes are polarization, drop scan
    if all([mode in ['Polarimetry','LISPOL'] for mode in detail['Mode']]):
        if all(mode == 'Polarimetry' for mode in detail['Mode']):
            detail.remove_columns(['Scan Dur','Scan Amp','Scan Rate'])
        footer = ['\t%i %s'%(int(ditho),dithu) if ditho and not np.isnan(ditho) else '--' for ditho,dithu in zip(dithoffset,dithunit)]
        #footer = ['\t%i %s'%(int(ditho),dithu) for ditho,dithu in zip(dithoffset,dithunit) if ditho and not np.isnan(ditho)]
        #footer = {b[0]:f for b,f in zip(detail['Band'],footer)}
        if len(set(footer)) == 1:
            footer = footer[0]
        else:
            footer = '\\\\\n'.join(footer)
        footer = 'dither:%s\\\\'%footer
        detail.meta['footer'] = footer
        
        '''
            if 'footer' in detail.meta:
                detail.meta['footer'] = '%s\n\t%i %s' % (detail.meta['footer'],int(dithoffset[0]),dithunit[0])
            else:
                detail.meta['footer'] = 'dither:\t%i %s' % (int(dithoffset[0]),dithunit[0])
        '''

    # if all modes are imaging, drop slit
    if all(['IMAGING' in mode for mode in detail['Mode']]):
        if 'Slit' in detail.colnames:
            detail.remove_columns('Slit')


    if latex_format:
        if '#Iter' in detail.colnames:
            detail.rename_column('#Iter','\#Iter')

    if 'Mode' in detail.colnames:
        if latex_format:
            detail['Mode'] = [x.replace('_','\_') for x in detail['Mode']]
        detail['Mode'] = [x.replace('IMAGING','IMG') for x in detail['Mode']]
        detail['Mode'] = [x.replace('GRISM','GRI') for x in detail['Mode']]
        detail['Mode'] = [x.replace('ACQUISITION','ACQ') for x in detail['Mode']]
        
    if 'Slit' in detail.colnames:
        detail['Slit'] = [x.replace('FOR_','') if x and x not in ('NONE','None') else '' for x in detail['Slit']]

    if 'Nod Dwell' in detail.colnames:
        # FORCAST mode, so remove Nod Dwell and Obs Type
        detail.remove_columns(('Nod Dwell','Obs Type'))

    if 'Total Exp' in detail.colnames:
        texp = detail['Total Exp'].copy()
        detail.remove_column('Total Exp')
        detail.add_column(texp)
        
    if 'Sys' in detail.colnames:
        sys = detail['Sys'].copy()
        detail.remove_column('Sys')
        detail.add_column(sys)

    if 'Nod Throw' in detail.colnames:
        try:
            nt = [str(x)for x in detail['Nod Throw']]
            if all([x in (None,'None','--','') for x in nt]):
                detail.remove_column('Nod Throw')
                if 'NodT' in detail.colnames:
                    detail.remove_column('NodT')
        except (ValueError,TypeError):
            pass

    if 'Nod Ang' in detail.colnames:
        try:
            nt = [str(x)for x in detail['Nod Ang']]
            if all([x in (None,'None','--','') for x in nt]):
                detail.remove_column('Nod Ang')
                if 'NodA' in detail.colnames:
                    detail.remove_column('NodA')

        except (ValueError,TypeError):
            pass


    # add pupil to band
    bands = []
    for row in detail:
        if row['Mode'] in ['POL','LISPOL','Polarimetry']:
            band = '%s/%s' % (row['Band'],row['Band'])
        elif row['Mode'] in ['Lissajous','Box','LIS']:
            band = '%s/%s' % (row['Band'],'Open')
        else:
            band = row['Band']
        bands.append(band)
    detail.replace_column('Band',Column(bands,name='Band'))
            

    # if all LISPOL, remove columns
    if all([mode == 'LISPOL' for mode in detail['Mode']]):
        for col in ('Nod Time','Chop Throw','Chop Ang','Sys'):
            if col in detail.colnames:
                detail.remove_column(col)
        detail.meta['footer'] = ''
            

    if 'FAORfile' in tab.colnames:
        # faor mode selected, so build second table with faor info
        faors = [FAOR.read(ffile,aorids=tab['_AORID']) if ffile else None for ffile in tab['FAORfile']]

        # update name from FAOR file
        detail['Name'] = [faor.preamble['TARGET'].replace('_','\_') if faor is not None else row['Name'] for row,faor in zip(detail,faors)]
        
        idcs = [faor.index(aorid) if faor is not None else None for faor,aorid in zip(faors,tab['_AORID'])]  # indices in faors matching AORIDs
        params = {}
        for k in ('NODDWELL','REPEATS','loop','rewind'):
            # get params from each run block
            param = [faor.run[idx].get(k) if idx is not None else None for idx,faor in zip(idcs,faors)]
            params[k] = param
        ptable = Table(data=params)
        ptable.rename_column('NODDWELL','Nod Dwell')
        ptable['Nod Dwell'].unit = 's'
        ptable.rename_column('REPEATS','Repeats')
        ptable.rename_column('rewind','Rewind')
        ptable.replace_column('Rewind',Column([str(x) if x not in (None,'None') else '' for x in ptable['Rewind']],name='Rewind'))

        # get comment tables
        ctables = [faor.get_comments_as_tables()[idx] if idx is not None else None for idx,faor in zip(idcs,faors)]

        ctables = list(filter(lambda x:x is not None,ctables))
        dithers = [ctable.meta.get('Dithers') if ctable else 0 for ctable in ctables]
        dithscale = [ctable.meta.get('DITHSCALE') if ctable else '' for ctable in ctables]
        dithscale = [str(d) if d else '' for d in dithscale]
        #span = ctables[0].meta.get('TLSPN') if ctables[0] else None
        #tlos = ctables[0]['TLOS'][0] if ctables[0] else None
        span = [ctable.meta.get('TLSPN') if ctable else None for ctable in ctables]
        #tlos = [ctable.meta.get('TLOS') if ctable else None for ctable in ctables]
        obsmode = [ctable.meta.get('Obs Type') if ctable else None for ctable in ctables]
        faorstart = [ctable.meta.get('Start') if ctable else None for ctable in ctables]

        try:
            ctable = vstack(ctables,metadata_conflicts='silent')
        except ValueError:
            # no faors
            return cloudpickle.dumps(detail)

        ctable.add_column(Column([s.to(u.deg).value for s in span],name='TLSPN'))
        ctable['TLSPN'].unit = u.deg
        ctable.add_column(Column(faorstart,name='CFG'),index=1)

        #ctable.remove_columns(('EFF','TLOS'))
        ctable.remove_column('EFF')
        #ctable.add_column(Column(dithers,name='Dithers'))
        ptable = hstack((ptable,ctable))
        
        ptable.add_column(Column(dithers,name='Dithers'),index=1)
        ptable.add_column(Column(obsmode,name='Type'),index=1)

        ptable.rename_column('DURPERREW','TREW')
        ptable.rename_column('INTTIME','Int Time')

        # remove TAOR
        if 'TAOR' in ptable.colnames:
            ptable.remove_column('TAOR')
        if 'TREQ' in ptable.colnames:
            ptable.remove_column('TREQ')

        # add dither offset
        if any([dith not in (0,None) for dith in dithers]):
            #ptable.add_column(Column([ditho if ditho and not (np.isnan(ditho) or ditho == 'None') else None \
            #                          for ditho in dithoffset],name='Scale'),
            #                  index=ptable.index_column('Dithers')+1)
            #ptable.replace_column('Scale',Column([str(np.int(dith)) if dith not in (None,0,'0','None') else '' for dith in ptable['Scale']],name='Scale'))
            #ptable.replace_column('Scale',Column(dithscale,name='Scale'))
            ptable.add_column(Column(dithscale,name='Scale'),index=ptable.index_column('Dithers')+1)
            ptable['Scale'].unit = u.arcsec
        
        '''
        # add footer for dithscale
        footer = ['\t%i %s'%(int(ditho),'arcsec') if ditho and not np.isnan(ditho) else '?' for ditho in dithoffset]
        if len(set(footer)) == 1:
            if footer[0] == '?':
                footer = ''
            else:
                footer = footer[0]
        else:
            footer = '\\\\\n'.join(footer)
        footer = 'dither:%s\\\\'%footer if footer else ''
        '''

        # add TLOS info to footer
        if any([mode == 'C2NC2' for mode in ptable['Type']]):
            # leave TLOS and TLSPN in there
            ptable.meta['footer'] = ''
        else:
            tl = ptable['TLOS'][0]
            if tl == 'inf' or np.isinf(tl):
                tl = '--'
            losdet = 'TLOS = %s s @ %.1f deg'%(tl,span[0].to(u.deg).value)
            ptable.meta['footer'] = '%s\t\\hspace{2in}' % losdet
            ptable.remove_columns(('TLSPN','TLOS'))

        # add units to metadata
        ptable.meta['units'] = {col:ptable[col].unit if ptable[col].unit else '' for col in ptable.colnames}

        # set int formatter
        for col in ('Nod Dwell','Int Time','FDUR','TAOR','TREW','Scale'):
            try:
                if all([np.float(x).is_integer() if x not in ['','None','--',None] else False for x in ptable[col]]):
                    ptable[col].format = INT_FORMATTER
            except KeyError:
                continue

        ptable['Int Time'].format = INT_FORMATTER

        # remove repeats or loop if all are None
        if all([x in (None,'None','--','') for x in ptable['Repeats']]):
            ptable.remove_column('Repeats')
        if all([x in (None,'None','--','') for x in ptable['Rewind']]):
            ptable.remove_column('Rewind')
        if all([x in (None,'None','--','') for x in ptable['loop']]):
            ptable.remove_column('loop')        

        aorids = [x.replace('_','\_') for x in tab['AORID']]
        ptable.add_column(Column(aorids,name='AORID'),index=0)

        # if all the same, compress
        origlen = len(ptable)
        try:
            pctable = unique(ptable,keys=('Nod Dwell','Type','Dithers','Repeats','Int Time','FDUR','TREW'))
            if len(pctable) == 1 and origlen != 1:
                pctable['AORID'][0] = '_'.join(aorids[0].split('_')[0:2]) + '_*'
                ptable = pctable
        except (ValueError,KeyError,TypeError):
            try:
                pctable = unique(ptable,keys=('Nod Dwell','Type','Dithers','Scale','loop','Int Time','FDUR','TREW','TLOS','TLSPN'))
                if len(pctable) == 1 and origlen != 1:
                    pctable['AORID'][0] = '_'.join(aorids[0].split('_')[0:2]) + '_*'
                    ptable = pctable
            except (ValueError,KeyError,TypeError):
                pass

        if 'TREW' in ptable.colnames:
            if all([x in (0,0.,'0','None','--',None) for x in ptable['TREW']]):
                ptable.remove_column('TREW')

        if 'CFG' in ptable.colnames:
            if all([x in (None,'None','') for x in ptable['CFG']]):
                ptable.remove_column('CFG')
            else:
                col = ptable['CFG'].copy()
                ptable.remove_column('CFG')
                ptable.add_column(col,index=1)
        detail = (detail,ptable)

    return cloudpickle.dumps(detail)

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


def make_positions(tab):
    '''Make position tables'''
    #position = tab[['AORID','Name','POSname','Band/Bore','RA_aor','DEC_aor']]
    position = tab[['AORID','Name','POSname','RA_aor','DEC_aor','Order']]
    names = [row['POSname'] if row['POSname'] else row['Name'] for row in position]
    position['Name'] = names
    position.remove_column('POSname')


    #position.meta['caption'] = 'Positions:'
    position.meta['caption'] = '\\captionline{Positions}{}'

    position.rename_column('RA_aor','RA')
    position.rename_column('DEC_aor','DEC')

    # set formatter to replace '_' with '\_'
    #for col in ['AORID','Name','Band/Bore']:
    for col in ['AORID','Name']:
        position[col].format = COL_FORMATTER

    # if all RA/DEC ARE blank, drop columns
    if all([ra in ['',None] for ra in position['RA']]):
        position.remove_columns(['RA','DEC'])

    # if all positions are the same, remove duplicates
    origlen = len(position)
    try:
        position = unique(position,keys=('RA','DEC'))
        position.sort(['Order','AORID'])
        if len(position) == 1 and origlen != 1:
            position['AORID'][0] = '_'.join(position['AORID'][0].split('_')[0:2]) + '_*'
    except KeyError:
        # issues with NAIF ID objects
        if 'RA' not in position.colnames:
            #position.remove_rows(slice(1,None))
            #position['AORID'][0] = '_'.join(position['AORID'][0].split('_')[0:2]) + '_*'
            return None

    #position.sort('AORID')

    # if faor, grab names from there
    if 'FAORfile' in tab.colnames:
        faors = [FAOR.read(ffile,aorids=tab['_AORID']) if ffile else None for ffile in tab['FAORfile']]
        # update name from FAOR file
        names = {aorid:faor.preamble['TARGET'].replace('_','\_') for aorid,faor in zip(tab['AORID'],faors) if faor}
        position['Name'] = [names[row['AORID']] if row['AORID'] in names else row['Name'] for row in position]

    position.replace_column('Name',Column([name.replace('_','\_') if '\_' not in name else name for name in position['Name']],name='Name'))
    position.remove_column('Order')
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
    posnames = tab['POSname']
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


def generate_overlay(tab,TARFoffset=5*u.deg,idx=0,nod=True,roll=None,dithers=True):
    #tab = unique(tab,keys=['RA_aor','DEC_aor'])
    try:
        coords = SkyCoord(ra=tab['RA_aor'],dec=tab['DEC_aor'],unit=(u.hourangle,u.deg))
    except ValueError:
        # likely a solar system object
        return None

    # we generally only want to plot the first one
    coord = coords[idx]
    row = tab[idx]

    # get roll angle
    rolls = tab.meta['Leg%i'%row['Leg']]['ROF']
    rolls = ROLL_RE.findall(rolls)[0]

    # override displayed roll if specified
    if isinstance(roll,bool):
        roll = float(rolls[0])*u.deg
    if isinstance(roll,(float,int,np.float,np.int)):
        roll = roll*u.deg
    elif isinstance(roll,u.Quantity):
        roll = roll
    elif 'FAORfile' in row.colnames:
        faor = FAOR.read(row['FAORfile'])

        if 'ROF' in faor.preamble:
            roll = u.Quantity(faor.preamble.get('ROF'))
        else:
            roll = float(rolls[0])*u.deg            
    else:
        #roll = np.mean([float(r) for r in roll])*u.deg
        roll = float(rolls[0])*u.deg


    # get band and mode for FOV
    band = row['Band']
    mode = row['Mode']
    mode = 'POL' if mode == 'Polarimetry' else 'TOT'

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
    #roll = 0*u.deg

    # override names from faorfiles
    if 'FAORfile' in row.colnames and row['FAORfile']:
        faor = FAOR.read(row['FAORfile'],aorids=tab['_AORID'])
        # update name from FAOR file
        #names = {aorid:faor.preamble['TARGET'].replace('_','\_') for aorid,faor in zip(row['AORID'],faors) if faor}
        row['POSname'] = faor.preamble['TARGET']

    name = '%s %s' % (row['POSname'],coord.to_string('hmsdms',precision=2,sep=':'))

    
    #name = name.replace('_','\_')
    #overlay = Aladin.generate_overlay(coord,width,height,angle=roll)

    split = True if mode == 'TOT' and label in FOV else False

    if label not in FOV:
        # FORCAST
        label = label.replace('_TOT','')
        ## NOTE: WE HAVE HARDCODED TARFOFFSET HERE!!!!!
        TARFoffset = u.Quantity(40.6,u.deg)
    
    overlay = make_box(coord,width,height,roll,TARFoffset,label=label,name=name,color=COLORS[idx%len(COLORS)],split=split,aorid=row['_AORID'])
    
    overlay['roll'] = (float(rolls[0])*u.deg,float(rolls[1])*u.deg)
    overlay['obsblk'] = row['ObsBlk']

    # this is simply for the np.isnan check below
    if 'Nod Dwell' in row.colnames:
        nodtime = row['Nod Dwell']
    elif 'Nod Time' in row.colnames:
        nodtime = row['Nod Time']
    else:
        nodtime = np.nan    
        
    if (row['Mode'] in ('Polarimetry','ACQUISITION','IMAGING_LWC','IMAGING_DUAL','IMAGING_SWC')) and (row['Chop Throw'] != 0) and (not np.isnan(nodtime) or ('IMAGING' in row['Mode'])):
        # we are chop/nod dithering
        if row['Sys'] == 'ERF':
            scale = row['Dith Offset']*u.arcsec
        else:
            try:
                scale = row['Dith Offset']*PIXSIZE[label]
            except KeyError:
                if 'G' in label:
                    scale = row['Dith Offset']*PIXSIZE['FORCAST_GSM']
                else:
                    scale = row['Dith Offset']*PIXSIZE['FORCAST_IMG']


        if dithers:
            diths = make_dithers(overlay['center'],scale=scale,angle=roll)
        
            overlay['dithers'] = [make_box(dith, width, height, angle=roll, TARFoffset=TARFoffset, label=label, split=split, color=overlay['color']) for dith in diths]

        if nod and ('Obs Type' in row.colnames and row['Obs Type'] == 'C2NC2'):
            chopthrow = row['Chop Throw']*u.arcsec
            chopangle = row['Chop Angle']*u.deg
            nodthrow = row['Nod Throw']*u.arcsec
            nodangle = row['Nod Angle']*u.deg

            if row['Sys'] == 'SIRF':
                chopangle += roll

            nodAchopB,nodBchopA,nodBchopB = make_C2NC2(overlay['center'],
                                                       chopthrow=chopthrow,chopangle=chopangle,
                                                       nodthrow=nodthrow,nodangle=nodangle)
            nodAchopBtab = Table(rows=[tab[idx]],names=tab.colnames,meta=tab.meta)
            nodBchopAtab = Table(rows=[tab[idx]],names=tab.colnames,meta=tab.meta)
            nodBchopBtab = Table(rows=[tab[idx]],names=tab.colnames,meta=tab.meta)

            for ntab,n in zip((nodAchopBtab,nodBchopAtab,nodBchopBtab),(nodAchopB,nodBchopA,nodBchopB)):
                ntab.replace_column('RA_aor',[n.to_string('hmsdms').split()[0]])
                ntab.replace_column('DEC_aor',[n.to_string('hmsdms').split()[1]])

            overlay['nods'] = [generate_overlay(nod,TARFoffset=TARFoffset,
                                                idx=0,nod=False,roll=roll,dithers=False) \
                               for nod in (nodAchopBtab,nodBchopAtab,nodBchopBtab)]

        elif nod:
            chopthrow = row['Chop Throw']*u.arcsec
            chopangle = row['Chop Angle']*u.deg

            if row['Sys'] == 'SIRF':
                chopangle += roll
                
            nodA,nodB = make_NMC(overlay['center'],
                                 chopthrow=chopthrow,
                                 #chopangle=roll+chopangle)
                                 chopangle=chopangle)

            nodtabA = Table(rows=[tab[idx]],names=tab.colnames,meta=tab.meta)
            nodtabB = Table(rows=[tab[idx]],names=tab.colnames,meta=tab.meta)

            nodtabA.replace_column('RA_aor',[nodA.to_string('hmsdms').split()[0]])
            nodtabA.replace_column('DEC_aor',[nodA.to_string('hmsdms').split()[1]])
            nodtabB.replace_column('RA_aor',[nodB.to_string('hmsdms').split()[0]])
            nodtabB.replace_column('DEC_aor',[[nodB.to_string('hmsdms').split()[1]]])
            
            overlay['nods'] = [generate_overlay(nod,TARFoffset=TARFoffset,
                                                idx=0,nod=False,roll=roll,dithers=False) \
                               for nod in (nodtabA,nodtabB)]


    else:
        pass
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

def get_overlay_params(tab,TARFoffset=5*u.deg):
    overlays = deque()
    for idx,row in enumerate(tab):
        if 'IMGOVERRIDES' in tab.meta:
            key = 'Leg%02d__%s'%(row['Leg'],row['_AORID'])
            if key in tab.meta['IMGOVERRIDES']:
                roll = u.Quantity(tab.meta['IMGOVERRIDES'][key],u.deg)
            elif tab.meta['IMGOVERRIDES'].get('roll',False):
                roll = u.Quantity(row['Angle'],u.deg)
            else:
                roll = None
        else:
            roll = None
        overlay = generate_overlay(tab,TARFoffset,idx=idx,roll=roll)
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

def make_figures(tup,cfg,guidestars,fdir,reg=False,irsurvey=None,savefits=False,**kwargs):
    '''Generate figure from each overlay'''
    overlay,tab = tup
    if overlay is None or not overlay:
        return None
        
    if isinstance(overlay,dict):
        options = overlay.copy()
        blk = overlay['obsblk']
    else:
        options = overlay[-1].copy()
        blk = overlay[-1]['obsblk']
        options['recenter'] = get_recenter_image(overlay)

    if cfg is not None and cfg.has_section(blk):
        # override options with cfg
        cfgoptions = get_cfgoptions(cfg,blk)
        options.update(**cfgoptions)

    # override options with tab metadata
    if 'IMGOVERRIDES' in tab.meta:
        options.update(**tab.meta['IMGOVERRIDES'])

    if 'nofigure' in options and options['nofigure']:
        return None

    # add irsurvey
    if irsurvey is not None:
        options['irsurvey'] = irsurvey

    try:
        fig,regs,hdu = get_image(overlay,**options)
    except TypeError:
        warnings.warn('Issue querying SkyView. Skipping figure.',RuntimeWarning)
        return None

    if fig is None:
        return None


    if guidestars is not None:
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
                fig.show_polygons(nod['box'],edgecolor=o['color'],lw=1.5,
                                  linestyle='dotted',alpha=0.7)


    #outfile = fdir/('Leg%02d.png'%tab['Leg'][0])
    outfile = fdir/('Leg%02d.pdf'%tab['Leg'][0])
    
    fig.savefig(outfile,dpi=300)
    fig.close()

    if savefits:
        fitsdir = fdir/'images'
        fitsdir.mkdir(exist_ok=True)
        aorid = tab['AOR'][0]
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
    if 'ObsBlkComment' not in table.colnames:
        return ''
    comments = [c for c in table['ObsBlkComment'].filled('') if c not in ('','--')]

    comments = np.unique(comments)
    comments = '\n\n'.join(comments)

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

def write_tex_dossier(tables, name,title,filename,
                      template='template.tex',
                      config=None,
                      mconfig=None,
                      refresh_cache=False,
                      guide=None,
                      faor=False,
                      posfiles=False,
                      reg=False,
                      #TARFoffset=40.6*u.deg,
                      TARFoffset=3*u.deg,
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

        # attach options to rows, and remove rows if specified
        for idx,table in enumerate(tables):
            # if leg and program has aorids, then keep only the matching rows
            key = 'Leg%02d_%s'%(table['Leg'][0],table.meta['Flight Plan ID'])
            blk = table['ObsBlk'][0]
            if cfg.has_section(key) and 'aorids' in cfg[key]:
                keep = [x.strip() for x in cfg[key].get('aorids').split(',')]
                idcs = [jdx for jdx,row in enumerate(table) if not ((row['_AORID'] in keep) or (row['AORID'] in keep))]
                if len(idcs) == len(table):
                    warnings.warn('All aorids would be removed from this AOR. Check cfg file.',RuntimeWarning)
                    continue
                table.remove_rows(idcs)
                tables[idx] = table
            elif cfg.has_section(blk) and cfg.has_option(blk,'aorids'):
                keep = [x.strip() for x in cfg[blk].get('aorids').split(',')]
                idcs = [idx for idx,row in enumerate(table) if not ((row['_AORID'] in keep) or (row['AORID'] in keep))]
                table.remove_rows(idcs)
                tables[idx] = table
            else:
                pass

            # get other options specified
            if cfg.has_section(blk):
                table.meta['IMGOVERRIDES'] = get_cfgoptions(cfg,blk)
            if cfg.has_section(key):
                if 'IMGOVERRIDES' in table.meta:
                    table.meta['IMGOVERRIDES'].update(**get_cfgoptions(cfg,key))
                else:
                    table.meta['IMGOVERRIDES'] = get_cfgoptions(cfg,key)

            faordirs = {}
            imgoverrides = {}
            intervals = {}
            for row in table:
                # create chainmap for each aorid
                aorid = row['AORID']
                planid = '_'.join(aorid.split('_')[0:2])
                aoridmod = row['_AORID']
                miss = table.meta['Flight Plan ID']
                legmiss = 'Leg%02d_%s'%(row['Leg'],miss)
                mission = '_'.join((miss,aorid))
                missionmod = '_'.join((miss,aoridmod))
                missionplan = '_'.join((miss,planid))
                legsec = 'Leg%02d_%s'%(row['Leg'],mission)
                legsecmod = 'Leg%02d_%s'%(row['Leg'],missionmod)
                legsecplan = 'Leg%02d_%s'%(row['Leg'],missionplan)
                aorsec = (legsec,legsecmod,legsecplan,legmiss,mission,missionmod,missionplan,blk,aorid,aoridmod,planid,'DEFAULT')
                # chainmap checks each section for values, starting with aorid
                chmap = CMap(cfg,aorsec)
                #row.meta['CFG'] = chmap
                #row.meta['faordir'] = chmap.get('faordir')
                #print(row.meta['faordir'])
                key = 'Leg%02d__%s'%(row['Leg'],row['_AORID'])
                faordirs[key] = chmap.get('faordir')

                if 'interval' in chmap:
                    inter = chmap.getfloat('interval')
                    # get split_leg table
                    utctab = list(filter(lambda x:x.meta['Leg'] == row['Leg'], utctabs))[0]
                    iTab = split_leg(utctab, inter)
                    intervals[key] = iTab
                    
                if 'roll' in chmap:
                    imgoverrides[key] = chmap.get('roll')
                
            table.meta['FAORDIRS'] = faordirs
            table.meta['INTERVALS'] = intervals
            #print(intervals)
            if 'IMGOVERRIDES' in table.meta:
                table.meta['IMGOVERRIDES'].update(**imgoverrides)
            else:
                table.meta['IMGOVERRIDES'] = imgoverrides
                
    else:
        cfg = None


    if faor:
        # merge faor information
        #  locate faors
        print('Matching faors to AORIDs...')
        tables = match_FAORs(tables,dcs)

    # remove CFG from rows
    '''
    # unecessary with __deecopy__ implemented
    for table in tables:
        for row in table:
            if 'CFG' in row.meta:
                del row.meta['CFG']
    '''    
    # make overview tables
    #overviews = [make_overview(tab) for tab in tables]
    print('Generating overview...')
    over_func = partial(make_overview,tex=tex)
    overviews = ProgressBar.map(over_func,tables,multiprocess=True)
    #legnames = [o.meta['legName'].replace('_','\_') for o in overviews]
    #overviews = ProgressBar.map(generate_overview_tex, overviews, multiprocess=True)


    # make details tables
    #details = [make_details(tab) for tab in tables]
    print('Writing detail table...')
    # strip meta data for multiprocessing
    #metas = [table.meta.copy() for table in tables]
    #for table in tables:
    #    table.meta = None

    mp = not DEBUG
    details = ProgressBar.map(make_details,tables,multiprocess=False)

    exit()
    details = [pickle.loads(detail) for detail in details]
    #for table,meta in zip(tables,metas):
    #    table.meta = meta
    
    details = [generate_details_tex(detail) for detail in details]

    # make pos tables
    positions = [make_positions(tab) for tab in tables]
    #positions = [generate_pos_tex(pos) for pos in positions]
    print('Writing position table...')
    positions = ProgressBar.map(generate_pos_tex,positions,multiprocess=False)

    if posfiles:
        print('Copying posfiles...')
        pdir = Path(filename).parent/'pos'
        pdir.mkdir(exist_ok=True)
        for tab in tables:
            get_pos_bundle(tab,d,pdir)

    # get im
    #overlays = [get_overlay_params(tab) for tab in tables]
    print('Generating overlays...')
    overlay_func = partial(get_overlay_params,TARFoffset=TARFoffset)
    overlays = ProgressBar.map(overlay_func,tables,multiprocess=True)
    #overlays = [pickle.loads(overlay) for overlay in overlays]

    #print(overlays)
    #print(overlays[0]['reg'])
    #exit()

    #images = [get_image(overlay) for overlay in overlays]

    if guide is not None:
        guidestars = SkyCoord(ra=guide['RA'],dec=guide['DEC'],unit=(u.hourangle,u.deg))
    else:
        guidestars = None

    fdir = Path(filename).parent

    # define partial function for multiprocessing
    figfunc = partial(make_figures,cfg=cfg,guidestars=guidestars,fdir=fdir,reg=reg,
                      irsurvey=irsurvey,savefits=savefits)

    print('Generating figures...')
    imagefiles = ProgressBar.map(figfunc, list(zip(overlays,tables)),
                                 multiprocess=True)

    if sio:
        print('Generating comments...')
        comments = ProgressBar.map(hawc_sio_comments,tables,multiprocess=True)
        #comments0 = ProgressBar.map(hawc_sio_comments,tables,multiprocess=True)
        #comments = ProgressBar.map(make_comments,tables,multiprocess=True)
        #comments = ['%s\n\n%s'%(c0,c1) for c0,c1 in zip(comments0,comments)]
    else:
        print('Gathering comments...')
        comments = ProgressBar.map(make_comments,tables,multiprocess=True)

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
