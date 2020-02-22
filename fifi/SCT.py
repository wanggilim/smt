#! /usr/bin/env python
from configparser import ConfigParser
import argparse
from pathlib import Path
#from dcs import DCS
import json
from astropy.utils.console import ProgressBar
from functools import partial
import numpy as np
from obsmaker.grating import wavelength2inductosyn, inductosyn2wavelength
from contextlib import redirect_stdout
from io import StringIO

sct_DEFAULT = str(Path(__file__).parent.resolve()/'obsmaker.cfg')
mcfg_DEFAULT = str(Path(__file__).resolve().parent.parent/'dcs/DBmodels.cfg')

BADCHAR = ('/','\\',':','*',' ')

def replaceBadChar(string):
    """ replace some reserved characters with '_'
    """
    for bad in BADCHAR:
        string = string.replace(bad,'_')
    return string

def proc_row(aor, values, cfg):

    # all values are primitives, so this should be fine
    values = values.copy()
    # for Moving targets, set lat and lon to 0, equinox to J2000
    #if aor['lat'][0] == '':
    #    aor['lat'][0] = '0.0'
    #    aor['lon'][0] = '0.0'
    #    aor['equinoxDesc'][0] = 'J2000'

    values['TARGET_NAME'] = aor['target'].replace(" ", "_")
    values['AORID'] = aor['aorID']
    values['OBSID'] = '%s_%s' % (values['TARGET_NAME'].replace('@', ''),
                                 replaceBadChar(aor['title']))
    values['SRCTYPE'] = aor['SourceType'].upper()
    values['INSTMODE'] = aor['ObsPlanMode']

    if aor['RA'] is None:
        # likely a solar system object
        values['TARGET_LAMBDA'] = '0 0 0'
        values['TARGET_BETA'] = '0 0 0'
    else:
        # I hate this, but I can't figure out Dario's formatting
        ra = aor['RA'].split(':')
        sra = ' '.join((str(int(ra[0])),str(int(ra[1])),str(round(float(ra[2]), 2))))
        values['TARGET_LAMBDA'] = sra

        dec = aor['DEC'].split(':')
        sdec = '+' * (int(dec[0]) >= 0) + str(int(dec[0])) + ' ' + str(int(dec[1])) + ' ' + \
               str(round(float(dec[2]), 1))
        values['TARGET_BETA'] = sdec

    # from Dario's code.  no clue
    detang = ((aor['MapRotationAngle'] + 180 + 360) % 360) - 180
    values['DETANGLE'] = detang    # in degrees
    detang *= np.pi/180.0
    cosa = np.cos(detang)
    sina = np.sin(detang)
    r = np.array([[cosa, -sina], [sina, cosa]])

    dra = json.loads(aor['deltaRaV'])
    ddec = json.loads(aor['deltaDecW'])

    if len(dra) == 1:
        mapoffsets = np.array(dra + ddec)
    else:
        mapoffsets = np.array((dra,ddec))
    rot_mapoffsets = np.transpose(np.dot(np.transpose(r), mapoffsets))

    if len(dra) == 1:
        values['DITHMAP_NUMPOINTS'] = 1
    else:
        values['DITHMAP_NUMPOINTS'] = len(rot_mapoffsets)

    values['CHOPCOORD_SYSTEM'] = aor['ChopAngleCoordinate']
    values['CHOP_AMP'] = aor['ChopThrow']/2
    values['CHOP_POSANG'] = (aor['ChopAngle'] + 270 + 360) % 360
    # CCW from N in USpot, S of E in ObsMaker

    if values['CHOPCOORD_SYSTEM'] == 'HORIZON':
        values['CHOP_POSANG'] = 0

    if aor['ChopType'] == 'Sym':
        values['TRACKING'] = 'On'
        #values['OBSMODE'] = 'Beam switching'
        values['OBSMODE'] = 'Symmetric'
        values['OFFPOS'] = 'Matched'
        values['OFFPOS_LAMBDA'] = '0.0'
        values['OFFPOS_BETA'] = '0.0'
    elif aor['ChopType'] == 'Asym':
        values['TRACKING'] = 'Off'
        #values['OBSMODE'] = 'Asym chop + off'
        values['OBSMODE'] = 'Asymmetric'
        if aor['ReferenceType'] == 'RA_Dec':
            values['OFFPOS_LAMBDA'] = aor['RefRA']
            values['OFFPOS_BETA']   = aor['RefDec']
            values['OFFPOS'] = 'Absolute'
        elif aor['ReferenceType'] == 'Offset':
            values['OFFPOS_LAMBDA'] = aor['RAOffset']
            values['OFFPOS_BETA']   = aor['DecOffset']
            if aor['MapRefPos'] == 'true':
                values['OFFPOS'] = 'Relative to active pos'
            elif aor['MapRefPos'] == 'false':
                values['OFFPOS'] = 'Relative to target'

    values['PRIMARYARRAY'] = aor['PrimeArray'].upper()
    values['DICHROIC'] = aor['Dichroic'][0:3]  # 105 or 130

    #### updated for Cycle 4
    #  blue rest wavelength and species
    blue_lam = aor['WavelengthBlue']  # Wavelength in Cycle 3
    values['BLUE_LINE'] = 'Custom'
    values['BLUE_MICRON'] = blue_lam  # user-entered wavelength
    # red rest wavelength and species
    red_lam = aor['WavelengthRed']   # Wavelength2 in Cycle 3
    values['RED_LINE'] = 'Custom'
    values['RED_MICRON'] = red_lam  # user-entered wavelength


    # Cycle 3: offset in km/s
    # values['BLUE_OFFSET'] = diff[line_idx] / \
    #    obs_ref_blue_lambdas[idx_blue[line_idx]] * speed_of_light
    # Cycle 4: offset in either kmPerSec or z
    # ObsMaker line offset must be in kms or um
    # convert z to um: offset = obs_um - rest_um = z * um_rest
    if aor['RedshiftUnit'] == 'z':
        values['BLUE_OFFSET'] = aor['Redshift'] * blue_lam
        values['BLUE_OFFSET_TYPE'] = 'um'
        values['RED_OFFSET'] = aor['Redshift'] * red_lam
        values['RED_OFFSET_TYPE'] = 'um'
        values['REDSHIFT'] = aor['Redshift']
    else:  # other option is 'kmPerSec' (Cycle4), '' (Cycle5; kmPerSec implied)
        values['BLUE_OFFSET'] = aor['Redshift']
        values['BLUE_OFFSET_TYPE'] = 'kms'
        values['RED_OFFSET'] = aor['Redshift']
        values['RED_OFFSET_TYPE'] = 'kms'
        values['REDSHIFT'] = aor['Redshift']/299792.458  # (v/c)


    # File Group IDs for DPS
    # Target_wavelength - USpot allows 3 significant digits
    dot = str(blue_lam).find('.')
    values['FILEGP_B'] = values['TARGET_NAME'].replace('@', '') + '_' + \
        str(blue_lam)[: dot + 4]
    dot = str(red_lam).find('.')
    values['FILEGP_R'] = values['TARGET_NAME'].replace('@', '') + '_' +  \
        str(red_lam)[: dot + 4]

    #Order filter
    if blue_lam < 71:   # Wavelength in Cycle 3
        values['ORDER'] = '2'
        #blue_um_per_pix = poly(blue_lam, config['blue2_coef'])
        bluecoef = json.loads(cfg['CONSTANTS']['blue2_coef'])
        blue_um_per_pix = np.polyval(np.flip(bluecoef), blue_lam)
        #blue_um_per_pix = np.polyval(config['blue2_coef'], blue_lam)
    else:
        values['ORDER'] = '1'
        bluecoef = json.loads(cfg['CONSTANTS']['blue1_coef'])
        blue_um_per_pix = np.polyval(np.flip(bluecoef), blue_lam)
        #blue_um_per_pix = np.polyval(config['blue1_coef'], blue_lam)
    redcoef = json.loads(cfg['CONSTANTS']['red_coef'])
    red_um_per_pix = np.polyval(np.flip(redcoef), red_lam)

    dichroic = int(values['DICHROIC'])
    gratpos = wavelength2inductosyn(blue_lam, dichroic, 'BLUE', values['ORDER'], obsdate='')
    result, result_dwdp = inductosyn2wavelength(gratpos=gratpos, order=values['ORDER'],array='BLUE',
                                                dichroic=dichroic, obsdate='')
    blue_um_per_pix = np.mean(result_dwdp, axis=(1,2))

    gratpos = wavelength2inductosyn(blue_lam, dichroic, 'RED', values['ORDER'], obsdate='')
    result, result_dwdp = inductosyn2wavelength(gratpos=gratpos, order=values['ORDER'], array='RED',
                   dichroic=dichroic, obsdate='')
    red_um_per_pix = np.mean(result_dwdp, axis=(1,2))
    
    values['BLUE_FILTER'] = values['ORDER']

    nodcycles = aor['Repeat']
    values['NODCYCLES'] = nodcycles

    # Bandwidths  - updated for Cycle 4
    # if NodType or ObsPlanMode is SPECTRAL_SCAN, BandwidthBlue/Red is in
    # micron; other Modes are in kmPerSec.
    if aor['NodType'] == 'SPECTRAL_SCAN':
        bandwidthBlue_pix = max((aor['BandwidthBlue']/blue_um_per_pix, 6.))
        bandwidthRed_pix = max((aor['BandwidthRed']/red_um_per_pix, 6.))
    else:
        bandwidthBlue_pix = max((aor['BandwidthBlue'] * blue_lam / 299792.458 / blue_um_per_pix, 6.))
        bandwidthRed_pix = max((aor['BandwidthRed'] * red_lam / 299792.458 / red_um_per_pix, 6.))

    blue_pix_per_nod = bandwidthBlue_pix / nodcycles
    red_pix_per_nod = bandwidthRed_pix / nodcycles

    values['BLUE_POSUP'] = int(np.ceil(blue_pix_per_nod / cfg['CONSTANTS'].getfloat('max_stepsize_inPix')) * nodcycles)
    values['BLUE_SIZEUP'] = 0.5
    # values['BLUE_SIZEUP'] = bandwidthBlue_pix / values['BLUE_POSUP']
    values['RED_POSUP'] = int(np.ceil(red_pix_per_nod / cfg['CONSTANTS'].getfloat('max_stepsize_inPix')) * nodcycles)
    values['RED_SIZEUP'] = 0.5

    # Cycle 5: SCANDIST is always Up, SPLITS is always 1,
    # RED_LAMBDA and BLUE_LAMBDA are always Inward dither
    values['SCANDIST'] = 'Up'
    values['RED_LAMBDA'] = 'Inward dither'
    values['BLUE_LAMBDA'] = 'Inward dither'
    values['SPLITS'] = 1

    values['TIME_POINT'] = aor['TimePerPoint']
    chopCycles_per_nod = 2 * aor['TimePerPoint']
    values['BLUE_CHOPCYC'] = int(np.ceil(chopCycles_per_nod * nodcycles / \
                                           (values['BLUE_POSUP'] * values['SPLITS'])))
    values['RED_CHOPCYC'] = int(np.ceil(chopCycles_per_nod * nodcycles / \
                                        (values['RED_POSUP'] * values['SPLITS'])))

    skey = '_'.join((values['OBSID'],values['AORID']))
    #fname = '%s.sct' % skey
    mname = '%s.map' % skey

    values['MAPLISTPATH'] = mname

    # Add extra keywords
    values['PROPID'] = aor['planID']
    values['OBSERVER'] = aor['PI']

    # THESE NEED TO BE POPPED LATER
    #values['sctfile'] = fname
    values['rot_mapoffsets'] = json.dumps(rot_mapoffsets.tolist())
    values['skey'] = skey
    
    return values


def AOR_to_SCTdicts(aors, cfg=sct_DEFAULT, mcfg=mcfg_DEFAULT, mp=True):
    # SCT keys with empty values
    values = dict.fromkeys(json.loads(mcfg['SCT']['data_keys']))

    # update with defaults
    values.update(json.loads(cfg['DEFAULTS']['defaults']))

    proc_func = partial(proc_row,values=values,cfg=cfg)

    # suppress bar
    with redirect_stdout(StringIO()):
        scts = ProgressBar.map(proc_func,aors,multiprocess=mp)

    return scts


def write_files(sct, odir):

    mapfile = Path(odir)/sct['MAPLISTPATH']
    sctfile = Path(odir)/('%s.sct'%sct['skey'])

    #write map file
    rot_mapoffsets = json.loads(sct['rot_mapoffsets'])
    if isinstance(rot_mapoffsets[0],float):
        rot_mapoffsets = [rot_mapoffsets]

    with open(mapfile,'w') as f:
        f.write("%s%s" % (str(sct['TARGET_LAMBDA']).rjust(12),
                          str(sct['TARGET_BETA']).rjust(12) + "\n"))

        for line in rot_mapoffsets:
            rot = str(round(line[0], 4)).rjust(12) + \
                  str(round(line[1], 4)).rjust(
                      max([12, len(str(round(line[1], 4))) + 2]))
            f.write(rot + "\n")
        """
        if len(rot_mapoffsets[0]) > 1:
            for line in rot_mapoffsets:
                rot = str(round(line[0], 4)).rjust(12) + \
                      str(round(line[1], 4)).rjust(
                        max([12, len(str(round(line[1], 4))) + 2]))
                f.write(rot + "\n")
        else:
            # rot = str(round(rot_mapoffsets[0][0], 4)).rjust(12) +  \
            #       str(round(rot_mapoffsets[1][0], 4)).rjust(12)
            rot = str(rot_mapoffsets[0][0]).rjust(12) + \
                  str(rot_mapoffsets[1][0]).rjust(
                    max([12, len(str(rot_mapoffsets[1][0])) + 2]))
            f.write(rot)
        """

    
    #write sct file
    with open(sctfile,'w') as f:
        for k,v in sct.items():
            
            if k in ('skey','rot_mapoffsets','FILENAME','TIMESTAMP'):
                # skip these
                continue
            
            f.write("%s#%s\n" %
            (str(v).ljust(max([25, len(str(v)) + 2])), k))

    return sctfile,mapfile

def read_sctfile(filename,mdir=None):
    with open(filename,'r') as f:
        lines = f.readlines()
    tup = (line.split('#') for line in lines)
    params = {v[1].strip():v[0].strip() for v in tup}

    if mdir:
        m = Path(mdir)/params['MAPLISTPATH']
    else:
        m = Path(filename).parent/params['MAPLISTPATH']

    _,rot_mapoffsets = read_mapfile(m)
    params['rot_mapoffsets'] = rot_mapoffsets
    params['skey'] = '_'.join((params['OBSID'],params['AORID']))
    params['FILENAME'] = str(filename)

    return params


def read_mapfile(filename):
    with open(filename,'r') as f:
        lines = f.readlines()

    coord = lines[0].strip()
    lines = (line.split() for line in lines[1:])
    lines = [[float(line[0]),float(line[1])] for line in lines]
    
    return coord,json.dumps(lines)
    


def main():
    parser = argparse.ArgumentParser(description='Convert AOR to SCT')
    parser.add_argument('aorfile',type=str,help='AOR file')
    parser.add_argument('-o','--outdir',dest='outdir',
                        type=str,help='Output directory (default is aorfile dir)')
    parser.add_argument('-cfg',type=str,default=sct_DEFAULT,help='cfg file (default=/library/path/obsmaker.cfg)')
    parser.add_argument('-mcfg',type=str,default=mcfg_DEFAULT,help='Model config file (default=/library/path/DBmodels.cfg)')

    args = parser.parse_args()

    cfg = ConfigParser()
    cfg.read(args.cfg)

    mcfg = ConfigParser()
    mcfg.read(args.mcfg)

    from dcs.DBmodels import AOR_to_rows
    aors = AOR_to_rows(args.aorfile, mcfg['AOR'], convert_dtype=True)
    
    scts = AOR_to_SCTdicts(aors, cfg, mcfg)

    odir = args.outdir if args.outdir else Path(args.aorfile).parent.resolve()

    write_func = partial(write_files,odir=odir)

    print('Writing files...')
    res = ProgressBar.map(write_func,scts,multiprocess=False)
    print('SCT and MAP files written:')
    for sfile,mfile in res:
        print(sfile)
        print(mfile)
    

if __name__ == '__main__':
    main()
    #fname = '202002_FI/LANCELOT/NGC4258P3N_NGC4258P3N_07_0154_4.sct'
    #print(read_sctfile(fname))
