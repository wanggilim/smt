#! /usr/bin/env python
from bs4 import BeautifulSoup
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from astropy.io import registry
from collections import OrderedDict
from AOR import is_integer, get_first_tag, get_all_tags

def get_preamble(elem):
    '''Returns preamble properties (e.g. ProposalID, Title, PI)'''

    tags = ['ObsPlanID', 'ProposalTitle', 'Category','ScienceKeywords']

    preamble = {tag:get_first_tag(elem,tag,as_string=True) for tag in tags}
    investigators = get_all_tags(elem,'Investigator')
    #preamble['Investigators'] = investigators

    # PI
    primary = investigators[0].identity
    preamble['Investigator'] = primary
    
    return preamble

def get_obs_mode(obs):
    '''Get observation mode from request'''
    if obs.config.string == 'POLARIZATION':
        return 'Polarimetry'
    if obs.config.string == 'TOTAL_INTENSITY':
        return 'Total Intensity'
    '''
    if get_first_tag(obs,'ScanType') == 'Box':
        return 'Box'
    '''
    return None

def get_filter(elem):
    '''Get filters and boresights from elem'''
    filt = elem.spectralelement.string.split('_')[1]

    bore = '%s_main' % filt.lower()
    if filt in ['B','C']:
        # do we even use B?
        bore = 'bc_main'
    return filt,bore


def proc_observation(obs):
    '''Pull info from observation. Return dict'''
    ra =  float(obs.ra.string)
    dec = float(obs.dec.string)
    try:
        equ = obs.equinox.strip()
    except (AttributeError,TypeError):
        equ = 'J2000'
    coord = SkyCoord(ra=ra,dec=dec,equinox=equ,unit=(u.hourangle,u.deg))

    # convert to hexadecimal
    ra,dec = coord.to_string('hmsdms',sep=':',precision=2).split()
    name = obs.astroobjectname.string.replace('Name','').strip()
    aorid = obs.aorid.string

    # get mode/filter
    mode = get_obs_mode(obs)
    filt,bore = get_filter(obs)

    params = OrderedDict(zip(('AORID','Name','RA','DEC'),(aorid,name,ra,dec)))
    params['Mode'] = mode
    params['Filter'] = filt
    params['Boresight'] = bore
    return params

def get_observations(elem,
                     fields=['AORID','Mode','Name','RA','DEC','Filter','Boresight']):
    '''Returns observation table'''
    observations = get_all_tags(elem,'ProposedObservation')

    params = [proc_observation(obs) for obs in observations]
    
    observations = Table(params)
    observations = observations[fields]
    
    return observations


def read_OBS_file(filename):
    '''Parse OBS file into table'''

    # open xml
    with open(filename,'r') as f:
        obs = BeautifulSoup(f,'lxml')

    preamble = get_preamble(obs)

    observations = get_observations(obs)
    
    observations.meta.update(**preamble)

    return observations

# Register Table class reader
registry.register_reader('obs-tab', Table, read_OBS_file)


if __name__ == '__main__':
    fname = '/home/gordon/.astropy/cache/DCS/astropy/download/py3/50d85f94b22a8d542ed33027e23257ad'

    obstab = Table.read(fname,format='obs-tab')
    obstab.pprint()
