#! /usr/bin/env python
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import download_file, clear_download_cache, is_url_in_cache
from astropy.config import set_temp_cache,get_cache_dir
import urllib.parse
from pathlib import Path
from bs4 import BeautifulSoup
import re

IRSA_URL = 'https://irsa.ipac.caltech.edu/cgi-bin/MSX/nph-msx'
WORKSPACE_URL = 'https://irsa.ipac.caltech.edu'
MSX_BAND_MAP = {'A':'0','C':'1','D':'2','E':'3'}

FITS_RE = re.compile('\/workspace.*\.fits')

cachedir = Path(get_cache_dir())/'DCS'

def query_region(coordinate,size=1*u.deg,band='E',show_progress=True):
    if isinstance(coordinate,SkyCoord):
        coordinate = coordinate.to_string('hmsdms',sep=':')
        
    data = {'objstr':coordinate,'size':size.to(u.deg).value,
            'band':MSX_BAND_MAP[band],'submit':'submit'}

    data_str = urllib.parse.urlencode(data)
    url = '%s?%s' % (IRSA_URL, data_str)
    
    with set_temp_cache(cachedir):
        cfile = download_file(url,show_progress=show_progress,cache=True)

    with open(cfile,'r') as f:
        #soup = BeautifulSoup(f.read(),'lxml')
        match = FITS_RE.findall(f.read())

    if not match:
        return None

    durl = '%s%s' % (WORKSPACE_URL, match[0])
    with set_temp_cache(cachedir):
        dfile = download_file(durl,show_progress=show_progress,cache=True)
    return dfile



if __name__ == '__main__':
    coord = SkyCoord('20:39:02.9993 +42:25:30.011',unit=(u.hourangle,u.deg))
    query_region(coord)
    
