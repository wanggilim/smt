#! /usr/bin/env python
from bs4 import BeautifulSoup,NavigableString
from astropy.table import Table,Column
from pandas import read_html
from astropy.io import registry
import numpy as np
import astropy.units as u
from collections import deque


def read_obssearch_page(filename):
    '''Read AOR search page results'''
    with open(filename,'r') as f:
        soup = BeautifulSoup(f.read(),'html.parser')


    ths = soup.find_all('th')
    # this table should be unique
    try:
        htable = list(filter(lambda x: 'AOR ID [Target]' in x.text,ths))[0]
    except IndexError:
        return None

    htable = htable.parent.parent
    table = read_html(str(htable))
    table = Table.from_pandas(table[0])

    # get rid of 'View AORs'
    obsblks = Column([x.split('View')[0].strip() for x in table['ObsBlock ID']],name='ObsBlock ID')
    table.replace_column('ObsBlock ID',obsblks)

    # split targets and aor ids
    lines = [x.split('  ') for x in table['AOR ID [Target]']]
    aorids = deque()
    targets = deque()
    for line in lines:
        ais = [x.split()[0].strip() for x in line]
        ts = [' '.join(x.split()[1:]).strip() for x in line]
        ts = [target[1:-1] for target in ts if target[0] == '[' and target[-1] == ']']
        aorids.append(';'.join(ais))
        targets.append(';'.join(ts))
        
    aorids = Column(aorids,name='AOR ID')
    targets = Column(targets,name='Target')
    table.remove_column('AOR ID [Target]')
    table.add_columns((aorids,targets),indexes=(2,3))

    return table

# Register Table class reader
registry.register_reader('obssearch-tab', Table, read_obssearch_page)
