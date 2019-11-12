#! /usr/bin/env python
from bs4 import BeautifulSoup,NavigableString
from astropy.table import Table,Column
from pandas import read_html
from astropy.io import registry
import numpy as np
import astropy.units as u
from collections import deque


def read_aorsearch_page(filename):
    '''Read AOR search page results'''
    with open(filename,'r') as f:
        soup = BeautifulSoup(f.read(),'html.parser')

    ths = soup.find_all('th')
    # this table should be unique
    try:
        htable = list(filter(lambda x: 'NAIF_ID' in x.text,ths))[0]
    except IndexError:
        # no aorid found
        return Table()

    htable = htable.parent.parent
    table = read_html(str(htable))

    try:
        table = Table.from_pandas(table[0])
    except KeyError:
        table = Table.from_pandas(table)

    # get rid of 'ViewPlan'
    aorids = Column([x.split('View')[0].strip() for x in table['AORID']],name='AORID')
    table.replace_column('AORID',aorids)
    aors = Column(['_'.join(x.split('_')[:-1]) for x in aorids],name='AOR')
    table.add_column(aors,index=1)

    return table

# Register Table class reader
registry.register_reader('aorsearch-tab', Table, read_aorsearch_page)

