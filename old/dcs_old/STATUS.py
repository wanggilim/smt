#! /usr/bin/env python
from bs4 import BeautifulSoup
from astropy.table import Table,Column
from pandas import read_html
from astropy.io import registry
import numpy as np


def read_status_tab(filename):
    with open(filename,'r') as f:
        soup = BeautifulSoup(f.read(),'html.parser')

    # table is one with 'Completion' as header
    head = soup.find('th',string='Completion')
    htmltable = str(head.parent.parent)

    table = read_html(htmltable)
    table = Table.from_pandas(table[0])

    # remove links from first col
    table['PlanID'] = [x.split()[0] for x in table['PlanID']]

    # issue with CoI
    if not isinstance(table['Co-Is'][0],(str,np.str)):
        table.replace_column('Co-Is',Column([''],name='Co-Is'))
    return table

# Register Table class reader
registry.register_reader('stat-tab', Table, read_status_tab)
