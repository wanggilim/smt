#! /usr/bin/env python
from bs4 import BeautifulSoup,NavigableString
from astropy.table import Table
from pandas import read_html
from astropy.io import registry
import numpy as np
import astropy.units as u
from collections import deque

AORTABUNITDICT = {'Radius(arcmin)':u.arcmin,'ProperMotionRA(masyr)':u.mas/u.yr,'ProperMotionDec(masyr)':u.mas/u.yr,
                  'Magnitude (Type)':u.mag,'Magnitude2 (Type)':u.mag,'BminusV':u.mag}

def get_bs_attrs(bs):
    '''Get attrs from li bs'''
    attrs = {b.text.strip():b.next_sibling for b in bs if b}

    for k,v in attrs.copy().items():
        # get rid of :
        if k[-1] == ':':
            del attrs[k]
            k = k[:-1].strip()
            attrs[k] = v
        try:
            value = v.string.strip()
            attrs[k] = value
        except AttributeError:
            if v is None:
                continue
            if v.name == 'ul':
                # value is list
                value = v.find_all('li')
                value = [v.text.split('\n')[0].strip() for v in value]
                attrs[k] = value
        # get rid of '['
        if value[-1] == '[':
            attrs[k] = value[:-1].strip()

    return attrs

def get_li_attrs(lis):
    '''Get attrs from li'''
    names = deque()
    for li in lis:
        name = li.text.split('\n')[0]
        name = name.strip()
        #print(li,name)
        if name[-1] == '[':
            name = name[:-1]
        names.append(name)
    return list(names)

def get_obstable(htable):
    '''Get observation table'''
    table = read_html(str(htable))
    table = Table.from_pandas(table[0])
    table.pprint()
    if 'Requested Time(minutes)' in table.colnames:
        table.rename_column('Requested Time(minutes)','Requested Time')
        table['Requested Time'].unit = u.minute
    if 'Water Vapor(microns)' in table.colnames:
        table.rename_column('Water Vapor(microns)','Water Vapor')
        table['Water Vapor'].unit = u.micron
        
    return table

def read_obsdetails(filename):
    with open(filename,'r') as f:
        soup = BeautifulSoup(f.read(),'html.parser')

    # h3 has metadata title
    headers = soup.find_all('h3')
    metadata = {}
    for header in headers:
        htext = header.text.strip()
        if htext in ('Summary','Proposal Information'):
            bs = header.next_sibling.find_all('b')
            attrs = get_bs_attrs(bs)
        elif htext in ('Investigators','Support Scientists'):
            lis = header.next_sibling.find_all('li')
            attrs = get_li_attrs(lis)
        elif htext in ('Proposed Observations',):
            continue
            htable = header.next_sibling.find('table')
            table = get_obstable(htable)

        elif 'Note:' in htext:
            attrs = header.next_sibling.text
            htext = 'Note'

        else:
            attrs = {}

        metadata[htext] = attrs

    for header in headers:
        htext = header.text.strip()
        if htext != 'Proposed Observations':
            continue
        htable = header.next_sibling.find('table')
        table = get_obstable(htable)

    table.meta = metadata
    return table

def read_aordetails(filename):
    with open(filename,'r') as f:
        soup = BeautifulSoup(f.read(),'html.parser')

    # h3 has metadata title
    headers = soup.find_all('h3')
    metadata = {}
    for header in headers:
        htext = header.text.strip()
        if htext in ('Scheduling Concerns','Instrument Information','Target Information'):
            lis = header.findNext('ul').findChildren('li')
            attrs = {x.text.split(':')[0].strip():x.text.split(':')[1].strip() for x in lis}
            if 'Target Position' in attrs:
                for li in lis:
                    ltext = li.text.split(':')
                    l = ltext[0]
                    b = ltext[1:]
                    if l.strip() == 'RA':
                        attrs['RA'] = ':'.join(b)
                    if l.strip() == 'Dec':
                        attrs['Dec'] = ':'.join(b)

                attrs['Target Position'] = {'RA':attrs['RA'].strip(),'Dec':attrs['Dec'].strip()}
            if 'Dec' in attrs:
                del attrs['Dec']
            if 'RA' in attrs:
                del attrs['RA']
        else:
            attrs = {}

        metadata[htext] = attrs

    # get pos table
    ths = soup.find_all('th')
    for th in ths:
        if th.text.strip() == 'Imager':
            # this is unique
            htable = th.parent.parent
            break

    table = read_html(str(htable))
    table = Table.from_pandas(table[0])
    table.meta = metadata

    for col in table.colnames:
        if col in AORTABUNITDICT:
            table[col].unit = AORTABUNITDICT[col]
            table.rename_column(col,col.split('(')[0].strip())
    
    return table
        
            
    '''
    # ul has metadata
    uls = soup.find_all('ul')
    exit()
    for header in headers:
        print(header.find('ul'))

    exit()

    # table just has tac grade and propid
    #for li in soup.find('h3',string='Proposal Information').next_element().find_all('li'):
    #print(['Grade' in str(li) for li in soup.find_all('li')])
    try:
        b = soup.find('b',string='TAC Grade:  ')
        grade = float(b.next_sibling)
    except AttributeError:
        grade = np.nan

    try:
        b = soup.find('b',string='Priority: ')
        priority = b.next_sibling.strip()
    except:
        priority = ''
        
    b = soup.find('b',string='Proposal ID:  ')
    propid = b.next_sibling.string

    table = Table([{'PlanID':propid,'TAC Grade':grade,'Priority':priority}])
    
    return table
    '''

# Register Table class reader
registry.register_reader('obsdetails-tab', Table, read_obsdetails)
registry.register_reader('aordetails-tab', Table, read_aordetails)
