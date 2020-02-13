#! /usr/bin/env python
import configparser
import argparse
import eel
from dcs import DCS
from dossier import TEX
from dossier.dossier import parse_flightid
from astropy.table import Table
from astropy.table.jsviewer import DEFAULT_CSS
from io import StringIO
from bs4 import BeautifulSoup
from pathlib import Path
from jinja2 import Environment,FileSystemLoader
from functools import partial
from astropy.utils.console import ProgressBar
import cloudpickle
import pickle

env = None
#d = None
d = DCS.DCS()

def make_html_table(table,
                    css=DEFAULT_CSS, jskwargs={'use_local_files': True},
                    tableid=None,
                    table_class="table table-striped table-bordered table-hover",
                    addl_table_class=None,
                    show_row_index=None):
    if show_row_index:
        table = table._make_index_row_display_table(show_row_index)

    if addl_table_class:
        table_class = ' '.join((table_class,addl_table_class))

    buff = StringIO()
    with buff as f:
        table.write(f, format='jsviewer', css=css, jskwargs=jskwargs,
                    table_id=tableid, table_class=table_class)
        html_str = buff.getvalue()

    soup = BeautifulSoup(html_str,'lxml')
    return soup

def wrap_card(soup, content='table', title=None, style=None):
    card_body = soup.new_tag('div', **{'class':'card-body'})
    table = getattr(soup.body,content)
    if style:
        table['style'] = style
    table.wrap(card_body)

    card = soup.new_tag('div', **{'class':'card'})
    soup.body.div.wrap(card)

    if title:
        h5 = soup.new_tag('h5', **{'class':'card-header h5'})
        h5.string = title
        soup.body.div.div.insert_before(h5)

    card = soup.body.div
    return card


def proc_MIS(tup,d):
    flightid,name = tup
    mis = d.getFlightPlan(flightid)

    mis = [{k:v for k,v in row.items() if k in ('Leg','ObsBlkID','Target','Start')} for row in mis]

    mis = Table(rows=mis)
    tableid = 'table-%s' % (flightid)
    soup = make_html_table(mis,tableid=tableid,addl_table_class='table-MIS')

    flight = {'flightid':flightid,'name':name,'table':str(soup.body.table)}
    return flight


@eel.expose
def prompt_login():
    return 'login'

@eel.expose
def DCS_search(formID, flightid, template_file='panel.html'):

    #global d
    #d = DCS.DCS()
    
    flightids, names, series = parse_flightid(flightid, d, local=False)

    proc_func = partial(proc_MIS,d=d)
    flights = ProgressBar.map(proc_func,list(zip(flightids,names)),multiprocess=False)

    #card = wrap_card(soup, title=flightid) #style="width:100%"
    template = env.get_template(template_file)
    page = template.render(series=series,flights=flights)
    tableids = ['#table-%s'%flight['flightid'] for flight in flights]
    return '#flight-panel', tableids, page


def proc_AOR(tup,flightid,d):
    obsblk,leg = tup

    if obsblk in (None,'--','','None'):
        return None
    planid = '_'.join(obsblk.split('_')[1:3])
    aor = d.getAORs(obsblk)

    details = TEX.make_details(aor, tex=False)

    tableid = 'table-%s-leg%i-%s' % (obsblk,leg,flightid)
    soup = make_html_table(details,tableid=tableid,addl_table_class='table-AOR')

    # shorten a few things
    tab = str(soup.body.table)
    tab = tab.replace('Lissajous','LIS')
    tab = tab.replace('Polarimetry','POL')
    
    detail = {'obsblk':obsblk, 'planid':planid, 'table':tab, 'leg':leg,
              'tableid':tableid, 'tabid':'tab-%s'%tableid}
    return detail

    
@eel.expose
def getAORs(obsblks, legs, flightid, template_file='aor.html'):
    proc_func = partial(proc_AOR,d=d,flightid=flightid)
    legs = [int(x) for x in legs]
    details = ProgressBar.map(proc_func,list(zip(obsblks,legs)),multiprocess=False)

    template = env.get_template(template_file)
    page = template.render(details=details)

    tableids = ['#%s'%detail['tableid'] for detail in details]
    return '#aor-panel', tableids, page


def main():
    parser = argparse.ArgumentParser(description='Run server and startup SOFIA Mission Toolbox')
    parser.add_argument('-cfg',type=str,default='server.cfg',help='Specify config file (default=server.cfg)')
    args = parser.parse_args()

    # Read in options
    config = configparser.ConfigParser()
    config.read(args.cfg)
    eelcfg = config['EEL']  # eel config section

    eel.init(eelcfg['webdir'])

    if 'jinja' in eelcfg:
        global env
        env = Environment(
            loader=FileSystemLoader(str(Path(eelcfg['webdir'],eelcfg['jinja']))))
    
    eel.start(eelcfg['start'],
              host=eelcfg['host'],
              port=int(eelcfg['port']),
              mode=eelcfg['mode'],
              jinja=eelcfg['jinja'])
    
    
if __name__ == '__main__':
    #DCS_search('#DCS-search','201909_HA')
    main()
