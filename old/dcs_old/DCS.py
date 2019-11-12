#! /usr/bin/env python
import mechanicalsoup
from getpass import getpass
import requests
from bs4 import BeautifulSoup
from urlpath import URL
from pathlib import Path
from astropy.config import get_cache_dir
import json
from collections import deque
import tempfile
from astropy.utils.data import download_file, clear_download_cache, is_url_in_cache
from astropy.config import set_temp_cache
from astropy.table import Table,vstack
import uuid
import MIS
import OBS
import AOR
import POS
import STATUS
import DETAILS
import AORSEARCH
import OBSSEARCH
import eel

DEBUG = False
DCSURL = 'https://dcs.arc.nasa.gov'
#DCSURL = 'https://devweb2.arc.nasa.gov'

class DCS(object):
    def __init__(self, dcsurl=DCSURL, cachedir=get_cache_dir(), refresh_cache=False):
        self.browser = mechanicalsoup.StatefulBrowser()
        self.dcsurl = URL(dcsurl)

        self.cachedir = Path(cachedir)/'DCS'
        self.cachedir.mkdir(exist_ok=True)

        self.refresh_cache = refresh_cache
        
        #self._load_cookie()

        self.loggedin = False

        if refresh_cache:
            r = self.browser.open(dcsurl)
            if self._verify_login(r):
                self.loggedin = True

    def _clear_cache(self):
        with set_temp_cache(self.cachedir):
            clear_download_cache(None)
        

    def login(self,attempts=0):
        '''Login to dcs.  Prompts for username/password'''
        print()
        username = input('DCS User: ')
        password = getpass()
        
        self.browser.open(self.dcsurl)
        try:
            form = self.browser.select_form('form[action="./interface/doLogin.jsp"]')
        except mechanicalsoup.utils.LinkNotFoundError:
            raise RuntimeError('%s appears to be down.  Try again later.'%DCSURL)
            
        form['user_name'] = username
        form['user_password'] = password
        
        r = self.browser.submit_selected()

        if self._verify_login(r, username):
            print('DCS login successful')
            self.loggedin = True
            self._save_cookie()
            return r
        else:
            if attempts < 2:
                print('Invalid Username or Password')
                return self.login(attempts+1)
            else:
                raise RuntimeError('Invalid Username or Password')

    def _save_cookie(self, key='JSESSIONID',cookiefile='browser.cookie'):
        '''Save session cookie to cache'''
        jar = self.browser.get_cookiejar()
        value = jar[key]
        with open(self.cachedir/cookiefile,'w') as f:
            json.dump({key:value},f)
        return jar

    def _load_cookie(self, key='JSESSIONID',cookiefile='browser.cookie'):
        '''Load session cookie from cache'''
        jar = self.browser.get_cookiejar()
        try:
            with open(self.cachedir/cookiefile,'r') as f:
                cookies = json.load(f)
            jar.set(key,cookies[key])
        except (FileNotFoundError,json.decoder.JSONDecodeError,KeyError):
            return None
        return jar

    def _verify_login(self, response, username='none'):
        '''Return true if successfully logged in, else false'''
        if ('Login Failed:' in response.text) or ('Must Login' in response.text):
            return False
        else:
            soup = BeautifulSoup(response.text,'html.parser')
            logdiv = soup.find('div', class_='SignIn')
            try:
                user = logdiv.table.tr.td.b.string
                if user and (username == 'none'):
                    # we are logged in as someone
                    return True
            except:
                return False
            if user == username:
                return True
            else:
                return False

    def _get_query_tmp(self,query, submit=None, maxlen=100):
        '''Return tmpfile name and file:// prefix version for caching'''
        tmpdir = tempfile.gettempdir()

        if submit == 'TAR':
            tmpfile = str(query)
            tmpfile = tmpfile.split('posdownload/')[1]
        elif submit:
            # submit button vs. direct query
            tmpfile = query.name
            if hasattr(query,'query') and query.query:
                submit = submit.copy()
                submit['query'] = query.query
            substr = ['%s_%s' % (k,v) for k,v in submit.items()]
            substr = ';'.join(substr)
            tmpfile = '%s:%s'%(tmpfile,substr)
        else:
            tmpfile = query.query


        # truncate tmpfile by hashing in case filename is too long
        if len(tmpfile) > maxlen:
            tmpfile = str(uuid.uuid3(uuid.NAMESPACE_URL,tmpfile))
        tmpfile = Path(tmpdir)/tmpfile
        streamfile = 'file://%s' % str(tmpfile)

        return tmpfile,streamfile

    def _get_from_cache(self,query, submit=None):
        '''If file exists in cache, return it'''
        _,streamfile = self._get_query_tmp(query,submit)
        with set_temp_cache(self.cachedir):
            if is_url_in_cache(streamfile):
                return download_file(streamfile,show_progress=DEBUG,cache=True)
            else:
                return False

    def _copy_to_cache(self, response, query, submit=None):
        '''Copy response contents to cache'''
        tmpfile, streamfile = self._get_query_tmp(query,submit)

        with tmpfile.open('wb') as f:
            f.write(response.content)
                        
        with set_temp_cache(self.cachedir):
            if self.refresh_cache:
                clear_download_cache(streamfile)
            cfile = download_file(streamfile,show_progress=DEBUG,cache=True)

        # remove tmpfile
        tmpfile.unlink()
        return cfile


    def _queryDCS(self, query, form, submit=None):
        '''Perform DCS search'''
        if not self.refresh_cache:
            # if we are not refreshing the cache, and the cache exists, return plan
            cfile = self._get_from_cache(query, submit)
            if cfile:
                return cfile

        if not self.loggedin:
            self.login()

        if form == 'DIRECT':
            # direct download
            r = self.browser.open(query)
        else:
            # simulate form submission
            self.browser.open(query)
            form = self.browser.select_form(form)

            if submit:
                if 'FILL' in submit and submit['FILL']:
                    for k,v in submit.items():
                        if k in ('name','value','FILL'):
                            continue
                        form[k] = v
                form[submit['name']] = submit['value']

                
            
            r = self.browser.submit_selected()

        cfile = self._copy_to_cache(r, query, submit)
        return cfile


    def getAOR(self,aorid,
               rel='observationPlanning/DbProxy/getObsPlanAORs.jsp',
               form='DIRECT',
               obsblock=None,
               raw=False,
               include_meta=True):
        '''Download AOR xml.  If obsblock specified, return table of obsblock only'''
        query = (self.dcsurl/rel).with_query({'origin':'GI','obsplanID':aorid})
        cfile = self._queryDCS(query,form)
        if raw:
            return cfile

        #tab = Table.read(cfile,format='aor-tab',obsblock=obsblock,d=self)
        tab = Table.read(cfile,format='aor-tab',obsblock=obsblock,d=self,include_meta=include_meta)
        return tab
        
    def getObsPlan(self, aorid,
                   rel='observationPlanning/observingPlanDetails.jsp',
                   form='form[action="DbProxy/getObsPlan.jsp"]',
                   raw=False):
        '''Download obsplan xml'''
        query = (self.dcsurl/rel).with_query({'obspID':aorid})

        cfile = self._queryDCS(query, form)
        if raw:
            return cfile

        tab = Table.read(cfile,format='obs-tab')
        return tab

    def getFlightSeries(self, series,
                        rel='flightPlan/flightPlanSearch.jsp',
                        form='form[action="flightPlanSearch.jsp"]',
                        submit={'name':'flightNum'},
                        local=None):
        '''Get names of flight in series'''
        if local:
            local = Path(local).glob('*.mis')
            flightids = [f.stem.replace('_INIT','') for f in local]
            return flightids
        
        submit['value'] = series
        query = (self.dcsurl/rel)

        cfile = self._queryDCS(query,form,submit)

        # read flight names from series
        with open(cfile,'r') as f:
            soup = BeautifulSoup(f.read(),'html.parser')
            idtags = soup.find_all(attrs={'name':'fpID','type':'checkbox'})
            ids = [id['value'] for id in idtags]
        return ids
        
    def getFlightPlan(self, flightid,
                      rel='flightPlan/downloadFlightPlanFile.jsp',
                      form='DIRECT',
                      raw=False,
                      local=None,
                      utctab=False):
        '''Download mis file'''
        if local:
            try:
                cfile = list(Path(local).glob('%s*.mis'%flightid))[0]
            except IndexError:
                cfile = Path(local)
            if raw:
                return local
            if utctab:
                _,tab = MIS.get_legs(cfile)
            else:
                tab = Table.read(cfile,format='mis-tab')
            return tab
            
        
        query = (self.dcsurl/rel).with_query({'fpid':flightid,'fileType':'mis'})

        cfile = self._queryDCS(query, form)
        
        with open(cfile,'r') as f:
            text = f.read()
            if 'No Matching MIS found' in text:
                print('No Matching MIS found for Flight Plan %s' % flightid)

                # remove from cache
                _,streamfile = self._get_query_tmp(query)
                with set_temp_cache(self.cachedir):
                    clear_download_cache(streamfile)
                return None

        if raw:
            return cfile

        if utctab:
            _,tab = MIS.get_legs(cfile)
        else:
            tab = Table.read(cfile,format='mis-tab')
        return tab

    def getPOS(self,flightid,
               rel='observationPlanning/SaveTargetPos.jsp',
               form='DIRECT',
               raw=False,
               guide=False):
        '''Download pos file'''

        if guide:
            query = (self.dcsurl/rel).with_query({'flightPlanID':flightid,
                                                  'posFormat':'MCCS',
                                                  'includeGS':'YES'})
        else:
            query = (self.dcsurl/rel).with_query({'flightPlanID':flightid,
                                                  'posFormat':'MCCS'})
        cfile = self._queryDCS(query, form)
        if raw:
            return cfile

        tab = Table.read(cfile,format='pos-tab',guide=guide)
        return tab

    def getPOSBundle(self,flightid,
                     rel='flightPlan/getTargetPosBundle.jsp',
                     form='DIRECT'):
        '''Download individual pos files as tar bundle'''
        query = (self.dcsurl/rel).with_query({'fpID':flightid,
                                              'getSelectedOnly':'TRUE',
                                              'includeGS':'YES'})

        cfile = self._queryDCS(query, form)
        # cfile here is a download page
        #  need to parse for location.href
        with open(cfile) as f:
            soup = BeautifulSoup(f.read(),'html.parser')
            if 'location.href' in soup.script.text:
                fname = soup.script.text.strip().split('location.href')[1]
                fname = fname[3:-2]
        cfile = self._queryDCS(self.dcsurl/fname, form, 'TAR')
        return cfile                                                  

    def getProposal(self,propid,
                    rel='proposalAccess/viewProposalDoc.jsp',
                    form='DIRECT'):
        '''Download pdf'''
        query = (self.dcsurl/rel).with_query({'propid':propid})
        cfile = self._queryDCS(query, form)
        return cfile

    def getObsStatus(self,propid,
                     rel='observationPlanning/observingPlanSearch.jsp',
                     form='form[name="formSearchObsPlan"]',
                     submit={'name':'getStatus','value':'Submit'}):
        '''Get status of observations'''
        query = (self.dcsurl/rel).with_query({'planID':propid})
        cfile = self._queryDCS(query,form,submit)
        #print(cfile)

        tab = Table.read(cfile, format='stat-tab')
        return tab

    def getObsDetails(self,propid,
                      rel='observationPlanning/observingPlanDetails.jsp',
                      form='DIRECT'):
        '''Download obsplan detail metadata'''
        query = (self.dcsurl/rel).with_query({'obspID':propid})
        cfile = self._queryDCS(query,form)

        tab = Table.read(cfile,format='obsdetails-tab')
        return tab


    def getAORDetails(self,aorid,
                      rel='observationPlanning/AORDetails.jsp',
                      form='DIRECT'):
        '''Download AORID detail metadata and guide stars'''
        query = (self.dcsurl/rel).with_query({'aorID':aorid})
        cfile = self._queryDCS(query,form)

        tab = Table.read(cfile,format='aordetails-tab')
        return tab
    
    def getObsBlockSearch(self,blockid,
                          rel='flightPlan/observationBlockSearch.jsp',
                          form='DIRECT'):
        '''Download ObsBlock Leg Summary'''
        query = (self.dcsurl/rel).with_query({'obsBlkIDs':blockid,
                                              'resPerPage':500,
                                              'getObsBlkInfo':'Submit'})
        cfile = self._queryDCS(query,form)

        tab = Table.read(cfile,format='obssearch-tab')
        return tab
        

    def getAORSearch(self,instrument='ALL',cycle='-1',aorids='',
                     rel='observationPlanning/AORSearch.jsp',
                     #form='form[action="AORSearch.jsp"]',submit={'name':'getAORs','value':'Submit'},
                     form='DIRECT'):
        '''Get all AORs by instrument for a flight cycle'''
        if instrument in ('HAWC+','HAWC'):
            instrument = 'HAWC_PLUS'

        formfill = {'cycleID':str(cycle),
                    'instrument':instrument,
                    #'planID':'',
                    #'piFirstName':'',
                    #'piLastName':'',
                    #'newTarget':'',
                    #'inputRA':'',
                    #'inputDec':'',
                    #'inputRadius':'',
                    'aorIDs':aorids,
                    'resPerPage':str(500),
                    'spectel1':'ALL',
                    'spectel2':'ALL',
                    'mode':'ALL',
                    'targetType':'ALL',
                    'inputEquinox':'2000',
                    'aorState':'ALL',
                    'getAORs':'Submit'}
        '''
        submit['FILL'] = True
        formfill = {'cycleID':str(cycle),
                    'instrument':instrument,
                    'planID':'',
                    'piFirstName':'',
                    'piLastName':'',
                    'newTarget':'',
                    'inputRA':'',
                    'inputDec':'',
                    'inputRadius':'',
                    'aorIDs':'',
                    'resPerPage':str(500),
                    'spectel1':'ALL',
                    'spectel2':'ALL',
                    'mode':'ALL',
                    'targetType':'ALL',
                    'inputEquinox':'2000',
                    'aorState':'ALL'}
        submit.update(**formfill)
        '''
        
        query = (self.dcsurl/rel).with_query(formfill)
        cfile = self._queryDCS(query,form)
        
        tab = Table.read(cfile,format='aorsearch-tab')

        # check if multiple pages
        with open(cfile) as f:
            soup = BeautifulSoup(f.read(),'html.parser')

        try:
            pform = soup.find('form',attrs={'name':'formAORSearchNextTop'})
            options = pform.find('select',attrs={'name':'start'}).find_all('option')
        except AttributeError:
            # single page
            return tab

        if len(options) == 1:
            # only one table of 500 rows, so return
            return tab

        tabs = deque()
        tabs.append(tab)
        for option in options[1:]:
            formfill['start'] = option['value']
            query = (self.dcsurl/rel).with_query(formfill)
            cfile = self._queryDCS(query,form)
            ctab = Table.read(cfile,format='aorsearch-tab')
            tabs.append(ctab)
        tab = vstack(list(tabs))
        return tab

    
if __name__ == '__main__':
    #d = DCS(refresh_cache=False)
    #print(d.getProposal('07_0165'))

    #d = DCS()
    #d.getObsStatus('06_0011')

    #d = DCS()
    #d.getObsDetails('06_0011')

    #d = DCS()
    #d.getAORSearch('HAWC',cycle=7)

    #d = DCS()
    #d.getObsBlockSearch('OB_06_0058_03')

    d = DCS(refresh_cache=False)
    d.getFlightPlan('201909_HA_FRODO')
    #d.getAORDetails('07_0012_1')
    #d.getPOSBundle('201906_FO_DEBORAH')
