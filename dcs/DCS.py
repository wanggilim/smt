#! /usr/bin/env python
import mechanicalsoup
from getpass import getpass
from bs4 import BeautifulSoup
from urlpath import URL
from pathlib import Path
import json
from http.cookies import SimpleCookie
from collections import deque
import tempfile
from astropy.utils.data import download_file, clear_download_cache, is_url_in_cache
from astropy.config import set_temp_cache,get_cache_dir
#from astropy.table import Table,vstack
import datetime
import time
import uuid
from configparser import ConfigParser
from peewee import SqliteDatabase
from . import DBmodels

'''
import MIS
import OBS
import AOR
import POS
import STATUS
import DETAILS
import AORSEARCH
import OBSSEARCH
'''

DEBUG = False
DCSURL = 'https://dcs.arc.nasa.gov'
#DCSURL = 'https://devweb2.arc.nasa.gov'

HTMLPARSER = 'lxml'


class DCS(object):
    def __init__(self, username=None, password=None,
                 dcsurl=DCSURL,
                 cachedir=get_cache_dir(),
                 refresh_cache=False,
                 modelcfg = 'DBmodels.cfg',
                 models = ('AOR','MIS','FAOR','AORSEARCH')):
        
        self.browser = mechanicalsoup.StatefulBrowser()
        self.dcsurl = URL(dcsurl)

        self.cachedir = Path(cachedir)/'DCS'
        self.cachedir.mkdir(exist_ok=True)

        self.refresh_cache = refresh_cache

        self.username = username
        self.password = password
        
        #self._load_cookie()

        # setup database with models
        if isinstance(modelcfg,str):
            # process modelcfg
            self.mcfg = ConfigParser()
            self.mcfg.read(modelcfg)
        else:
            # assume configparser object
            self.mcfg = modelcfg
        self.db = self._initialize_database(self.mcfg, models)
        
        self.loggedin = False

        if refresh_cache:
            r = self.browser.open(dcsurl)
            if self._verify_login(r):
                self.loggedin = True

    def _clear_cache(self):
        with set_temp_cache(self.cachedir):
            clear_download_cache(None)

    def _force_db_sync(self):
        urlmap = self.cachedir/'astropy/download/py3/urlmap.dir'
        if not urlmap.exists():
            return False

        with set_temp_cache(self.cachedir):
            with open(urlmap,'r') as f:
                for line in f:
                    streamfile = line.split("', (")[0][1:]
                    if 'fileType=misxml&fpid=' in streamfile:
                        # MIS file
                        cfile = download_file(streamfile,show_progress=False,cache=True)
                        rows = MIS.to_rows(cfile, self.mcfg['MIS'])
                        if rows:
                            MIS.replace_rows(self.db, rows)
                    elif 'obsplanID=' in streamfile:
                        # AOR file
                        cfile = download_file(streamfile,show_progress=False,cache=True)
                        rows = AOR.to_rows(cfile, self.mcfg['AOR'])
                        if rows:
                            AOR.replace_rows(self.db, rows)
                    else:
                        continue
                    

    def _initialize_database(self, mcfg, models=('AOR','MIS','FAOR','AORSEARCH')):
        """Initialize database, create models dynamically, and register models globally"""
        db_file = mcfg['DEFAULT']['db_file']
        db = SqliteDatabase(db_file)

        # generate models
        mods = [DBmodels.ModelFactory(name, mcfg, db) for name in models]

        # create database
        with db:
            db.create_tables(mods,safe=True)

        self.ACTIVE_MODELS = {name:model for name,model in zip(models,mods)}
            
        # register model names globally
        for name,model in self.ACTIVE_MODELS.items():
            globals()[name] = model

        return db

    def get_active_models(self):
        if hasattr(self,'ACTIVE_MODELS'):
            return self.ACTIVE_MODELS
        else:
            return None
        

    def login(self,username=None,password=None,attempts=0):
        '''Login to dcs.  Prompts for username/password'''
        if username is None:
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
            self.username = username
            self.password = password
            #self._save_cookie()
            return r
        else:
            if attempts < 2:
                print('Invalid Username or Password')
                return self.login(username=None,attempts=attempts+1)
            else:
                raise RuntimeError('Invalid Username or Password')

    def _save_cookie(self, key='JSESSIONID',cookiefile='browser.cookie',
                     expkey='expires',expiry=20):
        '''Save session cookie to cache'''
        jar = self.browser.get_cookiejar()
        
        if key not in jar:
            return jar

        cookie = SimpleCookie()
        
        for c in jar:
            if c.name != key:
                continue
            cookie = c
        
        #value = jar[key]

        expires = datetime.datetime.utcnow() + datetime.timedelta(minutes=expiry)
        cookie.expires = expires.strftime("%a, %d %b %Y %H:%M:%S GMT")
        
        #cookie = [{key:value,expkey:expires.strftime("%a, %d %b %Y %H:%M:%S GMT")}]
        
        with open(self.cachedir/cookiefile,'w') as f:
            #json.dump(cookie,f)
            json.dump(cookie.__dict__,f)
        return jar

    def _load_cookie(self, key='JSESSIONID',cookiefile='browser.cookie',expkey='expires'):
        '''Load session cookie from cache'''
        jar = self.browser.get_cookiejar()
        try:
            with open(self.cachedir/cookiefile,'r') as f:
                cookie = json.load(f)
        except (FileNotFoundError,json.decoder.JSONDecodeError,KeyError):
            return None

        #print(cookie)
        '''
        print(c)
        cookie = SimpleCookie()
        cookie.name = key
        cookie[key] = c['value']
        cookie['value'] = c['value']
        for k in (expkey,'path','domain'):
            setattr(cookie,k,c[k])
            cookie[key][k] = c[k]

        print(cookie)
        exit()
        '''
        #jar.set(key,cookies[key])

        expires = datetime.datetime.strptime(cookie[expkey],"%a, %d %b %Y %H:%M:%S GMT")
        delta = expires - datetime.datetime.utcnow()

        if delta.total_seconds() > 0:
            # cookie not expired
            #  but must erase expires from cookie
            value = cookie['value']
            jar.set(key,value)

        return jar

    def _verify_login(self, response, username='none'):
        '''Return true if successfully logged in, else false'''
        if ('Login Failed:' in response.text) or ('Must Login' in response.text):
            return False
        else:
            soup = BeautifulSoup(response.text,HTMLPARSER)
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

    def _get_from_cache(self, query, submit=None):
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
            # if we are not refreshing the cache, and the cache exists, return hit
            cfile = self._get_from_cache(query, submit)
            if cfile:
                return cfile

        if not self.loggedin:
            self.login(username=self.username,password=self.password)

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



    def getAOR(self, aorID):
        """Get AOR row."""
        split = aorID.split('_')
        if len(split) == 2:
            # this is a planID
            return self.getAORs_by_planID(aorID)
        planID = '_'.join(split[:-1])
        
        if self.refresh_cache:
            # bypass DB
            return self.getObsPlan(planID, aorID=aorID)
        
        aor = DCS._query_AOR_table(aorID=aorID)
        if not aor:
            return self.getObsPlan(planID,aorID=aorID)
        return json.dumps(aor)

    def getAORs_by_planID(self, planID):

        if self.refresh_cache:
            # bypass DB
            return self.getObsPlan(planID)
        
        aors = DCS._query_AOR_table(planID=planID)
        if not aors:
            return self.getObsPlan(planID)

        return json.dumps(aors)

    
    def getAORs_by_ObsBlk(self, ObsBlk):
        planID = '_'.join(ObsBlk.split('_')[1:-1])
        
        if self.refresh_cache:
            # bypass DB
            return self.getObsPlan(planID,ObsBlk=ObsBlk)

        aors = DCS._query_AOR_table(ObsBlk=ObsBlk)
        if not aors:
            return self.getObsPlan(planID,ObsBlk=ObsBlk)

        return json.dumps(aors)

    @staticmethod
    def _query_AOR_table(aorID=None,ObsBlk=None,planID=None):
        """Perform DB lookup on AOR table"""
        if aorID:
            try:
                aor = AOR.get(AOR.aorID==aorID)
                aors = aor.__data__
            except:
                aors = None
        elif planID:
            aors = AOR.select().where(AOR.planID==planID)
            aors = [aor.__data__ for aor in aors]
            aors = aors if aors else None
        elif ObsBlk:
            aors = AOR.select().where(AOR.ObsBlk==ObsBlk).order_by(AOR.order,AOR.aorID)
            aors = [aor.__data__ for aor in aors]
            aors = aors if aors else None
        else:
            aors = None
        return aors
    
    
    def getObsPlan(self, planID,
                   rel='observationPlanning/DbProxy/getObsPlanAORs.jsp',
                   form='DIRECT',
                   raw=False,
                   aorID=None,
                   ObsBlk=None):
        '''Download AOR xml.'''
        query = (self.dcsurl/rel).with_query({'origin':'GI','obsplanID':planID})
        cfile = self._queryDCS(query,form)

        aorcfg = self.mcfg['AOR']
        rows = AOR.to_rows(cfile, aorcfg)
        if not rows:
            return None
        AOR.replace_rows(self.db, rows)
        
        if raw:
            return cfile

        aors = DCS._query_AOR_table(aorID=aorID,ObsBlk=ObsBlk,planID=planID)

        return json.dumps(aors)
        
    def getObsPlanXML(self, planID,
                      rel='observationPlanning/observingPlanDetails.jsp',
                      form='form[action="DbProxy/getObsPlan.jsp"]',
                      raw=False):
        '''Download obsplan xml'''
        query = (self.dcsurl/rel).with_query({'obspID':planID})
        cfile = self._queryDCS(query, form)

        if raw:
            return cfile
        
        return cfile

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
            soup = BeautifulSoup(f.read(),HTMLPARSER)
            idtags = soup.find_all(attrs={'name':'fpID','type':'checkbox'})
            ids = [id['value'] for id in idtags]
        return ids
        
    def getFlightPlan(self, flightid,
                      rel='flightPlan/downloadFlightPlanFile.jsp',
                      form='DIRECT',
                      local=None,
                      raw=False):
        '''Download mis file'''
        if local:
            try:
                cfile = list(Path(local).glob('%s*.mis'%flightid))[0]
            except IndexError:
                cfile = Path(local)

        else:
            query = (self.dcsurl/rel).with_query({'fpid':flightid,'fileType':'misxml'})
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

        miscfg = self.mcfg['MIS']
        rows = MIS.to_rows(cfile, miscfg)
        if not rows:
            return None
        MIS.replace_rows(self.db, rows)

        if raw:
            return cfile

        legs = MIS.select().where(MIS.FlightPlan==flightid).order_by(MIS.Leg)
        legs = (leg.__data__ for leg in legs)
        return json.dumps(list(legs))
        

    def getPOS(self,flightid,
               rel='observationPlanning/SaveTargetPos.jsp',
               form='DIRECT',
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

        return cfile
    
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
            soup = BeautifulSoup(f.read(),HTMLPARSER)
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

        #tab = Table.read(cfile, format='stat-tab')
        #return tab
        return cfile

    def getObsDetails(self,propid,
                      rel='observationPlanning/observingPlanDetails.jsp',
                      form='DIRECT'):
        '''Download obsplan detail metadata'''
        query = (self.dcsurl/rel).with_query({'obspID':propid})
        cfile = self._queryDCS(query,form)

        #tab = Table.read(cfile,format='obsdetails-tab')
        #return tab
        return cfile


    def getAORDetails(self,aorid,
                      rel='observationPlanning/AORDetails.jsp',
                      form='DIRECT'):
        '''Download AORID detail metadata and guide stars'''
        query = (self.dcsurl/rel).with_query({'aorID':aorid})
        cfile = self._queryDCS(query,form)

        #tab = Table.read(cfile,format='aordetails-tab')
        #return tab
        return cfile
    
    def getObsBlockSearch(self,blockid,
                          rel='flightPlan/observationBlockSearch.jsp',
                          form='DIRECT'):
        '''Download ObsBlock Leg Summary'''
        query = (self.dcsurl/rel).with_query({'obsBlkIDs':blockid,
                                              'resPerPage':500,
                                              'getObsBlkInfo':'Submit'})
        cfile = self._queryDCS(query,form)

        #tab = Table.read(cfile,format='obssearch-tab')
        #return tab
        return cfile
        

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

        #tab = Table.read(cfile,format='aorsearch-tab')

        # check if multiple pages
        with open(cfile) as f:
            soup = BeautifulSoup(f.read(),HTMLPARSER)

        try:
            pform = soup.find('form',attrs={'name':'formAORSearchNextTop'})
            options = pform.find('select',attrs={'name':'start'}).find_all('option')
        except AttributeError:
            # single page
            return tab

        if len(options) == 1:
            # only one table of 500 rows, so return
            return tab

        #tabs = deque()
        #tabs.append(tab)
        cfiles = deque()
        cfiles.append(cfile)
        for option in options[1:]:
            formfill['start'] = option['value']
            query = (self.dcsurl/rel).with_query(formfill)
            cfile = self._queryDCS(query,form)
            #ctab = Table.read(cfile,format='aorsearch-tab')
            #print(cfile)
            #tabs.append(ctab)
            cfiles.append(cfiles)
            
        #tab = vstack(list(tabs))
        #return tab

        if self.db:
            aorsearchcfg = self.mcfg['AORSEARCH']
            rows = AORSEARCH.to_rows(cfile, aorsearchcfg)
            AORSEARCH.replace_rows(self.db, rows)
        
        return cfiles

    
if __name__ == '__main__':
    #d = DCS(refresh_cache=False)
    #print(d.getProposal('07_0165'))

    #d = DCS()
    #d.getObsStatus('06_0011')

    #d = DCS()
    #d.getObsDetails('06_0011')

    d = DCS(database=True,refresh_cache=False)
    d.getAOR('07_0149')

    exit()
    d = DCS()
    t = d.getAORSearch('HAWC',cycle=-1)
    print(t)
    print(t.colnames)

    #d = DCS()
    #d.getObsBlockSearch('OB_06_0058_03')

    #d = DCS(refresh_cache=True)
    #d.getAORDetails('07_0012_1')
    #d.getPOSBundle('201906_FO_DEBORAH')
    #
    #d = DCS(refresh_cache=True)
    #print(d.getAOR('90_0071',raw=True))
