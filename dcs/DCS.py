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
from astropy.utils.data import download_file, clear_download_cache, is_url_in_cache, get_cached_urls
from astropy.config import set_temp_cache,get_cache_dir
from astropy.table import Table,Column,vstack,join
import datetime
import time
import uuid
from configparser import ConfigParser
from peewee import SqliteDatabase
from pandas import DataFrame
from astropy.utils.console import ProgressBar
try:
    from . import MIS as dcsMIS
except ImportError:
    import MIS as dcsMIS
try:
    from . import DBmodels
except ImportError:
    import DBmodels
try:
    from . import OBSSEARCH
except ImportError:
    import OBSSEARCH

DEBUG = False
DCSURL = 'https://dcs.arc.nasa.gov'
#DCSURL = 'https://devweb2.arc.nasa.gov'

HTMLPARSER = 'lxml'


def _aorID_to_planID(aorID):
    if 'OB' in aorID:
        return _ObsBlkID_to_planID(aorID)
    split = aorID.split('_')
    if len(split) == 2:
        planID = aorID
    else:
        planID = '_'.join(split[:-1])
    return planID

def _ObsBlkID_to_planID(ObsBlkID):
    return '_'.join(ObsBlkID.split('_')[1:-1])

def _is_aorID(idstr):
    if 'OB' in idstr:
        return False
    split = idstr.split('_')
    if len(split) == 2:
        return False
    elif len(split) == 3:
        return True
    else:
        return False

def _is_planID(idstr):
    if 'OB' in idstr:
        return False
    split = idstr.split('_')
    if len(split) == 2:
        return True
    return False

def _is_ObsBlkID(idstr):
    if 'OB' in idstr and len(idstr.split('_')) == 4:
        return True
    return False


class DCS(object):
    def __init__(self, username=None, password=None,
                 dcsurl=DCSURL,
                 cachedir=get_cache_dir(),
                 refresh_cache=False,
                 #modelcfg = 'dcs/DBmodels.cfg',
                 modelcfg = str(Path(__file__).parent.resolve()/'DBmodels.cfg'),
                 models = ('AOR','MIS','FAOR','SCT','POS','GUIDE','AORSEARCH')):
        
        self.browser = mechanicalsoup.StatefulBrowser()
        self.dcsurl = URL(dcsurl)

        self.cachedir = Path(cachedir)/'DCS'
        self.cachedir.mkdir(exist_ok=True)

        self.refresh_cache = refresh_cache

        self.username = username
        self.password = password
        
        #self._load_cookie()

        # get model params
        if isinstance(modelcfg,str):
            # process modelcfg
            self.mcfg = ConfigParser()
            self.mcfg.read(modelcfg)
        else:
            # assume configparser object
            self.mcfg = modelcfg
            
        # setup database with model config
        self.models = models
        self.db = self._initialize_database(models)
        
        self.loggedin = False

        if refresh_cache:
            r = self.browser.open(dcsurl)
            if self._verify_login(r):
                self.loggedin = True

    def _clear_cache(self):
        """Remove downloaded DCS files from cachedir (default=$HOME/.astropy/cache/DCS)"""
        with set_temp_cache(self.cachedir):
            clear_download_cache(None)

    def _get_db_file(self):
        return Path(self.cachedir)/self.mcfg['DEFAULT']['db_file']
            
    def _remove_db(self):
        db_file = self._get_db_file()
        print('Removing %s' % str(db_file))
        #db_file.unlink(missing_ok=True)
        db_file.unlink()
        return db_file
        

    def _force_db_sync(self):
        """Ingest all downloaded files in cache into database"""
        '''
        urlmap = self.cachedir/'astropy/download/py3/urlmap.dir'
        if not urlmap.exists():
            urlmap = self.cachedir/'astropy/download/py3/urlmap.db'
            if not urlmap.exists():
                return False
        '''

        # remove db_file and re-init
        self._remove_db()
        self.db = self._initialize_database(self.models)

        print('Resyncing local DB with cache')
        with set_temp_cache(self.cachedir):
            cached_urls = get_cached_urls()
            for line in ProgressBar(cached_urls):
                #streamfile = line.split("', (")[0][1:]
                streamfile = line
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
                    grows = GUIDE.to_rows(cfile, self.mcfg['GUIDE'])
                    if grows:
                        GUIDE.replace_rows(self.db, grows)
                elif 'pos_' in streamfile:
                    cfile = download_file(streamfile,show_progress=False,cache=True)
                    rows = POS.to_rows(cfile, self.mcfg['POS'])
                    if rows:
                        POS.replace_rows(self.db, rows)
                    #rows = GUIDE.to_rows(cfile, self.mcfg['GUIDE'])
                    #if rows:
                    #    GUIDE.replace_rows(self.db, rows)
                else:
                    continue
                    

    def _initialize_database(self, models=('AOR','MIS','FAOR','SCT','POS','GUIDE','AORSEARCH')):
        """Initialize database, create models dynamically, and register models globally"""
        db_file = self._get_db_file()
        db = SqliteDatabase(str(db_file))

        # generate models
        mods = [DBmodels.ModelFactory(name, self.mcfg, db) for name in models]

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
        

    def login(self,username=None,password=None,attempts=0,gui=False):
        """Login to dcs.  Prompts for username/password"""
        if gui:
            # wait for gui response
            username,password = gui()
        
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


    def _add_missing_obsblks(self,rows):
        """If ObsBlkID is missing from AOR rows
        (usually, very old programs), try updating from DCS"""
        aor_plan_blk = [(row['aorID'],row['planID'],row['ObsBlkID']) \
                        for row in rows]
        if any([not x[-1] for x in aor_plan_blk]):
            # missing an obsblk
            plans = {x[1] for x in aor_plan_blk}
            # query for obsblk table
            obstabs = [self.getObsBlockSearch(planid=p) for p in plans]
            obstabs = vstack(obstabs)
            # make lookup
            blkdict = dict(zip(obstabs['aorID'],obstabs['ObsBlkID']))
            # add blkid to rows
            for row in rows:
                row['ObsBlkID'] = blkdict.get(row['aorID'],row['ObsBlkID'])
        return rows
            
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
            if 'posFormat' in tmpfile:
                prefix = 'pos_'
            else:
                prefix = ''
            tmpfile = '%s%s'%(prefix,str(uuid.uuid3(uuid.NAMESPACE_URL,tmpfile)))
        tmpfile = Path(tmpdir)/tmpfile
        streamfile = 'file://%s' % str(tmpfile)

        return tmpfile,streamfile

    def _get_from_cache(self, query, submit=None):
        '''If file exists in cache, return it'''
        _,streamfile = self._get_query_tmp(query,submit)
        with set_temp_cache(self.cachedir):
            if is_url_in_cache(streamfile):
                cfile = download_file(streamfile,show_progress=DEBUG,cache=True)
                # just in case, check if cfile exists
                if Path(cfile).exists():
                    return cfile
                else:
                    clear_download_cache(streamfile)
                    return False
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


    @staticmethod
    def as_table(rows):
        return AOR.as_table(rows)

    @staticmethod
    def as_pandas(rows):
        return AOR.as_pandas(rows)

    @staticmethod
    def as_json(rows):
        return AOR.as_json(rows)
    
    def getAORs(self, search, *args, **kwargs):
        """Get AOR row.
        If aorID is a planID (e.g. 07_0225) or ObsBlk (e.g. OB_07_0225_04),
        all aorIDs in that obs plan will be returned."""
        '''
        if isinstance(search,str) or len(search) == 1:
            # single aor/plan/obsblk
            if isinstance(search,(Table,Column,list,tuple)):
                search = str(search[0])
            if self.refresh_cache:
                # bypass DB
                return self._getObsPlan(search,*args,**kwargs)
            
            aors = self._query_AOR_table(search, *args, **kwargs)
            if not aors:
                # search not found in DB
                return self._getObsPlan(search,*args,**kwargs)
            
        else:
            # many aors/plans
            if self.refresh_cache:
                # bypass DB, but don't perform query yet
                cfiles = [self._getObsPlan(s,raw=True) for s in search]
            aors = self._query_AOR_table(search,*args,**kwargs)
            if not aors:
                cfiles = [self._getObsPlan(search,raw=True) for s in search]
                aors = self._query_AOR_table(search,*args,**kwargs)
        return aors
        '''
        return self._get(search, 'AOR', *args, **kwargs)

    def getGuideStars(self, search, *args, **kwargs):
        """Get GuideStars for aorID"""
        # search MUST be an aorID, planID, or ObsBlkID
        return self._get(search, 'GUIDE', *args, **kwargs)
        
        
    def getFAORs(self, search, *args, **kwargs):
        if kwargs.get('match'):
            # return faors in format to match aorids in dossier
            key_map = {'NODDWELL':'Nod','REPEATS':'Repeat','AORID':'aorID',
                       'DITHER':'Dithers','DITHSCALE':'Scale','FILENAME':'FAORfile',
                       'INTTIME':'IntTime','FDUR':'FDUR','SLIT':'Slit','TREQ':'TREQ',
                       'DURPERREW':'TREW','TLOS':'TLOS','TLSPN':'TLSPN',
                       'rewind':'Rewind','loop':'Loop','TARGET':'Name'}
            faors = self._get(search, 'FAOR', *args, **kwargs)
            faors = [{newk:f[k] for k,newk in key_map.items()} for f in faors]
            faors = {f['aorID']:f for f in faors}
            return faors
        else:
            return self._get(search, 'FAOR', *args, **kwargs)

    def getSCTs(self, search, *args, **kwargs):
        if kwargs.get('match'):
            # return scts in format to match aorids in dossier
            #key_map = {'NODDWELL':'Nod','REPEATS':'Repeat','AORID':'aorID',
            #           'DITHER':'Dithers','DITHSCALE':'Scale','FILENAME':'FAORfile',
            #           'INTTIME':'IntTime','FDUR':'FDUR','SLIT':'Slit','TREQ':'TREQ',
            #           'DURPERREW':'TREW','TLOS':'TLOS','TLSPN':'TLSPN',
            #           'rewind':'Rewind','loop':'Loop','TARGET':'Name'}
            key_map = {'FILENAME':'SCTfile'}
            scts = self._get(search, 'SCT', *args, **kwargs)
            scts = [{newk:f[k] for k,newk in key_map.items()} for f in scts]
            scts = {s['AORID']:s for s in scts}
            return scts
        else:
            return self._get(search, 'SCT', *args, **kwargs)


    def _get(self, search, clsname, *args, **kwargs):
        """Unified interface to DB and DCS fallback"""
        DB_funcs  = {'AOR':self._query_AOR_table,
                     'MIS':self._query_MIS_table,
                     'GUIDE':self._query_GUIDE_table,
                     'POS':self._query_POS_table,
                     'FAOR':self._query_FAOR_table,
                     'SCT':self._query_SCT_table}
        DCS_funcs = {'AOR':self._getObsPlan,
                     'MIS':self._getFlightPlan,
                     'GUIDE':self._getObsPlan,
                     'POS':self._getPOS,
                     'FAOR':lambda x: None,
                     'SCT':lambda x: None}

        getDB = DB_funcs[clsname]
        getDCS = DCS_funcs[clsname]

        if isinstance(search,str) or len(search) == 1:
            # single aor/plan/obsblk
            if isinstance(search,(Table,Column,list,tuple,set)):
                search = str(search[0])
            if self.refresh_cache:
                # bypass DB
                if clsname == 'MIS':
                    getDB(search,delete=True)
                return getDCS(search,*args,**kwargs)
            res = getDB(search, *args, **kwargs)
            if not res:
                # search not found in DB
                return getDCS(search,*args,**kwargs)
            
        else:
            # many aors/plans
            if self.refresh_cache:
                # bypass DB, but don't perform query yet
                if clsname == 'MIS':
                    getDB(search,delete=True)
                cfiles = [getDCS(s,raw=True) for s in search]
            res = getDB(search,*args,**kwargs)
            if not res:
                cfiles = [getDCS(search,raw=True) for s in search]
                res = getDB(search,*args,**kwargs)
        return res

    def _proc_res(self, res, cls, *args, **kwargs):
        """Unified results from DB"""

        if kwargs.get('sql'):
            return res

        sql = str(res)
        
        rows = list(res.dicts())
        rows = rows if rows else None
        if rows is None:
            return None

        for row in rows:
            row['sql'] = sql

        if kwargs.get('raw'):
            fnames = {row['FILENAME'] for row in rows}
            if len(fnames) == 1:
                return fnames.pop()
            else:
                return fnames

        if kwargs.get('as_table'):
            return cls.as_table(rows)

        if kwargs.get('as_pandas'):
            return cls.as_pandas(rows)
            
        if kwargs.get('as_json'):
            return cls.as_json(rows)
            
        return rows

        

    def getPOS(self, search, *args, **kwargs):
        '''
        if isinstance(search,str):
            # single aor/plan/obsblk
            if self.refresh_cache:
                # bypass DB
                return self._getPOS(search, *args, **kwargs)

            pos = DCS._query_POS_table(search, *args, **kwargs)
            if not pos:
                return self._getPOS(planID, aorID=aorID, guide=guide, raw=raw,
                                    as_table=as_table,as_pandas=as_pandas,as_json=as_json)
            
        else:
            # many aors/plans
            if self.refresh_cache:
                # bypass DB, but don't perform query yet
                planIDs = (_aorID_to_planID(a) for a in aorID)
                cfiles = [self._getPOS(planID, raw=True) for planID in planIDs]

            pos = DCS._query_POS_table(aorID=aorID,guide=guide)
            if not pos:
                planIDs = (_aorID_to_planID(a) for a in aorID)
                cfiles = [self._getObsPlan(planID, raw=True) for planID in planIDs]
                pos = DCS._query_POS_table(aorID=aorID, guide=guide)
            pos = pos if pos else None

        if raw:
            pos = list(pos)
            fnames = {p['FILENAME'] for p in pos}
            if len(fnames) == 1:
                return fnames.pop()
            else:
                return fnames

        if as_table:
            return POS.as_table(pos)

        if as_pandas:
            return POS.as_pandas(pos)
            
        if as_json:
            return POS.as_json(pos)
        return pos
        '''
        if isinstance(search,str):
            search = _aorID_to_planID(search)
        else:
            search = [_aorID_to_planID(s) for s in search]
        # search must be a planID
        return self._get(search, 'POS', *args, **kwargs)
        

    def getFlightPlan(self, search, *args, **kwargs):
        """Get MIS dict"""
        '''
        if kwargs.get('utctab'):
            # must get .mis file
            return self._getFlightPlan(search,*args,**kwargs)
        '''
        '''
        if self.refresh_cache:
            # bypass DB
            return self._getFlightPlan(flightid,ObsBlk=ObsBlk,local=local,
                                       as_table=as_table,as_pandas=as_pandas, as_json=as_json)

        legs = DCS._query_MIS_table(flightid, ObsBlk=ObsBlk)
        if not legs:
            return self._getFlightPlan(flightid,ObsBlk=ObsBlk,local=local,
                                       as_table=as_table,as_pandas=as_pandas,as_json=as_json)

        if as_table:
            return MIS.as_table(legs)
        
        if as_pandas:
            return MIS.as_pandas(legs)
            
        if as_json:
            return MIS.as_json(legs)
        return legs
        '''
        return self._get(search, 'MIS', *args, **kwargs)

    def getFlightSeries(self, series, get_ids=False, as_json=True, local=None):
        flightids = self._getFlightIDs_in_Series(series,local=local)
        if get_ids:
            return flightids

        mis = [self.getFlightPlan(flightid, as_json=as_json, local=local) for flightid in flightds]
        return mis

    
    def _query_AOR_table(self, search, *args, **kwargs):
        if isinstance(search,str):
            search = [search]

        '''
        # query multiple, and allow for planIDs, ObsBlks
        if kwargs.get('pos'):
            # BROKEN!!!
            aors = AOR.select(AOR,POS.POSName).join(POS).where((AOR.aorID.in_(search))  |
                                                               (AOR.planID.in_(search)) |
                                                               (AOR.ObsBlk.in_(search))).order_by(AOR.ObsBlk,AOR.order,AOR.aorID)
        elif kwargs.get('guide'):
            print('sup')
            print(search)
            aors = AOR.select().where((AOR.aorID.in_(search))  |
                                      (AOR.planID.in_(search)) |
                                      (AOR.ObsBlk.in_(search))).order_by(AOR.ObsBlk,AOR.order,AOR.aorID)
            aors = list(aors)
            print(list(aors[0].guide.execute()))
            exit()
        else:
            aors = AOR.select().where((AOR.aorID.in_(search))  |
                                      (AOR.planID.in_(search)) |
                                      (AOR.ObsBlk.in_(search))).order_by(AOR.ObsBlk,AOR.order,AOR.aorID)
        '''
        if kwargs.get('delete'):
            q = AOR.delete().where((AOR.aorID.in_(search))  |
                                   (AOR.planID.in_(search)) |
                                   (AOR.ObsBlkID.in_(search)))
            q.execute()
            return None
        
        aors = AOR.select().where((AOR.aorID.in_(search))  |
                                  (AOR.planID.in_(search)) |
                                  (AOR.ObsBlkID.in_(search))).order_by(AOR.ObsBlkID,AOR.order,AOR.aorID)
        return self._proc_res(aors, AOR, *args, **kwargs)

    def _query_GUIDE_table(self, search, *args, **kwargs):
        if isinstance(search,str):
            search = [search]

        if kwargs.get('delete'):
            q = GUIDE.delete().where((GUIDE.aorID.in_(search)) |
                                     (GUIDE.planID.in_(search)) |
                                     (GUIDE.ObsBlkID.in_(search)))
            q.execute()
            
        guides = GUIDE.select().where((GUIDE.aorID.in_(search)) |
                                      (GUIDE.planID.in_(search)) |
                                      (GUIDE.ObsBlkID.in_(search))).order_by(GUIDE.Radius)
        return self._proc_res(guides, GUIDE, *args, **kwargs)


    def _query_FAOR_table(self, search, *args, **kwargs):
        if isinstance(search,str):
            search = [search]

        if kwargs.get('delete'):
            q = FAOR.delete().where(FAOR.AORID.in_(search))
            q.execute()
            return None
            
        faors = FAOR.select().where(FAOR.AORID.in_(search))
        return self._proc_res(faors,FAOR, *args, **kwargs)

    def _query_SCT_table(self, search, *args, **kwargs):
        if isinstance(search,str):
            search = [search]

        if kwargs.get('delete'):
            q = SCT.delete().where(SCT.AORID.in_(search))
            q.execute()
            return None
            
        scts = SCT.select().where(SCT.AORID.in_(search))
        return self._proc_res(scts, SCT, *args, **kwargs)

    def _query_AOR_table_old(self,aorID=None,planID=None,ObsBlkID=None,pos=False,guide=False):
        """Perform DB lookup on AOR table"""
        if aorID:
            if isinstance(aorID,str):
                # single aorID
                if _is_aorID(aorID):
                    try:
                        aor = AOR.get(AOR.aorID==aorID)
                        aors = aor.__data__

                        if pos:
                            # get POSname
                            prows = self.getPOS(aorID)
                            if prows:
                                try:
                                    prow = list(filter(lambda x:not x['isNod'],prows))[0]
                                    posname = prow['POSName']
                                except IndexError:
                                    posname = ''
                            else:
                                posname = ''
                        aors['POSname'] = posname
                                    
                            
                        if guide:
                            # attach guide stars to aors
                            grows = self.getPOS(aorID,guide=guide)
                            if grows:
                                aors['GUIDESTARS'] = grows
                    except RuntimeError:
                        aors = None
                elif _is_planID(aorID):
                    return self._query_AOR_table(planID=aorID, guide=guide, pos=pos)
                elif 'OB' in aorID:
                    return self._query_AOR_table(ObsBlkID=aorID, guide=guide, pos=pos)
                else:
                    return None
                    
            else:
                # query multiple, and allow for planIDs, ObsBlkIDs
                aors = AOR.select().where((AOR.aorID.in_(aorID))  |
                                          (AOR.planID.in_(aorID)) |
                                          (AOR.ObsBlkID.in_(aorID))).order_by(AOR.ObsBlkID,AOR.order,AOR.aorID).dicts()
                #aors = [aor.__data__ for aor in aors]
                aors = list(aors)
                aors = aors if aors else None
                if pos:
                    # get POSname
                    prows = self.getPOS(aorID)#,guide=guide)
                    if prows:
                        prow = list(filter(lambda x:not x['isNod'],prows))
                        for aor in aors:
                            p = list(filter(lambda x:x['AORID'] == aor['aorID'],prow))
                            if p:
                                posname = p[0]['POSName']
                            else:
                                posname = ''
                            aor['POSName'] = posname
                    else:
                        for aor in aors:
                            aor['POSName'] = ''

                if guide:
                    # attach guide stars to aors
                    grows = self.getPOS(aorID,guide=guide)
                    if grows:
                        for aor in aors:
                            g = list(filter(lambda x:x['AORID'] == aor['aorID'],grows))
                            if g:
                                aor['GUIDESTARS'] = g
                
        elif planID:
            if isinstance(planID,str):
                aors = AOR.select().where(AOR.planID==planID).order_by(AOR.ObsBlkID,AOR.order,AOR.aorID).dicts()
            else:
                aors = AOR.select().where(AOR.planID.in_(planID)).order_by(AOR.ObsBlkID,AOR.order,AOR.aorID).dicts()
            aors = list(aors)
            aors = aors if aors else None

            if pos:
                # get POSname
                prows = self.getPOS(planID)
                if prows:
                    prow = list(filter(lambda x:not x['isNod'],prows))
                    for aor in aors:
                        p = list(filter(lambda x:x['AORID'] == aor['aorID'],prow))
                        if p:
                            posname = p[0]['POSName']
                        else:
                            posname = ''
                        aor['POSName'] = posname
                else:
                    for aor in aors:
                        aor['POSName'] = ''

            if guide:
                # attach guide stars to aors
                grows = self.getPOS(planID,guide=guide)
                if grows:
                    for aor in aors:
                        g = list(filter(lambda x:x['AORID'] == aor['aorID'],grows))
                        if g:
                            aor['GUIDESTARS'] = g
            
        elif ObsBlkID:
            if isinstance(ObsBlkID,str):
                aors = AOR.select().where(AOR.ObsBlkID==ObsBlkID).order_by(AOR.ObsBlkID,AOR.order,AOR.aorID).dicts()
            else:
                aors = AOR.select().where(AOR.ObsBlkID.in_(ObsBlkID)).order_by(AOR.ObsBlkID,AOR.order,AOR.aorID).dicts()
            aors = list(aors)
            aors = aors if aors else None
            
            if pos:
                # get POSname
                prows = self.getPOS(ObsBlkID)
                if prows:
                    prow = list(filter(lambda x:not x['isNod'],prows))
                    for aor in aors:
                        p = list(filter(lambda x:x['AORID'] == aor['aorID'],prow))
                        if p:
                            posname = p[0]['POSName']
                        else:
                            posname = ''
                        aor['POSName'] = posname
                else:
                    for aor in aors:
                        aor['POSName'] = ''

            if guide:
                # attach guide stars to aors
                grows = self.getPOS(ObsBlkID,guide=guide)
                if grows:
                    for aor in aors:
                        g = list(filter(lambda x:x['AORID'] == aor['aorID'],grows))
                        if g:
                            aor['GUIDESTARS'] = g
            
        else:
            aors = None
        return aors

    def _query_MIS_table(self, search, *args, **kwargs):
        if isinstance(search,str):
            search = [search]
        '''
        if flightid:
            if isinstance(flightid,str):
                # single flightid
                legs = MIS.select().where(MIS.FlightPlan==flightid).order_by(MIS.Leg).dicts()
            else:
                legs = MIS.select().where(MIS.FlightPlan.in_(flightid)).order_by(MIS.Leg).dicts()
            #legs = [leg.__data__ for leg in legs]
            legs = list(legs)
            legs = legs if legs else None
                
        elif ObsBlkID:
            legs = MIS.select().where(MIS.ObsBlkID==ObsBlkID).order_by(MIS.FlightPlan,MIS.Leg).dicts()
            #legs = [leg.__data__ for leg in legs]
            legs = list(legs)
            legs = legs if legs else None
        elif flightname:
            if isinstance(flightname,str):
                # single flightname
                legs = MIS.select().where(MIS.FlightName==flightname).order_by(MIS.FlightPlan,MIS.Leg).dicts()
            else:
                legs = MIS.select().where(MIS.FlightPlan.in_(flightname)).order_by(MIS.FlightPlan,MIS.Leg).dicts()
            #legs = [leg.__data__ for leg in legs]
            legs = list(legs)
            legs = legs if legs else None
        else:
            legs = None
        return legs
        '''
        if kwargs.get('delete'):
            q = MIS.delete().where((MIS.FlightPlan.in_(search)) |
                                   (MIS.FlightName.in_(search)) |
                                   (MIS.ObsBlkID.in_(search)))
            q.execute()
            return None
        
        legs = MIS.select().where((MIS.FlightPlan.in_(search)) |
                                  (MIS.FlightName.in_(search)) |
                                  (MIS.ObsBlkID.in_(search))).order_by(MIS.FlightPlan,MIS.Leg)
        return self._proc_res(legs, MIS, *args, **kwargs)

    def _query_POS_table(self, search, guide, *args, **kwargs):
        if isinstance(search,str):
            search = [search]

        if guide:
            pos = POS.select(POS,GUIDE).join(GUIDE).switch(POS).where((POS.AORID.in_(search)) | (POS.planID.in_(search)))
        else:
            pos = POS.select().where((POS.AORID.in_(search)) | (POS.planID.in_(search)))


        return self._proc_res(pos,POS, *args, **kwargs)
        
            
            

    def _getObsPlan(self, search,
                    rel='observationPlanning/DbProxy/getObsPlanAORs.jsp',
                    form='DIRECT',
                    raw=False,
                    as_table=False,
                    as_pandas=False,
                    as_json=False,
                    sql=False,
                    insert=True):
        """Download AOR xml, and store in database."""
        planID = _aorID_to_planID(search)
        query = (self.dcsurl/rel).with_query({'origin':'GI','obsplanID':planID})
        cfile = self._queryDCS(query,form)

        aorcfg = self.mcfg['AOR']
        guidecfg = self.mcfg['GUIDE']
        
        rows = AOR.to_rows(cfile, aorcfg)
        rows = self._add_missing_obsblks(rows)
        if rows and insert:
            AOR.replace_rows(self.db, rows)

        grows = GUIDE.to_rows(cfile, guidecfg)
        if grows and insert:
            GUIDE.replace_rows(self.db, grows)
            
        if raw:
            return cfile

        aors = self._query_AOR_table(search, sql=sql, raw=raw,
                                     as_table=as_table,as_json=as_json,as_pandas=as_pandas)
        return aors
        
    def _getObsPlanXML(self, planID,
                      rel='observationPlanning/observingPlanDetails.jsp',
                      form='form[action="DbProxy/getObsPlan.jsp"]',
                      raw=False):
        '''Download obsplan xml'''
        query = (self.dcsurl/rel).with_query({'obspID':planID})
        cfile = self._queryDCS(query, form)

        if raw:
            return cfile
        
        return cfile

    def _getFlightIDs_in_Series(self, series,
                                rel='flightPlan/flightPlanSearch.jsp',
                                form='form[action="flightPlanSearch.jsp"]',
                                submit={'name':'flightNum'},
                                local=None):
        '''Get names of flight in series'''
        if local:
            lfiles = list(Path(local).glob('**/*.misxml'))
            if not lfiles:
                lfiles = list(Path(local).glob('**/*.mis'))
                
            flightids = [f.stem.replace('_INIT','') for f in lfiles]
            flightids = [f.replace('_SCI','') for f in flightids]
            flightids = [f.replace('_MOPS','') for f in flightids]
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
        
    def _getFlightPlan(self, flightid,
                       rel='flightPlan/downloadFlightPlanFile.jsp',
                       form='DIRECT',
                       local=None,
                       sql=False,
                       raw=False,
                       ObsBlkID=None,
                       as_table=False,
                       as_pandas=False,
                       as_json=False,
                       force_mis=False,
                       insert=True):
        '''Download mis file'''
        if local:
            try:
                #cfile = list(Path(local).glob('%s*.mis'%flightid))[0]
                cfile = list(Path(local).glob('**/%s*.misxml'%flightid))
                if not cfile:
                    cfile = list(Path(local).glob('%s*.mis'%flightid))[0]
                else:
                    cfile = cfile[0]
            except IndexError:
                cfile = Path(local)
            if raw:
                return local
            '''
            if utctab:
                _,tab = dcsMIS.get_legs(cfile)
                return tab
            else:
                #tab = Table.read(cfile,format='mis-tab')
                #return tab
                pass
            '''
                

        else:
            ftype = 'mis' if force_mis else 'misxml'
            query = (self.dcsurl/rel).with_query({'fpid':flightid,'fileType':ftype})
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

        '''
        if utctab:
            _,tab = dcsMIS.get_legs(cfile)
            return tab
        '''

        miscfg = self.mcfg['MIS']
        rows = MIS.to_rows(cfile, miscfg)
        if not rows:
            return None
        if insert:
            MIS.replace_rows(self.db, rows)

        if raw:
            return cfile

        legs = self._query_MIS_table(search=flightid, ObsBlkID=ObsBlkID,
                                     sql=sql, raw=raw,
                                     as_table=as_table,as_json=as_json,as_pandas=as_pandas)

        return legs


    def _getPOS(self,planID,
                rel='observationPlanning/SaveTargetPos.jsp',
                form='DIRECT',
                guide=False,
                raw=True,
                as_table=False,
                as_pandas=False,
                as_json=False,
                sql=False,
                insert=True):
        '''Download pos file and store in both POS and GUIDE tables'''

        formfill = {'cycleID':'-1',
                    'instrument':'ALL',
                    'planID':planID,
                    'resPerPage':'500',
                    'spectel1':'ALL',
                    'spectel2':'ALL',
                    'mode':'ALL',
                    'targetType':'ALL',
                    'inputEquinox':'2000',
                    'aorState':'ALL',
                    'posFormat':'MCCS'}
        
        query = (self.dcsurl/rel).with_query(formfill)
        cfile = self._queryDCS(query, form)

        poscfg = self.mcfg['POS']
        guidecfg = self.mcfg['GUIDE']

        prows = POS.to_rows(cfile, poscfg)
        if prows and insert:
            POS.replace_rows(self.db, prows)

        #grows = GUIDE.to_rows(cfile, guidecfg)
        #if grows and insert:
        #    GUIDE.replace_rows(self.db, grows)

        if not prows:
            return None

        if raw: 
            return cfile

        #pos = DCS._query_POS_table(aorID=aorID, planID=planID, guide=guide)
        pos = self._query_POS_table(planID, guide=guide, sql=sql, raw=raw,
                                    as_table=as_table,as_json=as_json,as_pandas=as_pandas)
        return pos

    def _getPOSOld(self,flightid,
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
        try:
            cfile = self._queryDCS(self.dcsurl/fname, form, 'TAR')
        except UnboundLocalError:
            cfile = None
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
    
    def getObsBlockSearch(self,blockid=None,planid=None,
                          rel='flightPlan/observationBlockSearch.jsp',
                          form='DIRECT',
                          raw=False):
        '''Download ObsBlock Leg Summary'''
        if blockid:
            query = (self.dcsurl/rel).with_query({'obsBlkIDs':blockid,
                                                  'resPerPage':500,
                                                  'getObsBlkInfo':'Submit'})
        elif planid:
            query = (self.dcsurl/rel).with_query({'obsPlanID':planid,
                                                  'resPerPage':500,
                                                  'cycleID':-1,
                                                  'instrument':'ALL',
                                                  'fpStatus':'ALL',
                                                  'inputEquinox':2000,
                                                  'getObsBlkInfo':'Submit'})
        else:
            return None
        cfile = self._queryDCS(query,form)

        if raw:
            return cfile
        
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

    #d = DCS(database=True,refresh_cache=False)
    #d.getAORs('07_0149')

    #exit()
    #d = DCS()
    #t = d.getAORSearch('HAWC',cycle=-1)
    #print(t)
    #print(t.colnames)

    d = DCS()
    #d.getObsBlockSearch('OB_06_0058_03')
    print(d.getObsBlockSearch(planid='90_0072'))

    #d = DCS(refresh_cache=True)
    #d.getAORDetails('07_0012_1')
    #d.getPOSBundle('201906_FO_DEBORAH')
    #
    #d = DCS(refresh_cache=True)
    #print(d.getAOR('90_0071',raw=True))
