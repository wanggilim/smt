#! /usr/bin/env python
# chainmap + configparser combo
from collections import ChainMap
import astropy.units as u
import pickle
from configparser import ConfigParser
from copy import deepcopy

class DeepChainMap(ChainMap):
    """Variant of ChainMap that allows direct updates to inner scopes"""
    def __setitem__(self, key, value):
        for mapping in self.maps:
            if key in mapping:
                mapping[key] = value
                return
        self.maps[0][key] = value

    def __delitem__(self, key):
        for mapping in self.maps:
            if key in mapping:
                del mapping[key]
                return
        raise KeyError(key)

class CMap(ChainMap):
    """Creates a chainmap from cfg and sections"""
    def __init__(self,cfg,sections,**kwargs):
        self.cfg = cfg
        self.sections = sections

        maps = filter(lambda sec:sec in cfg,sections)
        maps = [cfg[sec] for sec in maps]

        if kwargs:
            maps.insert(0,kwargs) # set command-line args first
            self.kwargs = kwargs
        else:
            # insert blank dict at front
            maps.insert(0,{})
            
        super().__init__(*maps)

    def _get_conv(self, option, typ=str, fallback=None):
        if option in self:
            return typ(self[option])
        else:
            return fallback
        
    def getfloat(self,option,fallback=None):
        return self._get_conv(option,typ=float,fallback=fallback)
    def getint(self,option,fallback=None):
        return self._get_conv(option,typ=int,fallback=fallback)
    def getints(self,option,fallback=None):
        try:
            val = self.getint(option,fallback)
            return val
        except ValueError:
            split = self[option].split(' ')
            if len(split) > 1:
                return [int(x) for x in split]
            else:
                return [int(x)]
            
    def getboolean(self,option,fallback=None):
        if option in self:
            if str(self[option]).lower() in ('1','yes','true','on'):
                return True
            elif str(self[option]).lower() in ('0','no','false','off'):
                return False
            else:
                raise ValueError('boolean %s not understood' % self[option])
        else:
            return fallback
    def getquantity(self,option,unit='',fallback=None):
        try:
            val = self.getfloat(option,fallback)
            if val is not None:
                return u.Quantity(val,unit)
            else:
                return fallback
        except ValueError as e:
            split = self[option].split(' ')
            if len(split) > 1:
                return u.Quantity([float(x) for x in split],unit)
            else:
                return u.Quantity(split[0],unit)

    def has_option(self,option):
        return option in self

    def has_section(self,section):
        return section in self

    def add_map(self, mapping, last=True):
        """Add map.  If last is true, map is checked last.  If last is false, map is checked first."""
        if last:
            self.maps.append(mapping)
        else:
            self.maps.insert(0,mapping)

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        cls = self.__class__

        # pickle cfg and read it back in
        repcfg = pickle.dumps(self.cfg)
        newcfg = pickle.loads(repcfg)

        sections = deepcopy(self.sections, memo)

        if hasattr(self,'kwargs') and self.kwargs:
            kwargs = deepcopy(self.kwargs,memo)
            return cls(newcfg,sections,**kwargs)
        else:
            return cls(newcfg,sections)

    def as_dict(self):
        return dict(self.items())


    @classmethod
    def from_AOR_dict(cls, cfg, aor, **kwargs):
        """Make CMap from aor dictionary.  Sections are looked up with decreasing level of specificity."""

        # these are listed from last to first in look up sequence
        
        planid = aor['planID']                             # 07_0225
        blkid = aor['ObsBlkID']                            # OB_07_0225_01
        aorid = aor['aorID']                               # 07_0225_1        
        
        flight = aor['FlightName']                         # GAVIN
        mission = aor['FlightPlan']                        # 201910_FO_GAVIN

        flight_planid = '_'.join((flight,planid))          # GAVIN_07_0225
        mission_planid = '_'.join((mission,planid))        # 201910_FO_GAVIN_07_0225

        flight_blkid = '_'.join((flight,blkid))            # GAVIN_OB_07_0225_01
        mission_blkid = '_'.join((mission,blkid))           # 201910_FO_GAVIN_OB_07_0225_01

        flight_aorid = '_'.join((flight,aorid))            # GAVIN_07_0225_1
        mission_aorid = '_'.join((mission,aorid))          # 201910_FO_GAVIN_07_0225_1
        
        legflight = '_'.join(('Leg%i'%aor['Leg'],flight))  # Leg7_GAVIN
        legmission = aor['fkey']                           # Leg7_201910_FO_GAVIN
        
        legfplan = '_'.join((legflight,planid))            # Leg7_GAVIN_07_0225
        legmplan = '_'.join((legmission,planid))           # Leg7_201910_FO_GAVIN_07_0225

        legfblkid = '_'.join((legflight,blkid))            # Leg7_GAVIN_OB_07_0225_01
        legmblkid = '_'.join((legmission,blkid))           # Leg7_201910_FO_GAVIN_OB_07_0225_01
        
        legfaorid = '_'.join((legflight,aorid))            # Leg7_GAVIN_07_0225_1
        legmaorid = '_'.join((legmission,aorid))           # Leg7_201910_FO_GAVIN_07_0225_1        

        sections = (legmaorid,legfaorid, legmblkid,legfblkid, legmplan,legfplan, legmission,legflight,
                    mission_aorid,flight_aorid, mission_blkid,flight_blkid, mission_planid,flight_planid,
                    mission,flight, aorid,blkid,planid, 'DEFAULT')
        return cls(cfg, sections, **kwargs)
