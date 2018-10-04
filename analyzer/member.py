#!/usr/bin/env python

####
#### By Elim Cheung (07/24/2018)
####
#### This script contain Member class whose purpose is to
#### do stuff that belongs to this member. It includes:
####
####    - loads / selects / stores events of a given memeber
####    - define weighter
####    - get weights
####    - create a histogram
####
####################################################################

from __future__ import print_function
import numpy as np
import os, cPickle, sys

from analyzer import misc, weightcalculator
from misc import Map, Toolbox

####################################################################
#### import constants needed
####################################################################
from analyzer.nuparams import discrete_parameters as dparams

from misc import datatypes, seconds_per_year, default_ranges, greco_nyears
from misc import default_nu_sysvalues, default_mu_sysvalues
toolbox = Toolbox ()

####################################################################
#### Member class
####################################################################
class Member (object):

        ''' A class to store all events for generating a MC template.

            Given whichever discrete parameters are included, this
            class create a library with all events needed.

            Example
            -------
            In [0]: from library import Member
            In [1]: pdictpath = '/data/condor_builds/users/elims/ezfits/clean_ezfit/pickled_files/'
            In [2]: numucc = Member ('numucc', pdictpath, baseline=True)
        '''

        def __init__ (self,
                      dtype    =None,
                      pdictpath=None,
                      ranges   =default_ranges,
                      baseline =False,
                      isdragon =False,
                      sysvalues={}):

            ''' Initialize a Member instance

                :param dtype     (str) : data type name
                                         nu*cc / nu*nc / noise / muon / data
                :param pdictpath (str) : path to all pickled dictionaries folder
                :param ranges    (dict): limits of template ranges
                                         events outside are removed
                :param baseline  (bool): If True, get baseline sets
                :param isdragon  (bool): If True, it is a DRAGON set
                :param sysvalues (dict): discrete systematic values:
                                         domeff / holeice / forward /
                                         absorption / scattering /
                                         coin / oversizing
            '''

            self._dtype     = dtype
            self._isbase    = baseline
            self._ppath     = pdictpath
            self._ranges    = ranges
            self._sysvalues = sysvalues
            self._isdragon  = isdragon
            
            self._pfile     = self._get_pfile ()
            self._events    = self._load_events ()

        ###################################################
        #### call and set/get state of Library
        ###################################################
        def __getstate__ (self):

                ''' get state for pickling '''

                return self.__dict__

        def __setstate__ (self, d):

                ''' set state for pickling '''

                self.__dict__ = d
            
        def __call__ (self):
                
                ''' return event dictionary of this member '''

                return self._events
        
        ###################################################
        #### functions to load / select events
        ###################################################
        def _get_pfile (self):

                ''' return pickled file path'''
                
                pname = self.dtype
                if 'data' in self.dtype:
                        pname += '.p'
                elif self.isbase:
                        pname += '_baseline.p'
                else:
                        forward = str (self.sysvalues['forward']).replace ('-', 'n')
                        if self.sysvalues['forward'] > 0: forward = 'p' + forward
                        extra = '_coin'   + str (self.sysvalues['coin']) if 'nu' in self.dtype else \
                                '_oversizing' + str (self.sysvalues['oversizing']) 
                        pname +=  '_domeff'  + str (self.sysvalues['domeff'])  + \
                                  '_holeice' + str (self.sysvalues['holeice']) + \
                                  '_forward' + forward + \
                                  '_absorption' + str (self.sysvalues['absorption']) + \
                                  '_scattering' + str (self.sysvalues['scattering']) + \
                                  extra + '.p'
                return self.ppath + '/' + pname
        
        def apply_cut (self, ddict, cut):

                ''' apply cut to a dictionary

                    :param ddict (`Map`): dictionary containing all events
                    :param cut   (np.array of bool): boolean to select events

                    :return cdict (`Map`): dictionary containing selected events
                '''
            
                cdict = Map({})
                for key in ddict.keys():
                        ### deal with arrays
                        cdict[key] = toolbox.chop (ddict[key], cut)
                        ### deal with arrays in dictionary
                        if toolbox.is_dict (ddict[key]):
                                cdict[key] = Map ({})
                                for skey in ddict[key].keys():
                                        cdict[key][skey] = toolbox.chop (ddict[key][skey], cut)
                return cdict

        def _select_dragon_events (self, ddict):

                e = ddict.reco.e
                bdt = ddict.L5.bdt_score if 'L5' in ddict else ddict.reco.bdt_score
                boolean = np.logical_and (bdt>0.2, np.logical_and (e>=0, e<60))
                boolean = np.logical_and (boolean, boolean>=-3.)
                return self.apply_cut (ddict, boolean)
        
        def _select_events (self, ddict):

                ''' define and apply final cuts

                    :param ddict (`Map`): dictionary containing all events
                
                    :return sdict (`Map`): dictionary containing all selected events
                '''

                if self.isdragon:
                        return self._select_dragon_events (ddict)
                
                ### define containment cuts
                recoX = np.array (ddict.reco.X)
                recoY = np.array (ddict.reco.Y)
                recoZ = np.array (ddict.reco.Z)
                rho = np.sqrt ((recoX-46.29)**2 + (recoY+34.88)**2)
                cut = recoZ < -230.
                cut *= rho < 140.
                vertex = cut*np.logical_or((recoZ+230.)/(rho-90)<-4.4, rho< 90.)
                boolean = np.logical_and (recoZ>-500, vertex)
                ### define misc cuts
                boolean = np.logical_and (boolean, ddict.hits.SRT_nCh<100.)
                boolean = np.logical_and (boolean, ddict.geo.charge_asym<0.85)
                ### define template range cuts
                e = np.array (ddict.reco.e)
                z = np.array (ddict.reco.z)
                p = np.array (ddict.reco.pid)
                ecut = np.logical_and (e>=self.ranges['e'][0], e<self.ranges['e'][1])
                zcut = np.logical_and (z>=self.ranges['z'][0], z<self.ranges['z'][1])
                pcut = np.logical_and (p>=self.ranges['p'][0], p<self.ranges['p'][1])
                tempcut = np.logical_and (np.logical_and (ecut, zcut), pcut)
                boolean = np.logical_and (boolean, tempcut)
                return self.apply_cut (ddict, boolean)
        
        def _load_events (self):

                ''' load events and apply standard cuts
            
                    :return sdict (`Map`): dictionary containing selected events
                '''

                try:
                        with open (self._pfile, "rb") as f:
                                ddict = cPickle.load (f)
                        f.close()
                        return self._select_events (ddict)
                except IOError:
                        raise IOError ('{0} does not exist ...'.format (self._pfile))

        ###################################################
        #### functions related to weights
        ###################################################
        def get_weighter (self, params, matter=True, oscnc=False, pmap=None):

                ''' calculate weights for each events (simulation only)

                    :param params (dict): values of floating parameters
                    :param matter (bool): If True, include matter effect
                    :param oscnc  (bool): If True, oscillate NC events   
                    :param pmap   (`ProbMap`): a ProbMap instance with maps of
                                               oscillation probabilities

                    :return weighter (`weightCalculator`): a weight calculater
                                                           instance for this dtype
                '''
                
                weighter = None if 'data' in self.dtype else \
                           weightcalculator.NoiseWeighter (params, self._events) if 'noise' in self.dtype else \
                           weightcalculator.MuonWeighter  (params, self._events) if 'muon'  in self.dtype else \
                           weightcalculator.NeutrinoWeighter (self.dtype, params, self._events,
                                                              matter=matter, oscnc=oscnc, pmap=pmap)
                return weighter

        def get_weights (self, params, weighter=None,
                         matter=True, oscnc=False, pmap=None):

                ''' get weights for each events (all in Hz)

                    :param params   (dict): values of floating parameters
                    :param weighter (`WeightCalculator`): If None, define one here
                    :param matter   (bool): If True, include matter effect
                    :param oscnc    (bool): If True, oscillate NC events   
                    :param pmap     (`ProbMap`): a ProbMap instance with maps of
                                                 oscillation probabilities

                    :return weights (np.array): weights for each event
                '''

                ## data weights
                if 'data' in self.dtype:
                        greco_livetime = seconds_per_year * greco_nyears
                        return np.ones (len (self._events.reco.e)) / greco_livetime

                ## define weighter if not parsed
                if not weighter:
                        weighter = self.get_weighter (params, matter=matter, oscnc=oscnc, pmap=pmap)

                ## get weights
                return weighter.reweight (params)
        
        ###################################################
        #### functions to get histogram
        ###################################################
        def get_histogram (self, edges, weights=[], params=None):

                ''' get histogram based on template

                    :param edges   (dict): histogram bin edges
                                           {'e': 10**np.linspace (0.75, 1.75, 9),
                                            'z': np.arccos(np.linspace(-1.,1.,11))[::-1],
                                            'p': np.array ([0., 50., 1000.])  }
                    :param weights (np.array): weight for each events
                                               If empty, params must not be None.
                    :param params  (dict): values of floating parameters

                    :return H  (np.array): histogram based on given edges (weight)
                            H2 (np.array): variance of H
                '''

                reco = self._events.reco
                e, z, pid = np.array (reco.e), np.array (reco.z), np.array (reco.pid)
                w = self.get_weights (params=params) if len (weights) == 0 else \
                                              np.array (weights)
        
                #### make sure event variables are finite
                finite = np.logical_and (np.logical_and (np.logical_and (np.isfinite (e), np.isfinite (z)),
                                                         np.isfinite (w)), np.isfinite (pid))
                
                #### define data/edges/ranges
                data = [ e[finite], z[finite], pid[finite] ]
                hedges = (edges['e'], edges['z'], edges['p'])
                ranges = [ [edges['e'][0], edges['e'][-1]],
                           [edges['z'][0], edges['z'][-1]],
                           [edges['p'][0], edges['p'][-1]] ]
                
                #### build histogram
                H, edge = np.histogramdd (np.array (data).T, hedges, range=ranges, weights=w[finite])
                H2, edge = np.histogramdd (np.array (data).T, hedges, range=ranges, weights=(w[finite])**2)
                
                return H, H2

        ###################################################
        #### properties of Member
        ###################################################
        @property
        def dtype (self):
                return self._dtype

        @dtype.setter
        def dtype (self, dtype):
                self._dtype = dtype
        
        @property
        def ranges (self):
                return self._ranges

        @ranges.setter
        def ranges (self, ranges):
                self._ranges = ranges

        @property
        def ppath (self):
                return self._ppath

        @ppath.setter
        def ppath (self, pdictpath):
                self._ppath = pdictpath

        @property
        def isbase (self):
                return self._isbase

        @isbase.setter
        def isbase (self, isbase):
                self._isbase = isbase

        @property
        def isdragon (self):
                return self._isdragon

        @isdragon.setter
        def isdragon (self, isdragon):
                self._isdragon = isdragon
                
        @property
        def sysvalues (self):
                return self._sysvalues

        @sysvalues.setter
        def sysvalues (self, sysvalues):
                default = default_mu_sysvalues if 'muon' in self.dtype else \
                          default_nu_sysvalues
                ### make sure all discrete parameters are in
                for param in default:
                        ## if not in sysvalues, set to default value
                        if param not in sysvalues:
                                sysvalues[param] = default[param]
                self._sysvalues = sysvalues
