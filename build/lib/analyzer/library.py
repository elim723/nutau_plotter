#!/usr/bin/env python

####
#### By Elim Cheung (07/24/2018)
####
#### This script contain `Library` class, whose main purpose
#### is to perform manipulations involving different members
#### in the MC template. The key functions of this class are:
####
####     - load baseline events
####     - set weighters (given floating parameter values)
####     - create baseline histograms
####     - create hyperplane for a given data type
####
####################################################################

from __future__ import print_function
from copy import deepcopy
import numpy as np
import time

from weightcalculator import WeightCalculator
from hyperplane import HyperPlane
from misc import get_sets, Map
from member import Member

####################################################################
#### import constants needed
####################################################################
from misc import datatypes, seconds_per_year, default_edges, discrete_parameters
from misc import default_ranges, default_nu_sysvalues, default_mu_sysvalues

####################################################################
#### Library class
####################################################################
class Library (object):

        ''' A class that can
                  - load baseline events
                  - set weighters (given floating parameter values)
                  - create baseline histograms
                  - create hyperplane for a given data type

            Example to create Library object
            --------------------------------
            In [0] : import library
            In [1] : dtype = 'numucc'
            In [2] : lib   = library.Library ([dtype], 'pickled_files/')

            Example to set weights and get baseline histograms
            --------------------------------------------------
            In [3] : import nuparams
            In [4] : nufile     = 'nuisance_textfiles/nuparams_template.txt'
            In [5] : nuParams   = nuparams.Nuparams (nufile, isinverted=False)
            In [6] : params     = nuParams.extract_params ('seeded')
            In [7] : lib.set_weighters (params, oscnc=False)
            In [8] : basehistos = lib.collect_base_histograms (params)
            NOTE: weighters for baseline members are stored in Library.
            
            Example to get all systematic histograms for a specific data type
            -----------------------------------------------------------------
            In [9] : hparams   = nuParams.get_hplaned_dparams (dtype)
            In [10]: syshistos = lib.collect_sys_histograms (dtype, params, hparams)
            NOTE: everytime collect_sys_histograms () is called, it creates
                  a one-time-use weighter for each member called.

            Example to get HyperPlane objects for all data types
            ----------------------------------------------------
            In [10]: hplanes = lib.get_hplanes (nuParams)
        '''

        def __init__ (self,
                      members  =None,
                      pdictpath=None,
                      edges    =default_edges,
                      isdragon =False,
                      verbose  =1):

                ''' Initialize a `Library` instance

                    :param dtypes    (list): all members included for templates
                    :param pdictpath (str) : path to all pickled dataset dictionaries
                    :param edges     (dict): template bin edges 
                    :param isdragon  (bool): If True, it is a DRAGON set.
                    :param verbose   (int) : If 0 or 1, no printout
                                             If 2, printout progress
                '''
                #### define properties
                self._members = members
                self._edges   = edges
                self._ppath   = pdictpath
                self._isdragon = isdragon
                self._verbose = verbose
                self._ranges = {k:(v[0], v[-1]) for k, v in self.edges.items ()}
                
                ## a library must have baseline sets
                self._baseline = self._collect_base_members ()

        ###################################################
        #### call and set/get state of Library
        ###################################################
        def __getstate__ (self):

                ''' get state for pickling '''

                return self.__dict__

        def __setstate__ (self, d):

                ''' set state for pickling '''

                self.__dict__ = d

        def __call__ (self, dtype):

                ''' return baseline event info in this library '''
                
                return self._baseline [dtype]

        ###################################################
        #### internal functions related to dataset IDs
        ###################################################
        def _check_setid (self, setid, default):

                ''' check if setid is the default ID

                    :param default (dict): keys:default values of discrete systematics
                    :param setid   (str) : systematic values separated
                                           by '_' of this set

                    :return isref  (bool): If True, this is a reference set
                '''
                
                return setid == self._get_setid (default)

        def _get_setid (self, sysvalues):

                ''' get the set ID for the given sysvalues dictionary 

                    :param sysvalues (dict): keys:values of this discrete systematic set

                    :return setid (str): systematic values separated by '_' of this set
                '''
                
                setid = ''
                for i, dp in enumerate (sorted (sysvalues)):
                        setid += str (float (sysvalues[dp]))
                        if not i == len (sysvalues)-1: setid += '_'
                return setid

        def _get_setvalues (self, set_value, dtype, hparam, default):

                ''' get the values for the specific discrete set

                    :param set_value (float): value of one of the discrete param sets
                    :param hparam    (str)  : name of the discrete systematics
                    :param dtype     (str)  : name of this member
                    :param default   (dict) : keys:default values of
                                              discrete systematics
                    :return setvalues (dict): keys:values of discrete
                                              systematics for the set
                    :return setid     (str) : systematic values separated by
                                              '_' of this set
                '''

                ### copy default values
                setvalues = deepcopy (default)
                ### change the given hparam by the set_value
                setvalues[hparam] = set_value
                ### modification for muons
                if 'muon' in dtype:
                        ## if absorption / scattering: oversizing = 3
                        if hparam in ['absorption', 'scattering']:
                                setvalues['oversizing'] = 3
                        ## if domeff / forward: holeice = 30
                        if hparam in ['domeff', 'forward']: setvalues['holeice'] = 30
                ### modification for off axis bulkice sets
                if hparam == 'absorption' and set_value in [0.929, 1.142]:
                        setvalues['scattering'] = set_value

                ### define ID for this set
                setid = self._get_setid (setvalues)
                return setid, setvalues

        ###################################################
        #### internal functions related to collect members
        ###################################################
        def _collect_base_members (self):

                ''' collect all baseline members

                    :return members (`Map`): dictionary with all `Members` objects
                '''
                
                members = Map ({})
                for dtype in self.members:
                        members[dtype] = Member (dtype, self.ppath,
                                                 ranges=self._ranges,
                                                 isdragon=self.isdragon, 
                                                 baseline=True)
                return members

        def _collect_sys_members (self, dtype, hparam, default, has_bulkice):

                ''' collect all set members for a specific data type
                    and a specific discrete parameter

                    :param hparam  (str): name of the discrete systematics
                    :param dtype   (str): name of the given member
                    :param default (dict): keys:default values of discrete systematics
                    :param has_bulkice (bool): If True, include bulk ice off axis points

                    :return members (`Map`): dictionary with all `Members` objects
                '''
                
                set_members = Map ({})
                sets = get_sets (dtype, hparam, has_bulkice=has_bulkice)
                for s in sets:
                        setid, setvalues = self._get_setvalues (s, dtype,
                                                                hparam, default)
                        isdef = self._check_setid (setid, default)
                        ## don't waste time if it is the default set and already stored
                        if setid in set_members:
                                if not defid == setid:
                                        message = 'Library:' + fname + \
                                                  ' : WARNING : setid (' + setid + \
                                                  ') already exist in the dictionary !'
                                        print ('{0}'.format (message))
                                continue
                        set_members [setid] = Member (dtype, self.ppath,
                                                      ranges=self._ranges,
                                                      baseline=False,
                                                      isdragon=self.isdragon, 
                                                      sysvalues=setvalues)
                return set_members

        ###################################################        
        #### functions related to merge nunc
        ###################################################
        def merge_nunc (self, dtype, histos, ncdtypes):

                ''' merge all nc members into one histogram

                    :param dtype    (str)  : name of this data type
                    :param histos   (`Map`): histograms of all discrete
                                             sets from all dtypes
                    :param ncdtypes (list) : names of all NC members

                    :return nchistos (`Map`): merged NC histograms
                '''
                
                ## return if not nunc
                if not 'nc' in dtype: return histos[dtype]

                ## massage nc
                nchistos = Map ({})
                for setid in histos[ncdtypes[0]]:
                        nchistos[setid] = {'H': sum ([ histos[ncdtype][setid]['H']
                                                       for ncdtype in ncdtypes ]),
                                           'H2': sum ([ histos[ncdtype][setid]['H2']
                                                        for ncdtype in ncdtypes ]) }
                return nchistos
        
        ###################################################
        #### functions to set weighters
        ###################################################
        def set_weighters (self, params,
                           matter=True, oscnc=False):

                ''' set weighters for each member into self.
                    weighters per data type (for all systematic sets)
                    probmaps per numu / nue / nutau (for both CC and
                    NC and all systematic sets)

                    :param params (dict): values of floating parameters
                    :param matter (bool): If True, include matter effect
                    :param oscnc  (bool): If True, oscillate NC events
                '''
                
                weighters = Map ({})
                pmaps = self.probmaps if hasattr (self, 'probmaps') else Map ({})
                for dtype in self.members:
                        pmap = pmaps [dtype[:-2]] if 'nu' in dtype and \
                               dtype[:-2] in pmaps \
                               else None
                        weighters[dtype]=self._baseline[dtype].get_weighter(params,
                                                                            matter=matter,
                                                                            oscnc=oscnc,
                                                                            pmap=pmap)
                        if 'nu' in dtype: pmaps [dtype[:-2]] = weighters[dtype].probmap 
                self.weighters = weighters
                self.probmaps  = pmaps

        ###################################################
        #### functions related to obtain histograms
        ###################################################
        def _get_histogram (self, member, params,
                            isbaseline=False, matter=True, oscnc=False):

                ''' get one histogram from one member

                    :param member (`Member`): a member instance
                    :param params     (dict): values of floating parameters
                    :param isbaseline (bool): If True , use self.weighters
                                              If False, redefine weighters for sys sets
                    :param matter     (bool): If True, include matter effect
                    :param oscnc      (bool): If True, oscillate NC events

                    :return hdict (dict): histogram from this member
                                          {'H'  = histogram weighted by weights,
                                           'H2' = variance of H                 }
                '''

                dtype = member.dtype
                if isbaseline:
                        ## BASELINE: weighters in self
                        if not hasattr (self, 'weighters'):
                                self.set_weighters (params, matter=matter, oscnc=oscnc)
                        weighter = self.weighters [dtype]
                else:
                        ## NOT BASELINE: one-time weighter
                        pmap = self.probmaps[dtype[:-2]] if 'nu' in dtype else None
                        weighter = member.get_weighter (params, matter=matter,
                                                        oscnc=oscnc, pmap=pmap)
                        
                weights = member.get_weights (params, weighter=weighter)
                H, H2 = member.get_histogram (self.edges, weights=weights)
                return Map ({'H':H, 'H2':H2})

        def _get_sys_histograms (self, dtype, params, hparams, default,
                                 matter=True, oscnc=False):

                ''' get all systematic histograms of a data type.
                    
                    :param dtype   (str) : name of the given member
                    :param params  (dict): values of floating parameters
                    :param hparams (list): discrete parameters included in hplane
                    :param default (dict): keys:default values of discrete systematics
                    :param matter  (bool): If True, include matter effect
                    :param oscnc   (bool): If True, oscillate NC events
                
                    :return histos (dict): all systematic histograms for this
                                           member for this discrete parameter
                '''
                
                histos = Map ({})
                has_bulkice = 'scattering' in hparams and 'absorption' in hparams
                for hp in sorted (hparams):
                        if self.verbose > 1: print ('####    -- {0}'.format (hp))
                        ## get member of this discrete param
                        members = self._collect_sys_members (dtype, hp, default,
                                                             has_bulkice)
                        ## get histogram from each set
                        for setid in sorted (members):
                                if self.verbose > 1:
                                        print ('####        -- {0}'.format (setid))
                                # check if it is default set
                                isdef = self._check_setid (setid, default)
                                # don't waste time on redoing default set
                                if isdef and setid in histos: continue
                                # get histogram of this setid
                                histos [setid] = self._get_histogram (members[setid],
                                                                      params,
                                                                      isbaseline=False, 
                                                                      matter=matter,
                                                                      oscnc=oscnc)
                                # define normalization factor
                                if isdef: norm = np.sum (histos[setid]['H'])
                                # apply normalization factor if coin set
                                if hp=='coin':
                                        hsum = np.sum (histos [setid]['H'])
                                        factor = norm/hsum if hsum else 1.0
                                        histos [setid]['H'] *= factor
                        if self.verbose > 1: print ('####')
                return histos
       
        def collect_base_histograms (self, params,
                                     matter=True, oscnc=False):

                ''' get all baseline histograms from all members

                    Note: You might want to have weighters set defined
                          otherwise, it will do it here.

                    :param params (dict): values of floating parameters
                    :param matter (bool): If True, include matter effect
                    :param oscnc  (bool): If True, oscillate NC events

                    :return histos (dict): all baseline histograms
                '''
                
                histos = Map ({})
                for dtype in self.members:
                        if self.verbose > 1: print ('####    -- {0}'.format (dtype))
                        member = self._baseline [dtype]
                        histos [dtype] = self._get_histogram (member, params,
                                                              isbaseline=True,
                                                              matter=matter,
                                                              oscnc=oscnc)
                if self.verbose > 1: print ('####')
                return histos

        def collect_sys_histograms (self, dtype, params, hparams,
                                    matter=True, oscnc=False):

                ''' collect all systematic histograms of a data type

                    Note: You might want to have weighters set defined
                          otherwise, it will do it here.

                    :param dtype   (str) : name of this member
                    :param params  (dict): values of floating parameters
                    :param hparams (list): discrete parameters included in hyperplane
                    :param matter  (bool): If True, include matter effect
                    :param oscnc   (bool): If True, oscillate NC events

                    :return histos (dict): all systematic histograms for this member
                '''

                default = default_mu_sysvalues if 'muon' in dtype else \
                          default_nu_sysvalues
                hparams = sorted ([ param for param in default if param in hparams ])
                return self._get_sys_histograms (dtype, params, hparams, default,
                                                   matter=matter, oscnc=oscnc)

        ###################################################
        #### functions related to hyperplane
        ###################################################
        def get_hplanes (self, nuparams, matter=True, oscnc=False, verbose=1):

                ''' collect hyperplane objects from all data types

                    Note: hyperplanes are based upon systematic histograms
                          weighted by the seeded MC values in nuparams

                    :param nuparams (`Nuparams`): user settings of floating parameters
                    :param matter   (bool)      : If True, include matter effect
                    :param oscnc    (bool)      : If True, oscillate NC events
                    :param verbose  (int)       : If 0, no printout
                                                  If 1, print out basic info
                                                  If 2, print out info within chi2 fit

                    :return hplanes (dict): hyperplane objects for all members
                '''
                
                ## collect systematic histograms
                params = nuparams.extract_params ('seeded')
                ## only for neutrinos and muons
                dtypes = [ dtype for dtype in self.members
                           if dtype[:2] in ['nu', 'mu'] ]
                ## collect systematic information
                syshistos = Map ({})
                if self.verbose > 1: print ('####  collecting sys histograms')
                for dtype in dtypes:
                        if self.verbose > 1: print ('####  -- {0}'.format (dtype))
                        hparams = nuparams.get_hplaned_dparams (dtype)
                        syshistos[dtype] = self.collect_sys_histograms (dtype, params,
                                                                        hparams,
                                                                        matter=matter,
                                                                        oscnc=oscnc)
                if self.verbose > 1: print ('####')
                        
                ## special treatment for nunc
                nc = sorted ([ dtype for dtype in dtypes if 'nc' in dtype ])
                if len (nc) > 0:
                        dtypes = sorted ([ dtype for dtype in dtypes
                                           if not 'nc' in dtype ] + ['nunc'])

                ## collect hyperplane objects
                hplanes = Map ({})
                ## get hplanes
                if self.verbose > 1: print ('####  collecting hplanes')
                for dtype in dtypes:
                        if self.verbose > 1: print ('####    -- {0}'.format (dtype))
                        hparams = nuparams.get_hplaned_dparams (dtype)
                        ## hyperplane is built only if at least one discrete parameter
                        if len (hparams) == 0:
                                hplanes[dtype] = None; continue
                        expparams = nuparams.get_exp_dparams (dtype)
                        histos = self.merge_nunc (dtype, syshistos, nc)
                        hplanes[dtype] = HyperPlane (dtype, histos, expparams,
                                                     hparams, verbose=verbose)
                if self.verbose > 1: print ('####')
                return hplanes        

        def apply_hplanes (self, bhistos, hplanes, params):

                ''' multiply hyperplane factor to each member template

                    :param bhistos (dict): all baseline histograms
                    :param hplanes (dict): hyperplane objects for all members
                    :param params  (dict): values of floating parameters

                    :return mhistos (dict): modified histograms
                '''

                mhistos = Map ({})
                for dtype in self.members:
                        hplane = hplanes['nunc'] if 'nc' in dtype else \
                                 hplanes[dtype] if dtype in hplanes else None
                        ## if no hyperplane, histo same base histo
                        if not hplane:
                                mhistos[dtype] = bhistos[dtype]; continue
                        factors = hplane.apply (params)
                        mhistos[dtype] = {'H' : bhistos[dtype]['H'] * factors,
                                          'H2': bhistos[dtype]['H2'] * factors**2}
                return mhistos

        def scale_histos (self, histos, params):

                ''' scale histograms by appropriate factors

                    :param histos (dict): histograms to be scaled
                    :param params (dict): values of floating parameters

                    :return shistos (dict): scaled histograms
                '''

                shistos = Map ({})
                for dtype in histos:
                        ## get the normalization factor
                        key = 'atmmu' if 'muon' in dtype else \
                              'noise' if 'noise' in dtype else 'numu'
                        norm = params['norm_nutau'] if 'nutau' in dtype else \
                               params['norm_nc'] if 'nc' in dtype else 1.
                        norm *= params['norm_'+key]*seconds_per_year*params['nyears']
                        ## scale this histogram
                        shistos[dtype] = {'H' : histos[dtype]['H'] * norm,
                                          'H2': histos[dtype]['H2'] * norm**2}
                return shistos
        
        ###################################################
        #### properties of Library
        ###################################################
        @property
        def members (self):
                return self._members

        @members.setter
        def members (self, members):
                self._members = members

        @property
        def isdragon (self):
                return self._isdragon

        @isdragon.setter
        def isdragon (self, isdragon):
                self._isdragon = isdragon
                
        @property
        def edges (self):
                return self._edges

        @edges.setter
        def edges (self, edges):
                self._edges = edges

        @property
        def ppath (self):
                return self._ppath

        @ppath.setter
        def ppath (self, pdictpath):
                self._ppath = pdictpath

        @property
        def verbose (self):
                return self._verbose
        
        @verbose.setter
        def verbose (self, verbose):
                self._verbose = verbose
