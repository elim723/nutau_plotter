#!/usr/bin/env python

####
#### By Elim Cheung (07/24/2018)
####
#### This script includes basic functions to
#### define a template object with MC and data
#### histograms.
####
#### See generate_template.py to see how this
#### class is used.
#### 
####################################################################

from __future__ import print_function
import numpy as np
import cPickle, os

from misc import Map, seconds_per_year
from member import Member
from library import Library

####################################################################
#### Template class
####################################################################
class Template (object):

    ''' A class to create a template instance by
             - getting user's inputs
             - loading baseline histograms from library
             - creating hyperplane instance from library
             - adding up contributions from all members
               to obtain a final MC template
             - creating data histogram from either pseudo
               or real data

        See `generate_template.py` for how to use `Template`.
    '''
    
    def __init__(self,
                 members         =None ,
                 edges           =None ,
                 pdictpath       =None ,
                 nuparam_textfile=None ,
                 matter          =True ,
                 oscnc           =False,
                 inverted        =False,
                 verbose         =1    ):

        ''' initialize a template instance

            :param members         (list): data types included in MC template
            :param edges           (dict): histogram bin edges for energy,
                                           zenith, and pid
                                           {'e': [5.62, 10, 17.78, 31.62, 56.23],
                                            'z': [0, 0.93, 1.98, 3.1415],
                                            'p': np.array ([0., 50., 1000.])     } 
            :param pdictpath        (str): path to the folder where dictionaries
                                           of data sets are located
            :param nuparam_textfile (str): full filename of nuisance parameter textfile
            :param matter          (bool): If True, include matter effect
            :param oscnc           (bool): If True, oscillate NC
            :param inverted        (bool): If True, inverted hierarchy is assumed
            :param verbose          (int): If 0: no print out
                                           If 1: print histogram cascade / track rates
                                           If 2: print histogram counts per bin
                                           If 3: print progress in Hplane / Likelihood
        '''
        
        print ('#### ##################################################')
        print ('#### ############# Generate Template  #################')
        print ('####')

        #### define properties
        ## histograms
        self._members = members
        self._edges   = edges
        ## paths
        self._ppath  = pdictpath
        self._nufile = nuparam_textfile
        ## oscprob
        self._matter   = matter
        self._oscnc    = oscnc
        self._inverted = inverted
        ## misc
        self._fitdata = None
        self._hplanes = None
        self._verbose = verbose
        
    ###################################################
    #### string representation for Template
    ###################################################
    def __repr__ (self):

        ''' how to create a template in code '''

        return "Template ()"

    ###################################################
    #### set/get state of Template
    ###################################################
    def __getstate__ (self):

        ''' get state for pickling '''
        
        return self.__dict__

    def __setstate__ (self, d):

        ''' set state for pickling '''

        self.__dict__ = d

    ###################################################
    #### internal function to print rates
    ###################################################
    def _print_rates (self, htype, histos, nyears=None):

        ''' print histogram rate

            If nyears is provided, scale histos by livetime
        '''

        ## return if verbose == 0
        if not self._verbose: return 
        ## print header
        print ('#### ##################################################')
        print ('#### ############## {0:8} histogram ################'.format (htype))
        print ('####')
        line = '####  {0:9} | {1:9} | {2:9} | {3:9} |'
        print (line.format ('dtypes'.center (9), 'nevents'.center (9),
                            'cascade'.center (9), 'track'.center (9) ))
        print ('#### {0}'.format ('-'*48))
        ## print line
        if htype in ['data', 'template']:
            self._print_line (htype, histos, nyears=nyears)
        else: ## per data type
            for dtype in self.members:
                self._print_line (dtype, histos[dtype],
                                  nyears=nyears)
        ## print end
        print ('#### {0}'.format ('-'*48))
        print ('####')
    
    def _print_line (self, dtype, histo, nyears=None):

        ''' print rate info of a given histogram

            If nyears is provided,
            scale histos by livetime
        '''

        line = '####  {0:9} | {1:5} {2:3} | {3:5} {4:3} | {5:5} {6:3} |'
        numbers = self._collect_numbers (histo, nyears=nyears)
        print (line.format (dtype.center (7),
                            int (numbers['nevents'][0]),
                            int (numbers['nevents'][1]),
                            int (numbers['cascade'][0]),
                            int (numbers['cascade'][1]),
                            int (numbers['track'][0]),
                            int (numbers['track'][1])  ))

    def _collect_numbers (self, histo, nyears=None):

        ''' collect rates (and error bars)

            :param histo   (dict): dictionary of histogram and variance
                                   {'H':[], 'H2':[]}
            :param nyears (float): number of years
                                   If provided, scale histos by livetime

            :return numbers (dict): total rates / cascade /
                                    track counts and variances
        '''

        factor = nyears * seconds_per_year if nyears else 1.0
        cascade  = np.sum (histo['H'][:,:,0]) * factor
        cascade2 = np.sqrt (np.sum (histo['H2'][:,:,0])) * factor
        track  = np.sum (histo['H'][:,:,1]) * factor
        track2 = np.sqrt (np.sum (histo['H2'][:,:,1])) * factor
        nevents  = cascade + track
        nevents2 = np.sqrt (cascade2**2 + track2**2)
        return {'nevents': (nevents, nevents2),
                'cascade': (cascade, cascade2),
                'track'  : (track  , track2  )}
        
    ###################################################
    #### misc functions
    ###################################################
    def get_ranges (self):

        ''' get ranges of histogram axes '''
        
        return {k:(v[0], v[-1]) for k, v in self.edges.items ()}

    ###################################################
    #### external function to get baseline histograms
    ###################################################
    def get_baseline_histograms (self, params):

        ''' obtain all baseline histograms for all data types

            :param params (dict): values of floating parameters
                                  for weight calculations

            :return  lib     (`Library`): a library instance with weighters set
            :retrun  bhistos (dict)     : baseline histograms from all members
        '''
        
        ### set up a library to handle baseline members
        lib = Library (self.members, self.ppath,
                       edges  =self.edges,
                       verbose=self.verbose)

        ### set up weighters for baseline members
        lib.set_weighters (params,
                           matter=self.matter,
                           oscnc =self.oscnc )
        ### store oscillation probability maps
        self.probmaps = lib.probmaps

        ### collect baseline histograms
        bhistos  = lib.collect_base_histograms (params)
        self._print_rates ('baseline' , bhistos,
                           nyears=params['nyears'])
        
        return lib, bhistos

    ###################################################
    #### external function to get MC template
    ###################################################
    def get_template (self, params, lib, bhistos):

        ''' obtain a template from all members

            :param params  (dict)     : values of floating parameters
                                        for weight calculations
            :param lib     (`Library`): a Library instance for
                                        member manipulation
            :param bhistos (dict)     : baseline histogram and variances
                                        from all members in template (Hz)
       
            :retrun template (dict): count and variance of final template
                                     histogram (in counts)
        '''

        ### apply hyperplanes
        mhistos = lib.apply_hplanes (bhistos,
                                     self.hplanes,
                                     params)
        ### apply normalization
        mhistos = lib.scale_histos (mhistos, params)
        
        ### sum up histograms and variances
        mc = np.array (sum ([ mhistos[dtype]['H']
                              for dtype in self.members ]))
        var = np.array (sum ([ mhistos[dtype]['H2']
                               for dtype in self.members ]))

        ### store and print information
        template = {'H':mc, 'H2':var}
        self._print_rates ('hplaned' , mhistos)
        self._print_rates ('template', template)
        return mhistos, template

    ###################################################
    #### external function to get data histogram
    ###################################################
    def get_data (self, params, fitdata, diff):

        ''' obtain data histogram

            :param params  (dict): values of floating parameters
                                   for weight calculation
            :param fitdata (bool): If True, data histogram from real data
            :param diff    (bool): If True, injected values are different
                                            from seeded;
                                            build pseudo data histogram

            :retrun H  (np.array): data histogram in counts
            :return H2 (np.array): variance of data histogram
        '''

        if fitdata:
            ## If fitdata, fit to real data
            data = Member ('data', self.ppath, ranges=self.get_ranges ())
            weights = data.get_weights (params)
            H, H2 = data.get_histogram (self.edges, weights=weights)
            ## scaled by factors (in counts)
            norm = seconds_per_year * params['nyears']
            H *= norm; H2 *= norm**2

        elif diff:
            ## If injected and seeded are different, build pseudo data histogram
            ## library and baseline histograms from the injected parameters
            lib, bhistos = self.get_baseline_histograms (params)
            ## get template with injected data
            mhisto, temp = self.get_template (params, lib, bhistos)
            H, H2 = temp['H'], temp['H2']

        else:
            ## If none of the above, copy mc template from baseline
            H, H2 = self.template['H'], self.template['H2']

        ## store and print data histogram
        self.dhisto = Map ({'H':H, 'H2':H2})
        self._print_rates ('data', self.dhisto)
        return self.dhisto

    ###################################################
    #### properties of Template
    ###################################################
    ### histogram relateted
    @property
    def members (self):
        return self._members

    @members.setter
    def members (self, members):
        self._members = members

    @property
    def edges (self):
        return self._edges

    @edges.setter
    def edges (self, edges):
        self._edges = edges

    @property
    def bhistos (self):
        return self._bhistos

    @bhistos.setter
    def bhistos (self, bhistos):
        self._bhistos = bhistos

    @property
    def mhistos (self):
        return self._mhistos

    @mhistos.setter
    def mhistos (self, mhistos):
        self._mhistos = mhistos

    @property
    def template (self):
        return self._template

    @template.setter
    def template (self, template):
        self._template = template
        
    ### path related
    @property
    def ppath (self):
        return self._ppath

    @ppath.setter
    def ppath (self, pdictpath):
        self._ppath = pdictpath

    @property
    def nufile (self):
        return self._nufile

    @nufile.setter
    def nufile (self, nuparam_textfile):
        self._nufile = nuparam_textfile
    
    ### oscprob related
    @property
    def matter (self):
        return self._matter

    @matter.setter
    def matter (self, matter):
        self._matter   = matter

    @property
    def oscnc (self):
        return self._oscnc

    @oscnc.setter
    def oscnc (self, oscnc):
        self._oscnc    = oscnc

    @property
    def inverted (self):
        return self._inverted

    @inverted.setter
    def inverted (self, inverted):
        self._inverted = inverted
    
    ### misc
    @property
    def hplanes (self):
        return self._hplanes

    @hplanes.setter
    def hplanes (self, hplanes):
        self._hplanes = hplanes
    
    @property
    def verbose (self):
        return self._verbose

    @verbose.setter
    def verbose (self, verbose):
        self._verbose = verbose

