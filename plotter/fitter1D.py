#!/usr/bin/env python

####
#### By Elim Thompson (09/27/2018)
####
#### This script contains the fitter class
#### which performs a simple fit to a 1D
#### histogram with normalizations allowed
#### to float.
####
#### This is a very rough fitter, which
#### does not take into account systematics.
#### The purpose of this fitter is simply to
#### make 1D distributions at lower level
#### of event selection look better.
#### Only atmmu and numu and noise are allowed
#### to float.
####
############################################

from __future__ import print_function
from scipy.optimize import minimize
from copy import deepcopy
import numpy as np

############################################
#### Fitter1D class
############################################
class Fitter1D (object):

    def __init__ (self, events, scalar):

        ''' initialize fitter1D instance
        '''
        
        self._events  = events
        self._scalar  = scalar

        ### set parameters
        self._xranges = None
        self._nbins = None
        self._dtypes = None

        ### set norms
        self._norm_numu = 1.0
        self._norm_muon = 1.0
        self._norm_noise = 1.0

    @property
    def xranges (self):
        return self._xranges
        
    @xranges.setter
    def xranges (self, xranges):
        self._xranges = xranges

    @property
    def nbins (self):
        return self._nbins

    @nbins.setter
    def nbins (self, nbins):
        self._nbins = nbins

    @property
    def dtypes (self):
        return self._dtypes

    @dtypes.setter
    def dtypes (self, dtypes):
        self._dtypes = dtypes

    def _build_ahisto (self, array, weight):

        ### define boolean within histogram range
        selected = np.logical_and (array>=self.xranges[0],
                                   array<self.xranges[1] )

        ### select events within range
        xvalue = array [selected]
        weight = weight[selected] * self._scalar
        
        ### build histogram
        H = np.histogram (xvalue,
                          bins=self.nbins,
                          range=self.xranges,
                          weights=weight)[0]
        H2 = np.histogram (xvalue,
                           bins=self.nbins,
                           range=self.xranges,
                           weights=(weight)**2)[0]
        return H, H2

    def _build_all_histos (self, level, var):

        events = deepcopy (self._events[level])
        ### initialize holders
        counts, variances = {}, {}

        ### loop through each dtype
        for dtype, event in events.items ():
            ## skip data
            if 'data' in dtype: continue
            ## build this histogram
            hist, hist2 = self._build_ahisto (event[var],
                                              event['weight'])
            ## store
            counts[dtype] = hist
            variances[dtype] = hist2

        ### get mc dictionary
        self._mc = counts
        self._mc2 = variances
        ### get data hist
        self._data = self._build_ahisto (events['data'][var],
                                         events['data']['weight'])[0]
        return

    def _chi2 (self, norms):

        ### no negative norms
        for norm in norms:
            if norm < 0: return np.inf
        
        ### initialize holders
        counts    = np.zeros (self._nbins)
        variances = np.zeros (self._nbins)

        ### loop through each data type
        for dtype in self.dtypes:
            ## skip data
            if 'data' in dtype: continue
            ## get factor
            norm = norms[1] if 'muon' in dtype else \
                   norms[2] if 'noise' in dtype else \
                   norms[0]
            ## counts of this dtype
            count = norm * self._mc[dtype]
            counts += count
            ## variances of this dtype
            variance = norm**2 * self._mc2[dtype]
            variances += variance
        
        ### chi2 between data and mc
        chi2 = (counts - self._data)**2 / \
               (self._data.astype (float) + variances)
        chi2[~np.isfinite(chi2)] = 0.
        return np.sum (chi2)
    
    def _fit (self):

        ''' fit '''

        ### seed values
        norm_numu = self._norm_numu
        norm_muon = self._norm_muon
        norm_noise = self._norm_noise

        result = minimize (self._chi2,
                           [norm_numu, norm_muon, norm_noise],
                           options={'disp':False},
                           method = 'BFGS',
                           tol = 1e-16 )
        return result

    def fit_events (self, level, var):

        ### where self._mc / mc2 / data are defined
        ### for this variable at this level
        self._build_all_histos (level, var)

        ### perform fit
        result = self._fit ()
        norms = result.x
        
        ### initialize scalted events
        event = deepcopy (self._events[level])

        ### scale event weights
        for dtype in event.keys ():
            ## skip data
            if 'data' in dtype: continue
            norm = norms[1] if 'muon' in dtype else \
                   norms[2] if 'noise' in dtype else \
                   norms[0]
            event[dtype]['weight'] *= norm

        ### return event
        return event

