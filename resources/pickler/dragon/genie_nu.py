#!/usr/bin/env python

####
#### By Elim Thompson (09/25/2018)
####
#### This script convert DRAGON genie nue HDF5
#### files to pickled files with selection
#### variables at L4, 5, 6 and a pre-calculated
#### weight.
####
#### Of course, the same script can be used for
#### other flavor. But, for laziness (and the
#### potential that no one would use this code),
#### this script is dedicated to nue only.
##################################################

#### import packages
from __future__ import print_function
from optparse import OptionParser
import tables, cPickle, sys
import numpy as np

### for genie weighting
from icecube import weighting
from icecube.weighting.fluxes import GaisserH4a
from icecube.weighting.weighting import from_simprod

############################################
### parse options 
############################################
usage = "usage: python %prog --flavor <nu>"
parser = OptionParser (usage=usage)
parser.add_option ('--flavor', type='str', default=None,
                   help = "either numu / nue / nutau")
(options, args) = parser.parse_args ()
flavor = options.flavor

############################################
#### define output
############################################
outdir = '/data/user/elims/for_nutau_paper/resources/pickled_data/dragon/preweighted/'

############################################
#### define data
############################################
### where the files are
nu = '/data/user/fhuang/Matt_level3_mc/events__deepcore__IC86__runs_126001-126003,146001-146003,166001-166003__proc_v5digit__unjoined_with_fluxes_L3.hdf5'
nul6 = '/data/user/peller/level5p_pisa_hdf5/events__deepcore__IC86__runs_126001-126003,146001-146003,166001-166003__proc_v5digit__unjoined_with_fluxes_GENIE_Barr.hdf5'

### open it
nu = tables.open_file (nu, 'r')
nul6 = tables.open_file (nul6, 'r')

### WHERE THIS FLAVOR IS DEFINED
flavor = flavor
nucc, nubarcc = getattr (nu.root, flavor).cc, getattr (nu.root, flavor+'_bar').cc
nunc, nubarnc = getattr (nu.root, flavor).nc, getattr (nu.root, flavor+'_bar').nc
nul6cc, nubarl6cc = getattr (nul6.root, flavor).cc, getattr (nul6.root, flavor+'_bar').cc
nul6nc, nubarl6nc = getattr (nul6.root, flavor).nc, getattr (nul6.root, flavor+'_bar').nc 

############################################
#### Start pickling L4, 5, 6
############################################
### quick function to numpify and apply cut
func = lambda array, cut: np.array (array[:])[cut.astype (bool)]

### a function to apply everything to
### all current type and nu vs anti-nu
def merge (key, isL6=False, cut=None):
    array = []
    ns = [nul6cc, nubarl6cc, nul6nc, nubarl6nc] if isL6 else \
         [nucc, nubarcc, nunc, nubarnc]
    for n in ns:
        array += list (getattr (n, key)[:])
    array = np.array (array)
    if cut==None: return array
    return array [cut.astype (bool)]

############################################
#### function to get oscillated weights
############################################
#### import prob3 related packages
sys.path.append ('../../analysistools/')
from probmap import ProbMap

#### default scillation parameters
params = {'dm21': 7.49e-5,
          'dm31': 2.526e-3,
          'theta23': np.arcsin (np.sqrt (0.440))  ,
          'theta12': np.arcsin (np.sqrt (0.308))  ,
          'theta13': np.arcsin (np.sqrt (0.02163)),
          'deltacp': 0                             }

#### define oscillation probability map
probmap = ProbMap (matter=True, params=params)

#### a function to get pdg once and for all
apdg = 12 if flavor=='nue'  else \
       14 if flavor=='numu' else \
       16
def get_pdg (isL6=False, cut=None):
    pdgs = []
    ns = [nul6cc, nubarl6cc, nul6nc, nubarl6nc] if isL6 else \
         [nucc, nubarcc, nunc, nubarnc]
    for i, n in enumerate (ns):
        ## set global pdg based on the flavor
        length = len (n.weighted_aeff)
        p = np.ones (length) * apdg
        ## if index = 1 / 3, its nubar
        if i in [1, 3]: p*=-1
        ## append (as a list) to pdgs
        pdgs += list (p)
    ### numpify
    pdgs = np.array (pdgs)
    if cut==None: return pdgs
    ### apply cut if available
    return pdgs [cut.astype (bool)]

#### function to get oscillation weights
def get_osc_weights (e, cz, pdg, fluxes, aeff):

    ### get oscillation probability
    prob = probmap.get_prob (e, cz, pdg, params)

    ### get the two components
    atm_nue  = fluxes[:,0] * aeff * prob[:,0]
    atm_numu = fluxes[:,1] * aeff * prob[:,1]

    ### get total weights
    return atm_nue + atm_numu

### define true info for weights
nue_flux  = merge ('neutrino_nue_flux')
numu_flux = merge ('neutrino_numu_flux')
fluxes    = np.array ([nue_flux, numu_flux]).T
aeff      = merge ('weighted_aeff')
energy    = merge ('true_energy')
coszen    = merge ('true_coszen')
ptypes    = get_pdg (isL6=False, cut=None)
weights   = get_osc_weights (energy, coszen, ptypes,
                             fluxes, aeff)

############################################
#### pickle L4 from nu
####
#### 1. when pickling L4, L4 cut is not
####    applied. L3 cut is already applied.
#### 2. santa hits was in L3, but now moved
####    to L4.
############################################
### define L3 cut (i.e. all)

### define level 4 boolean
L4result = merge ('dunkman_L4').astype (bool)
L4santa  = merge ('dunkman_direct_doms') >= 3.0
level4_bool = np.logical_and (L4result, L4santa)

### selected variables/events
nuL4 = {'rt_fid_charge': merge ('rt_fid_charge')      ,
        'num_hit_doms' : merge ('num_hit_doms')       ,
        'santa_hits'   : merge ('dunkman_direct_doms'),
        'level4_bool'  : level4_bool.astype (float)   ,  
        'weight'       : weights                       }

### store output
with open (outdir + 'level4_genie_'+flavor+'.pckl', 'wb') as f:
    cPickle.dump (nuL4, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L4rate = nuL4['weight'][nuL4['level4_bool'].astype (bool)]
L4rate = round (np.sum (L4rate)*1000., 4)
print ('+-------------------------------------')
print ('| After L4, {0} rate is {1} mHz'.format (flavor, L4rate))

############################################
#### pickle L5 from nu
####
#### 1. when pickling L5, L5 cut is not
####    applied. L4 cut is applied.
############################################
### define L4 cuts
l4cut = nuL4['level4_bool'].astype (bool)

### define level 5 boolean
level5_bool = merge ('dunkman_L5', cut=l4cut) > 0.2

### selected variables/events
nuL5 = {'separation'  : merge ('separation'  , cut=l4cut),
        'total_charge': merge ('total_charge', cut=l4cut),
        'bdt_score'   : merge ('dunkman_L5'  , cut=l4cut),
        'level5_bool' : level5_bool.astype (float)       ,
        'weight'      : func  (weights, cut=l4cut) }

### store output
with open (outdir + 'level5_genie_'+flavor+'.pckl', 'wb') as f:
    cPickle.dump (nuL5, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L5rate = nuL5['weight'][nuL5['level5_bool'].astype (bool)]
L5rate = round (np.sum (L5rate)*1000., 3)
print ('| After L5, {0} rate is {1} mHz'.format (flavor, L5rate))

############################################
#### pickle L6 from nuL6
####
#### 1. when pickling L6, L6 cut is not
####    applied. L4 and L5 cuts are applied.
############################################
### define true info for weights
nue_flux  = merge ('neutrino_nue_flux', isL6=True)
numu_flux = merge ('neutrino_numu_flux', isL6=True)
fluxes    = np.array ([nue_flux, numu_flux]).T
aeff      = merge ('weighted_aeff', isL6=True)
energy    = merge ('true_energy', isL6=True)
coszen    = merge ('true_coszen', isL6=True)
ptypes    = get_pdg (isL6=True, cut=None)
weights   = get_osc_weights (energy, coszen, ptypes,
                             fluxes, aeff)

### define L4 and L5 cuts 
##  nuL6 doesn't have L4 result
##  safely assume all events passed L4 result
# L4result = merge ('dunkman_L4', isL6=True).astype (bool) 
l4cut = merge ('santa_direct_doms', isL6=True) >= 3.0
l5cut = merge ('dunkman_L5', isL6=True) > 0.2
l5cut = np.logical_and (l4cut, l5cut)

### define level 6 boolean
l6startcontain = merge ('mn_start_contained', isL6=True, cut=l5cut).astype (bool)
l6stopcontain  = merge ('mn_stop_contained' , isL6=True, cut=l5cut).astype (bool)
l6contain      = np.logical_and (l6startcontain, l6stopcontain)
## nuL6 does not have L6 corridor
#  assume all neutrino pass the corridor cut
l6corridor     = np.ones (len (l6contain)).astype (bool)
level6_bool    = np.logical_and (l6contain, l6corridor)

### define vertex position
startx = merge ('x_prime', isL6=True, cut=l5cut)
starty = merge ('y_prime', isL6=True, cut=l5cut)
startrho = np.sqrt (startx**2, starty**2)

stopx = merge ('stop_x_prime', isL6=True, cut=l5cut)
stopy = merge ('stop_y_prime', isL6=True, cut=l5cut)
stoprho = np.sqrt (stopx**2, stopy**2)

### selected variables/events
nuL6 = {'start_rho'   : startrho                                               ,
        'start_z'     : merge ('z_prime', isL6=True, cut=l5cut)-350.           ,
        'stop_rho'    : stoprho                                                ,
        'stop_z'      : merge ('stop_z_prime', isL6=True, cut=l5cut)-350.      ,
        'deltaLLH'    : merge ('pid', isL6=True, cut=l5cut)                    ,
        'energy'      : merge ('reco_energy', isL6=True, cut=l5cut)            ,
        'zenith'      : np.arccos (merge ('reco_coszen', isL6=True, cut=l5cut)),
        'x'           : merge ('reco_x', isL6=True, cut=l5cut)                 ,
        'y'           : merge ('reco_y', isL6=True, cut=l5cut)                 ,
        'z'           : merge ('reco_z', isL6=True, cut=l5cut)                 ,
        'level6_bool' : level6_bool.astype (float)                             ,
        'weight'      : func  (weights, cut=l5cut)                              }

### store output
with open (outdir + 'level6_genie_'+flavor+'.pckl', 'wb') as f:
    cPickle.dump (nuL6, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L6rate = nuL6['weight'][nuL6['level6_bool'].astype (bool)]
L6rate = round (np.sum (L6rate)*1000., 3)
print ('| After L6, {0} rate is {1} mHz'.format (flavor, L6rate))
print ('+-------------------------------------')

############################################
#### close files
############################################
nu.close ()
nul6.close ()
