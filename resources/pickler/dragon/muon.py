#!/usr/bin/env python

####
#### By Elim Thompson (09/25/2018)
####
#### This script convert DRAGON muon HDF5
#### files to pickled files with selection
#### variables at L4, 5, 6 and a pre-calculated
#### weight.
####
#### L4 and L5 are from corsika files
#### L6 is from icc data
##################################################

#### import packages
from __future__ import print_function
import tables, cPickle
import numpy as np

### for corsika weighting
##  GaisserH4a is assumed here
##  H3a and H4a are more or less the same
from icecube import weighting
from icecube.weighting.fluxes import GaisserH4a
from icecube.weighting.weighting import from_simprod

############################################
#### define output
############################################
outdir = '/data/user/elims/for_nutau_paper/resources/pickled_data/dragon/preweighted/'

############################################
#### define data
############################################
### where the files are
corsika = '/home/fhuang/scratch/fhuang/corsika/Matt_L3_corsika_11058.hdf5'
iccdata = '/data/ana/LE/unblinded/DRAGON_NuTau_appearance/muons_icc/Matt_L5b_icc_data_IC86_2_3_4_new_reco_def2.hdf5'

### open them
cor = tables.open_file (corsika, 'r')
icc = tables.open_file (iccdata, 'r')

############################################
#### calculate all the weights (may as well)
############################################
### corsika weighting
##  flux model
flux = GaisserH4a()
##  generator
generator = from_simprod(11058)
generator *= 32000.
##  truth info
energy = np.array (cor.root.CorsikaWeightMap.cols.PrimaryEnergy[:])
ptype  = np.array (cor.root.CorsikaWeightMap.cols.PrimaryType[:])
##  get weights
cor_weights = flux(energy, ptype)/generator(energy, ptype)

### icc data weighting
##  all icc background are scaled by 0.146
dragon_livetime = 2.5 * 3600. * 365.24 * 24.
icc_length  = len (icc.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.energy)
icc_weights = np.ones (icc_length) / 0.146 / dragon_livetime

############################################
#### Start pickling L4, 5, 6
############################################
### quick function to numpify and apply cut
func = lambda array, cut: np.array (array[:])[cut.astype (bool)]

############################################
#### pickle L4 from corsika
####
#### 1. when pickling L4, L4 cut is not
####    applied. L3 cut is already applied.
#### 2. santa hits was in L3, but now moved
####    to L4.
############################################
### define L3 cut (i.e. all)
l3cut = np.ones (cor.root.IC86_Dunkman_L4.nrows)

### define level 4 boolean
L4result = func (cor.root.IC86_Dunkman_L4.cols.result, l3cut).astype (bool)
L4santa  = func (cor.root.IC86_Dunkman_L6_SANTA_DirectDOMs.cols.value, l3cut) >= 3.0
level4_bool = np.logical_and (L4result, L4santa)

### selected variables/events
muonL4 = {'rt_fid_charge': func (cor.root.IC86_Dunkman_L4.cols.rt_fid_charge, l3cut),
          'num_hit_doms' : func (cor.root.IC86_Dunkman_L4.cols.num_hit_doms , l3cut),
          'santa_hits'   : func (cor.root.IC86_Dunkman_L6_SANTA_DirectDOMs.cols.value, l3cut),
          'level4_bool'  : level4_bool, 
          'weight'       : func (cor_weights, l3cut) }

### store output
with open (outdir + 'level4_muon.pckl', 'wb') as f:
    cPickle.dump (muonL4, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L4rate = muonL4['weight'][muonL4['level4_bool'].astype (bool)]
L4rate = round (np.sum (L4rate)*1000., 0)
print ('+-------------------------------------')
print ('| After L4, muon rate is {0} mHz'.format (L4rate))

############################################
#### pickle L5 from corsika
####
#### 1. when pickling L5, L5 cut is not
####    applied. L4 cut is applied.
############################################
### define L4 cuts
l4cut = muonL4['level4_bool']

### define level 5 boolean
level5_bool = func (cor.root.IC86_Dunkman_L5.cols.bdt_score   , l4cut) >= 0.2

### selected variables/events
muonL5 = {'separation'  : func (cor.root.IC86_Dunkman_L5.cols.separation  , l4cut),
          'total_charge': func (cor.root.IC86_Dunkman_L5.cols.total_charge, l4cut),
          'bdt_score'   : func (cor.root.IC86_Dunkman_L5.cols.bdt_score   , l4cut),
          'level5_bool' : level5_bool.astype (float),
          'weight'      : func (cor_weights, l4cut) }

### store output
with open (outdir + 'level5_muon.pckl', 'wb') as f:
    cPickle.dump (muonL5, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L5rate = muonL5['weight'][muonL5['level5_bool'].astype (bool)]
L5rate = round (np.sum (L5rate)*1000., 4)
print ('| After L5, muon rate is {0} mHz'.format (L5rate))

############################################
#### pickle L6 from icc data
####
#### 1. when pickling L6, L6 cut is not
####    applied. L4 and L5 cuts are applied.
############################################
### define L4 and L5 cuts
L4result = icc.root.IC86_Dunkman_L4.cols.result[:].astype (bool)
L4santa  = icc.root.IC86_Dunkman_L6_SANTA_DirectDOMs.cols.value >= 3.0
L4vich   = icc.root.IC86_Dunkman_L4.cols.result_invertedVICH[:].astype (bool)
l4cut = np.logical_and (L4santa, np.logical_or (L4result, L4vich))
l5cut = np.logical_and (l4cut, icc.root.IC86_Dunkman_L5.cols.bdt_score[:] >= 0.2)

### define level 6 boolean
l6startcontain = icc.root.IC86_Dunkman_L6.cols.mn_start_contained[:].astype (bool)
l6stopcontain  = icc.root.IC86_Dunkman_L6.cols.mn_stop_contained[:].astype (bool)
l6contain      = np.logical_and (l6startcontain, l6stopcontain)
l6corridor     = icc.root.IC86_Dunkman_L6.cols.corridor_doms_over_threshold[:]>=2
level6_bool    = np.logical_and (l6contain, l6corridor)

### define vertex position
startx = icc.root.IC86_Dunkman_L6.cols.x_prime[:]
starty = icc.root.IC86_Dunkman_L6.cols.y_prime[:]
startrho = np.sqrt (startx**2, starty**2)

stopx = icc.root.IC86_Dunkman_L6.cols.stop_x_prime[:]
stopy = icc.root.IC86_Dunkman_L6.cols.stop_y_prime[:]
stoprho = np.sqrt (stopx**2, stopy**2)

### selected variables/events
muonL6 = {'start_rho'   : func (startrho, l5cut),
          'start_z'     : func (icc.root.IC86_Dunkman_L6.cols.z_prime[:]-350., l5cut),
          'stop_rho'    : func (stoprho, l5cut),
          'stop_z'      : func (icc.root.IC86_Dunkman_L6.cols.stop_z_prime[:]-350., l5cut),
          'deltaLLH'    : func (icc.root.IC86_Dunkman_L6.cols.delta_LLH, l5cut),
          'reco_logen'  : np.log10 (func (icc.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.energy, l5cut)),
          'reco_coszen' : np.cos (func (icc.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.zenith, l5cut)),
          'reco_x'      : func (icc.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.x     , l5cut),
          'reco_y'      : func (icc.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.y     , l5cut),
          'reco_z'      : func (icc.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.z     , l5cut),
          'level6_bool' : func (level6_bool.astype (float), l5cut),
          'weight'      : func (icc_weights, l5cut) }

### store output
with open (outdir + 'level6_final_muon.pckl', 'wb') as f:
    cPickle.dump (muonL6, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L6rate = muonL6['weight'][muonL6['level6_bool'].astype (bool)]
L6rate = round (np.sum (L6rate)*1000., 4)
print ('| After L6, muon rate is {0} mHz'.format (L6rate))
print ('+-------------------------------------')

############################################
#### close files
############################################
cor.close ()
icc.close ()
