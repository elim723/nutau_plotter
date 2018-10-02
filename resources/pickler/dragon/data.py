#!/usr/bin/env python

####
#### By Elim Thompson (09/25/2018)
####
#### This script convert DRAGON data HDF5
#### files to pickled files with selection
#### variables at L4, 5, 6 and a pre-calculated
#### weight.
####
##################################################

#### import packages
from __future__ import print_function
import tables, cPickle
import numpy as np

############################################
#### define output
############################################
outdir = '/data/user/elims/for_nutau_paper/resources/pickled_data/dragon/preweighted/'

############################################
#### define data
############################################
### where the files are
lower = '/data/user/fhuang/Matt_level4_data/Matt_L4_data_IC86_2_3_4_new_reco.hdf5'
final = '/data/user/peller/nutau/data/Matt_L5b_data_BS_excluded_IC86_2_3_4_new_reco.hdf5'

### open them
lower = tables.open_file (lower, 'r')
final = tables.open_file (final, 'r')

############################################
#### calculate all the weights 
############################################
dragon_livetime = 2.5 * 3600. * 365.24 * 24.
lower_length  = lower.root.IC86_Dunkman_L4.nrows
lower_weights = np.ones (lower_length) / dragon_livetime
final_length  = final.root.IC86_Dunkman_L4.nrows
final_weights = np.ones (final_length) / dragon_livetime

############################################
#### Start pickling L4, 5, 6
############################################
### quick function to numpify and apply cut
func = lambda array, cut: np.array (array[:])[cut.astype (bool)]

############################################
#### pickle L4 from lower
####
#### 1. when pickling L4, L4 cut is not
####    applied. L3 cut is already applied.
#### 2. santa hits was in L3, but now moved
####    to L4.
############################################
### define L3 cut (i.e. all)
l3cut = np.ones (lower.root.IC86_Dunkman_L4.nrows)

### define level 4 boolean
L4result = func (lower.root.IC86_Dunkman_L4.cols.result, l3cut).astype (bool)
L4santa  = func (lower.root.IC86_Dunkman_L6_SANTA_DirectDOMs.cols.value, l3cut) >= 3.0
level4_bool = np.logical_and (L4result, L4santa)

### selected variables/events
dataL4 = {'rt_fid_charge': func (lower.root.IC86_Dunkman_L4.cols.rt_fid_charge, l3cut),
          'num_hit_doms' : func (lower.root.IC86_Dunkman_L4.cols.num_hit_doms , l3cut),
          'santa_hits'   : func (lower.root.IC86_Dunkman_L6_SANTA_DirectDOMs.cols.value, l3cut),
          'level4_bool'  : level4_bool, 
          'weight'       : func (lower_weights, l3cut) }

### store output
with open (outdir + 'level4_data.pckl', 'wb') as f:
    cPickle.dump (dataL4, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L4rate = dataL4['weight'][dataL4['level4_bool'].astype (bool)]
L4rate = round (np.sum (L4rate)*1000., 0)
print ('+-------------------------------------')
print ('| After L4, data rate is {0} mHz'.format (L4rate))

############################################
#### pickle L5 from lower
####
#### 1. when pickling L5, L5 cut is not
####    applied. L4 cut is applied.
############################################
### define L4 cuts
l4cut = dataL4['level4_bool']

### define level 5 boolean
level5_bool = func (lower.root.IC86_Dunkman_L5.cols.bdt_score   , l4cut) > 0.2

### selected variables/events
dataL5 = {'separation'  : func (lower.root.IC86_Dunkman_L5.cols.separation  , l4cut),
          'total_charge': func (lower.root.IC86_Dunkman_L5.cols.total_charge, l4cut),
          'bdt_score'   : func (lower.root.IC86_Dunkman_L5.cols.bdt_score   , l4cut),
          'level5_bool' : level5_bool.astype (float),
          'weight'      : func (lower_weights, l4cut) }

### store output
with open (outdir + 'level5_data.pckl', 'wb') as f:
    cPickle.dump (dataL5, f, protocol=2)
f.close ()


### print rates to see if it makes sense
L5rate = dataL5['weight'][dataL5['level5_bool'].astype (bool)]
L5rate = round (np.sum (L5rate)*1000., 3)
print ('| After L5, data rate is {0} mHz'.format (L5rate))

############################################
#### pickle L6 from final
####
#### 1. when pickling L6, L6 cut is not
####    applied. L4 and L5 cuts are applied.
############################################
### define L4 and L5 cuts
L4result = final.root.IC86_Dunkman_L4.cols.result[:].astype (bool)
L4santa  = final.root.IC86_Dunkman_L6_SANTA_DirectDOMs.cols.value >= 3.0
l4cut = np.logical_and (L4result, L4santa)
l5cut = np.logical_and (l4cut, final.root.IC86_Dunkman_L5.cols.bdt_score[:] > 0.2)

### define level 6 boolean
l6startcontain = final.root.IC86_Dunkman_L6.cols.mn_start_contained[:].astype (bool)
l6stopcontain  = final.root.IC86_Dunkman_L6.cols.mn_stop_contained[:].astype (bool)
l6contain      = np.logical_and (l6startcontain, l6stopcontain)
l6corridor     = final.root.IC86_Dunkman_L6.cols.corridor_doms_over_threshold[:]<=1
level6_bool    = np.logical_and (l6contain, l6corridor)

### define vertex position
startx = final.root.IC86_Dunkman_L6.cols.x_prime[:]
starty = final.root.IC86_Dunkman_L6.cols.y_prime[:]
startrho = np.sqrt (startx**2, starty**2)

stopx = final.root.IC86_Dunkman_L6.cols.stop_x_prime[:]
stopy = final.root.IC86_Dunkman_L6.cols.stop_y_prime[:]
stoprho = np.sqrt (stopx**2, stopy**2)

### selected variables/events
dataL6 = {'start_rho'   : func (startrho, l5cut),
          'start_z'     : func (final.root.IC86_Dunkman_L6.cols.z_prime[:]-350., l5cut),
          'stop_rho'    : func (stoprho, l5cut),
          'stop_z'      : func (final.root.IC86_Dunkman_L6.cols.stop_z_prime[:]-350., l5cut),
          'deltaLLH'    : func (final.root.IC86_Dunkman_L6.cols.delta_LLH, l5cut),
          'energy'      : func (final.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.energy, l5cut),
          'zenith'      : func (final.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.zenith, l5cut),
          'x'           : func (final.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.x     , l5cut),
          'y'           : func (final.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.y     , l5cut),
          'z'           : func (final.root.IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC.cols.z     , l5cut),
          'level6_bool' : func (level6_bool.astype (float), l5cut),
          'weight'      : func (final_weights, l5cut) }

### store output
with open (outdir + 'level6_data.pckl', 'wb') as f:
    cPickle.dump (dataL6, f, protocol=2)
f.close ()

### print rates to see if it makes sense
L6rate = dataL6['weight'][dataL6['level6_bool'].astype (bool)]
L6rate = round (np.sum (L6rate)*1000., 3)
print ('| After L6, data rate is {0} mHz'.format (L6rate))
print ('+-------------------------------------')

############################################
#### close files
############################################
lower.close ()
final.close ()
