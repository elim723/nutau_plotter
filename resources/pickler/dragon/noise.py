#!/usr/bin/env python

####
#### By Elim Thompson (09/25/2018)
####
#### This script convert DRAGON noise *I3*
#### files to pickled files with selection
#### variables at L4, 5 and a pre-calculated
#### weight.
####
#### L6 noise is negligible.
##################################################

#### import packages
from __future__ import print_function
from icecube import dataio, dataclasses, icetray
from glob import glob
import numpy as np
import cPickle

############################################
#### define output
############################################
outdir = '/data/user/elims/for_nutau_paper/resources/pickled_data/dragon/preweighted/'

############################################
#### define data
############################################
### where the files are
i3files = sorted (glob ('/data/user/elims/dragon_data/noise/Level5nc_IC86.2012_VuvuzelaPureNoise.990015*'))
nfiles  = len (i3files)

### to count corrupted file
##  corrupted files are defined as the one
##  where processing failed due to segmentation
##  fault. Segmentation fault arises when it
##  comes to bad memory often due to conflicts
##  in software. 
stream = []

############################################
#### calculate all the weights 
############################################
aweight = (1-2800*30e-6)/(nfiles * 1000.0 * (100e-3-20e-6))

############################################
#### initialize holders
############################################
levels = {'level4': {'rt_fid_charge': [], 'num_hit_doms' : [], 'santa_hits'   : [],
                     'level4_bool'  : [], 'weight'       : []                      },
          'level5': {'separation'   : [], 'total_charge' : [], 'bdt_score'    : [],
                     'level5_bool'  : [], 'weight'       : []                      } }

############################################
#### Start pickling L4, 5, 6
############################################
### function to apply cuts
def apply_cuts (frame, level):
    ## must pass L3
    if not frame.Has ('IC86_Dunkman_L3'): return False
    passL3 = frame['IC86_Dunkman_L3'].value
    ## check L4
    if not frame.Has ('IC86_Dunkman_L4'): return False
    if not frame.Has ('IC86_Dunkman_L6_SANTA_DirectDOMs'): return False
    passL4 = frame['IC86_Dunkman_L4']['result']==1. and \
             frame['IC86_Dunkman_L6_SANTA_DirectDOMs'].value >= 3.0
    if level=='level4':
        print ('+--------------------')
        print ('| Level 4')
        #print ('| passL3   : {0}'.format (passL3))
        print ('| L4 result: {0}'.format (frame['IC86_Dunkman_L4']['result']))
        print ('| L4 santa : {0}'.format (frame['IC86_Dunkman_L6_SANTA_DirectDOMs'].value))
        print ('| passL4   : {0}'.format (passL4))
        return passL3 and passL4
    ## check L5
    if not frame.Has ('IC86_Dunkman_L5'): return False
    passL5 = frame['IC86_Dunkman_L5']['bdt_score'] > 0.2
    #print ('+--------------------')
    #print ('| Level 5')
    #print ('| passL3: {0}'.format (passL3))
    #print ('| passL4: {0}'.format (passL4))
    #print ('| passL5: {0}'.format (passL5))
    return passL3 and passL4 and passL5

### function to append level 4
def append_level4 (frame):
    
    ## variables to be pickled and plotted
    levels['level4']['rt_fid_charge'].append (frame['IC86_Dunkman_L4']['rt_fid_charge'])
    levels['level4']['num_hit_doms'].append (frame['IC86_Dunkman_L4']['num_hit_doms'])
    levels['level4']['santa_hits'].append   (frame['IC86_Dunkman_L6_SANTA_DirectDOMs'].value)
    ## if you get to this point, you passed level4
    levels['level4']['level4_bool'].append (True)
    ## weight is all the same
    levels['level4']['weight'].append (aweight)
    return 

### function to append level 5
def append_level5 (frame):
    
    ## variables to be pickled and plotted
    levels['level5']['separation'].append   (frame['IC86_Dunkman_L5']['separation'])
    levels['level5']['total_charge'].append (frame['IC86_Dunkman_L5']['total_charge'])
    levels['level5']['bdt_score'].append    (frame['IC86_Dunkman_L5']['bdt_score'])
    ## if you get to this point, you passed level5
    levels['level5']['level5_bool'].append (True)
    ## weight is all the same
    levels['level5']['weight'].append (aweight)
    return 

### function to append events
def append_event (frame, level):
    ## append level4
    if level=='level4':
        return append_level4 (frame)
    ## append level5
    if level=='level5':
        return append_level5 (frame)
    return

###################################################
#### looping over all of the runs
###################################################
index = 0
#### loop through each i3 file
for infile in i3files:
    print ('infile: {0}'.format (infile))
    try:
        ### open the file
        i3file = dataio.I3File (infile, 'r')
        while i3file.more ():
            frame = i3file.pop_physics ()
            if frame:
                ## if available physics frame
                ## loop through each level
                for level in levels.keys ():
                    if apply_cuts (frame, level):
                        # append this event if it
                        # passes the cut
                        append_event (frame, level)
        i3file.close ()
    except:
        ### if can't open the file,
        ### this file is corrupted.
        ### weights will be adjusted at the end
        stream.append (infile)
        print (' THIS FILE IS CORRUPTED ...')
        pass

    ### print rates
    line = '  pass L4, L5 rates: {0}, {1}'
    print (line.format (np.sum (levels['level4']['weight']),
                        np.sum (levels['level5']['weight']) ))

    index += 1
    #if index == 100: break

###################################################
#### rescale rate based on # corrupted files
###################################################
for level, ddict in levels.items ():
    print ('{0} weight: {1}'.format (level, ddict['weight']))
    if len (ddict['weight']) > 0:
        ddict['weight'] *= nfiles / float (nfiles - len (stream))

### print final rates
print (line.format (np.sum (levels['level4']['weight']),
                    np.sum (levels['level5']['weight']) ))

### print number of corrupted files
print ('number of corrupted files {0}'.format (len (stream)))

###################################################
#### save pickled files
###################################################
for level, ddict in levels.items ():
    with open (outdir + '/'+level+'_noise.p', 'wb') as f:
        cPickle.dump (ddict, f, protocol=2)
    f.close ()

