#!/usr/bin/env python

####
#### By Elim Thompson (09/27/2018)
####
#### This is a script to pickle DRAGON level 6 
#### event (numu, nue, nutau, muon, data) from
#### i3 files.
####
#### This is needed because the hdf5 files from
#### Philipp already has level 6 level cuts
#### applied, which can be used to plot variables
#### before level 6 cuts are applied
####
#### This file only pickle corridor.
#### 
#### command to run:
#### $ python pickle_L6.py
####        --outdir <output folder>
####        --ptype  'genie_nue/mu/tau or data_201X or muon'
####
##########################################################

#### import packages
from __future__ import print_function
from optparse import OptionParser
from glob import glob
import cPickle, sys
import numpy as np

#### import icecube packages
from icecube import dataclasses, icetray, dataio
from icecube.dataclasses import I3Particle
from icecube.dataio import I3File

#### import corsika weighter packages
##   GaisserH4a is assumed here
##   H3a and H4a are more or less the same at ~ GeV
from icecube import weighting
from icecube.weighting.fluxes import GaisserH4a
from icecube.weighting.weighting import from_simprod
##  flux model
muflux = GaisserH4a()
##  generator
generator = from_simprod (11058)
generator *= 32000.

#### import neutrino weighter packages
from icecube import NuFlux, prob3
sys.path.append ('../../analysistools/')
from probmap import ProbMap
##   honda flux
nuflux = NuFlux.makeFlux ("IPhonda2014_spl_solmin")
##   osc prob 
params = {'dm21': 7.49e-5,
          'dm31': 2.526e-3,
          'theta23': np.arcsin (np.sqrt (0.440))  ,
          'theta12': np.arcsin (np.sqrt (0.308))  ,
          'theta13': np.arcsin (np.sqrt (0.02163)),
          'deltacp': 0                             }

####################################################
#### define where the i3 files are
####################################################
indir   = "/data/ana/LE/published/DRAGON_3y/level5p/"
datadir = "/data/user/elims/dragon_data/"

folder = {'genie_nue'  :sorted (glob (indir+"nue/12600/*"  )),
          'genie_numu' :sorted (glob (indir+"numu/14600/*" )),
          'genie_nutau':sorted (glob (indir+"nutau/16600/*")),
          'data_2012'  :sorted (glob (datadir+'/IC86_2/*'  )), 
          'data_2013'  :sorted (glob (datadir+'/IC86_3/*'  )), 
          'data_2014'  :sorted (glob (datadir+'/IC86_4/*'  )), 
          'muon'       :sorted (glob (datadir+'/corsika/11058/*/*'))}

####################################################
#### define variables to be pickled
####
#### {keyname:[frameobj, name]}
#### level6_bool and weight will be added
####################################################
corridor_key = 'IC86_Dunkman_L6_CorridorDOMs'
dragon_livetime = 2.5 * 3600. * 365.24 * 24

variables = {'corridor'   :[], 'weight'     :[]}

####################################################
#### function to append events
####################################################
def get_data_weight (frame):
    
    return 1/dragon_livetime

def get_corsika_weight (frame):

    ### ignore if no weight map
    if not frame.Has ('CorsikaWeightMap'): return -1

    pritype = [frame['CorsikaWeightMap']['PrimaryType']]
    energy  = [frame['CorsikaWeightMap']['PrimaryEnergy']]
    return muflux(energy, pritype)/generator(energy, pritype)

def get_osc_weight (frame, nfiles, probmap):

    ### ignore if no tree / weightdict
    if not frame.Has ('I3MCTree'): return -1
    if not frame.Has ('I3MCWeightDict'): return -1

    ### define frame objc
    primary    = frame['I3MCTree'][0]
    weightdict = frame['I3MCWeightDict']

    ### define neutrino vs anti-neutrino
    sign = np.sign (primary.pdg_encoding)
    if sign > 0:
        ## neutrinos
        nue  = I3Particle.ParticleType.NuE
        numu = I3Particle.ParticleType.NuMu
        ratio = 0.7
    else:
        ## anti-neutrinos
        nue  = I3Particle.ParticleType.NuEBar
        numu = I3Particle.ParticleType.NuMuBar
        ratio = 1-0.7

    ### get nue and numu fluxes
    coszen = np.cos (primary.dir.zenith)
    energy = primary.energy
    nueflux  = nuflux.getFlux (nue, energy, coszen)
    numuflux = nuflux.getFlux (numu, energy, coszen)

    ### get oscillation prob
    prob = probmap.get_prob ([energy], [coszen],
                             [primary.pdg_encoding],
                             params)

    ### get norm factor
    one_weight = weightdict['OneWeight']
    nsimevents = weightdict['NEvents']
    norm = one_weight / nsimevents / ratio / nfiles

    ### total weight
    atm_nue  = nueflux * prob[:,0] * norm
    atm_numu = numuflux * prob[:,1] * norm
    return atm_nue + atm_numu

def get_weight (frame, dtype, nfiles, probmap=None):
    
    if 'data' in dtype:
        return get_data_weight (frame)

    if 'muon' in dtype:
        return get_corsika_weight (frame)

    ## remaining must be neutrino
    return get_osc_weight (frame, nfiles, probmap)

def append (frame, dtype, events, nfiles, probmap=None):

    ### remove if not have L6 frame object
    if not frame.Has (corridor_key): return events
    
    ### get weights
    weight = get_weight (frame, dtype, nfiles,
                         probmap=probmap)
    if weight<0: return events
    events['weight'].append (weight)
    
    ### get corridor
    corridor = frame[corridor_key].value
    events['corridor'].append (corridor)

    return events

####################################################
#### function to define pass ? and level6_bool
####################################################
def passed (frame):

    ### pass L3 ?
    if not frame.Has ("IC86_Dunkman_L3"): return False
    if not frame["IC86_Dunkman_L3"].value: return False

    ### pass L4 ?
    if not frame.Has ("IC86_Dunkman_L4"): return False
    if not frame['IC86_Dunkman_L4']['result']==1.: return False

    if frame.Has ("IC86_Dunkman_L6"): 
        santa = frame['IC86_Dunkman_L6']['santa_direct_charge']
    elif frame.Has ('IC86_Dunkman_L6_SANTA_DirectDOMs'):
        santa = frame ['IC86_Dunkman_L6_SANTA_DirectDOMs'].value
    else: return False
    if santa<3.: return False

    ### pass L5 ?
    if not frame.Has ("IC86_Dunkman_L5"): return False
    if frame['IC86_Dunkman_L5']['bdt_score'] < 0.2: return False
    return True

####################################################
#### function to collect events
####################################################
def collect_events (dtype):

    ### files to be looped through
    filenames = folder[dtype]
    nfiles    = float (len (filenames))
    print ('nfiles: {0}'.format (nfiles))

    ### initialize an event holder
    events = {key:[] for key in variables.keys ()}

    ### get prob map if genie
    ### also oscillate NC by default
    probmap = ProbMap (matter=True, params=params) \
              if 'genie' in dtype else None

    ### loop through files
    for filename in filenames:
        ## update me
        print ('+-- {0}'.format (filename))
        ## open file
        i3file = I3File (filename, 'r')
        ## loop through frames
        while i3file.more ():
            try:
                frame = i3file.pop_physics ()
                # append when physics frame
                if frame:
                    if passed (frame):
                        events = append (frame, dtype, events,
                                         nfiles,
                                         probmap=probmap)
                else:
                    print ('| this frame is empty ...')
                    continue
            except:
                print ('| this file is busted :/')
                pass
        ## close file
        print ('| unweighted nevents: {0}'.format (len (events['weight'])))
        print ('| weighted   nevents: {0}'.format (np.sum (events['weight'])))
        i3file.close()

    ### numpify
    events = {key:np.array (value)
              for key, value in events.items ()}

    ### update me
    weight = events['weight']
    print ('+-------------------------------------')
    print ('| unweighted before L6cuts : {0}'.format (len (weight)))
    print ('| weighted   before L6cuts : {0:4f} mHz'.format (np.sum (weight)*1000.))

    ### return output
    return events

###########################################
#### main function
###########################################
if __name__ == '__main__':

    #### parser
    usage = "%prog [options] <inputfiles>"
    parser = OptionParser(usage=usage)
    parser.add_option ("--outdir", type="string", default='~/',
                       help = "out directory for plotting")
    parser.add_option ("--ptype", type="string", default='genie_nue',
                       help = "genie_nue/mu/tau or corsika or data_2012/3/4")
    (options, args) = parser.parse_args()
    ptype   = options.ptype
    outfile = options.outdir + 'level6_'+ptype+'.pckl'

    ### loop through events
    events = collect_events (ptype)

    ### save it
    with open (outfile, 'wb') as f:
        cPickle.dump (events, f, protocol=2)
    f.close ()
