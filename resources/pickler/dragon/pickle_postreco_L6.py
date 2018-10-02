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
####        --ptype  'genie_nue/mu/tau or data_201X or muon_201X'
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

#### import neutrino weighter packages
from icecube import NuFlux, prob3
sys.path.append ('../../analysistools/')
from probmap import ProbMap
##   honda flux
nuflux = NuFlux.makeFlux("IPhonda2014_spl_solmin")
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

folder = {'genie_nue'  :sorted (glob (indir+"nue/12600/*"  )),
          'genie_numu' :sorted (glob (indir+"numu/14600/*" )),
          'genie_nutau':sorted (glob (indir+"nutau/16600/*")),
          'data_2012'  :sorted (glob (indir+'data/IC86_2/*'  )), 
          'data_2013'  :sorted (glob (indir+'data/IC86_3/*'  )), 
          'data_2014'  :sorted (glob (indir+'data/IC86_4/*'  )), 
          'muon_2012'  :sorted (glob (indir+'data/IC86_2/*'  )), 
          'muon_2013'  :sorted (glob (indir+'data/IC86_3/*'  )), 
          'muon_2014'  :sorted (glob (indir+'data/IC86_4/*'  )) }

####################################################
#### define variables to be pickled
####
#### {keyname:[frameobj, name]}
#### level6_bool and weight will be added
####################################################
corridor_key = 'IC86_Dunkman_L6_CorridorDOMs'
reco_key     = 'IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC'
L6_key       = 'IC86_Dunkman_L6'
dragon_livetime = 2.5 * 3600. * 365.24 * 24

variables = {'corridor'   :[], 'weight'     :[],
             'start_z'    :[], 'start_rho'  :[],
             'stop_z'     :[], 'stop_rho'   :[],
             'reco_z'     :[], 'reco_rho'   :[],
             'reco_logen' :[], 'reco_coszen':[],
             'logen'      :[], 'coszen'     :[],
             'deltaLLH'   :[], 'level6_bool':[] }

####################################################
#### function to get weight
####################################################
def get_data_weight (frame, icc=False):
    
    weight = 1/dragon_livetime
    if icc: weight /= 0.146
    return weight

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

def get_weight (frame, dtype, nfiles, probmap=None, icc=False):
    
    if 'data' in dtype:
        return get_data_weight (frame, icc=False)

    if 'muon' in dtype:
        return get_data_weight (frame, icc=icc)

    ## remaining must be neutrino
    return get_osc_weight (frame, nfiles, probmap)

####################################################
#### function to append events
####################################################
def append (frame, dtype, events, nfiles, probmap=None, icc=False):

    ### get weights
    weight = get_weight (frame, dtype, nfiles,
                         probmap=probmap, icc=icc)
    if weight<0: return events
    events['weight'].append (weight)
    
    ### get corridor
    corridor = frame[L6_key]['corridor_doms_over_threshold']
    events['corridor'].append (corridor)

    ### get delta_LLH
    deltaLLH = frame[L6_key]['delta_LLH']
    events['deltaLLH'].append (deltaLLH)

    ### get level6 bool
    level6_bool = frame[L6_key]['result']
    events['level6_bool'].append (level6_bool)

    ### get reco info
    reco_logen  = np.log10 (frame[reco_key].energy)
    reco_coszen = np.cos   (frame[reco_key].dir.zenith)
    reco_x      = frame[reco_key].pos.x
    reco_y      = frame[reco_key].pos.x
    reco_rho    = np.sqrt ((reco_x-46.29)**2 + (reco_y+34.88)**2)
    reco_z      = frame[reco_key].pos.z
    events['reco_logen'].append  (reco_logen)
    events['reco_coszen'].append (reco_coszen)
    events['reco_rho'].append    (reco_rho)
    events['reco_z'].append      (reco_z)

    ### get containment
    start_x   = frame[L6_key]['x_prime']
    start_y   = frame[L6_key]['y_prime']
    start_rho = np.sqrt (start_x**2 + start_y**2)
    stop_x   = frame[L6_key]['stop_x_prime']
    stop_y   = frame[L6_key]['stop_y_prime']
    stop_rho = np.sqrt (stop_x**2 + stop_y**2)
    start_z  = frame[L6_key]['z_prime'] - 350.
    stop_z   = frame[L6_key]['stop_z_prime'] - 350.
    events['start_rho'].append (start_rho)
    events['stop_rho'].append  (stop_rho)
    events['start_z'].append   (start_z)
    events['stop_z'].append    (stop_z)
    
    ### get true info
    ##  no true info if using data
    if 'data' in dtype or 'muon' in dtype:
        return events

    ##  neutrino true info from I3Tree
    primary = frame['I3MCTree'][0]
    logen   = np.log10 (primary.energy)
    coszen  = np.cos (primary.dir.zenith)
    events['logen'].append  (logen)
    events['coszen'].append (coszen)

    return events

####################################################
#### function to define pass ? and level6_bool
####################################################
def passed (frame, icc=False):

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

    ### must have reco frame and L6 frame
    if not frame.Has (reco_key): return False
    if not frame.Has (L6_key)  : return False

    ### apply corridor cut
    corridor = frame[L6_key]['corridor_doms_over_threshold']
    if icc:
        ## inverted corridor: keep > 1, remove <= 1
        if corridor <= 1: return False
    else:
        ## normal cut: remove > 1
        if corridor > 1: return False
    return True

####################################################
#### function to collect events
####################################################
def collect_events (dtype, icc=False):

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
            #try:
            frame = i3file.pop_physics ()
            # append when physics frame
            if frame:
                if passed (frame, icc=icc):
                    events = append (frame, dtype,
                                     events, nfiles,
                                     probmap=probmap,
                                     icc=icc)
            else:
                print ('| this frame is empty ...')
                continue
            #except:
            #    print ('| this file is busted :/')
            #    pass
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
    print ('| weighted   before L6cuts : {0:4f} mHZ'.format (np.sum (weight)*1000.))

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
    outfile = options.outdir + 'level6_postreco_'+ptype+'.pckl'
    icc = True if 'muon' in ptype else False

    ### loop through events
    events = collect_events (ptype, icc=icc)

    ### save it
    with open (outfile, 'wb') as f:
        cPickle.dump (events, f, protocol=2)
    f.close ()
