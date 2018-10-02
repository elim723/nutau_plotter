#!/usr/bin/env python
from __future__ import print_function

###################### PRINTOUTS ########################
## python  get_corsika_muongun_pckl_from_michael.py --infile --outpath

import socket
print ('########################################################################################')
print ('#### This job is running on {0}'.format(socket.gethostname()))
print ('########################################################################################')
print (' ')
#########################################################

from optparse import OptionParser
import numpy as np
import cPickle, os
from misc_tools import Map
from GoodRunSampler import getRuns, getTime

def get_cut (nch, recoE_nchannel, trms):
    cut =  recoE_nchannel < 3.0
    cut *= trms < 800.
    cut *= abs ((recoE_nchannel-1)/(trms-800.)) < abs ((3.-1)/(500-800.))
    cut *= nch <14
    return np.logical_or (nch>=14, cut)

##################################################################################
##### Parsing variables
##################################################################################
usage = "%prog [options] --infile <i3file> --outpath <path of out pickled file>"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--infile", type="string", default = '~/',
                  help = "oscfit pickled files")
parser.add_option("-o", "--outpath", type = "string", default = '/data/i3home/elims/oscillation/fitter/wisc_work/',
                  help = "output path")
(options, args) = parser.parse_args()
infile   = options.infile
outpath  = options.outpath

##################################################################################
##### Set up basic variables
##### NOTE: muongun and noise weights are calculated here
#####       only neutrino detector weight is calculated here (not total weight)
##################################################################################
#### set up out file name
outfile = outpath + 'noise_qIndPegleg.p'

#### set intereaction type kept for genie sets
events = Map ({})

#### open pickled file
with open (infile, 'rb') as f:
    original = cPickle.load(f)
f.close()
cut = get_cut (original['nchannel'], original['GeV_per_channel'], original['t_rms'])
print ('originally have {0} events'.format (len (original['reco_energy'])))
print ('after noise cut... I have {0} events'.format (len (original['reco_energy'][cut])))

#### for all datatype: reco / hits / w (empty for neutrinos)
events['w'] = original['weight'][cut]

events['reco'] = Map({ 'e'           :original['reco_energy'][cut]        , 
                       'z'           :original['reco_zenith'][cut]        , 
                       'cz'          :np.cos(original['reco_zenith'][cut]), 
                       'pid'         :original['tracklength'][cut]        ,
                       'X'           :original['reco_x'][cut]             ,
                       'Y'           :original['reco_y'][cut]             ,
                       'Z'           :original['reco_z'][cut]             ,
                       'dllh'        :original['deltaLLH'][cut]           ,
                       'pegleg_rllh' :original['rLLH_Pegleg'][cut]        ,
                       'monopod_rllh':original['rLLH_Monopod'][cut]       ,
                       'a'           :original['reco_azimuth'][cut]       })

events['hits'] = Map({ 'qTotal' :original['charge'][cut],
                       'SRT_nCh':original['nchannel'][cut]           })

events['geo']  = Map({ 'dist_nearest_dom'   :original['dist_nearest_om'][cut]      ,
                       'dist_nearest_string':original['dist_nearest_string'][cut]  ,
                       'charge_asym'        :original['charge_rms_normalized'][cut],
                       'sorted_hits'        :[],
                       'nearest_dom'        :original['nearest_om'][cut]           ,
                       'nearest_s'          :original['nearest_string'][cut]       })


for key in ['L6', 'L5', 'L4', 'FR']:
    events[key] = Map ({})
    subkeys = ['firstHLC_rho', 'acc_time', 'vich','cog_seperation', 'ztravel', 'spe11_zenith', 'bdt_score'] if key=='L5' else \
              ['finitereco_pho', 'FR_timewindow', 'tauL6_icc_bool', 'fill_ratio_nch', 'SRT_nch', 'tauL6_bool', 'finitereco_z', 'corridor_nch', 'tauL6_ifr_bool'] if key=='L6' else \
              ['c2qr6', 'ilf_speed', 'vertex_guess_z', 'bdt_score', 'qr6', 'nabove200', 'toi_ratio'] if key=='L4' else \
              ['eventID', 'CalibrationErrata', 'month', 'runID', 'year', 'duration', 'day'] if key=='date' else \
              ['radius_from_mean', 'hit_counts', 'pulse_dt', 'vertex', 'mean_distance', 'ratio_from_mean'] ## FR
    for subkey in subkeys:
        events[key][subkey] = original['fr_from_nch'][cut] if subkey=='fill_ratio_nch' else \
                              original['nchannel'][cut]    if subkey=='SRT_nch'        else []

print ('## Save {0} (weighted: {1}) events to {2}'.format(len(events.reco.e), np.sum (events.w), outfile))
##### writing out pickled file
with open(outfile, "wb") as f:
    cPickle.dump( events, f, protocol=2 )
f.close ()

print ('## Done.')
