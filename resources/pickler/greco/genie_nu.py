#!/usr/bin/env python
from __future__ import print_function

###################### PRINTOUTS ########################
## python  get_genie_pckl_from_michael.py --infile --outpath --datatype --datset (--merged)
## 

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

def get_cut (nch, recoE_nchannel, trms, isnugen=None):
    cut =  recoE_nchannel < 3.0
    cut *= trms < 800.
    cut *= abs ((recoE_nchannel-1)/(trms-800.)) < abs ((3.-1)/(500-800.))
    cut *= nch <14
    boolean = np.logical_or (nch>=14, cut)
    if not merged and not isnugen==None:
        boolean = np.logical_and (boolean, ~isnugen)
    return boolean

##################################################################################
##### Parsing variables
##################################################################################
usage = "%prog [options] --infile <i3file> --outpath <path of out pickled file>"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--infile", type="string", default = '~/',
                  help = "oscfit pickled files")
parser.add_option("-o", "--outpath", type = "string", default = '/data/i3home/elims/oscillation/fitter/wisc_work/',
                  help = "output path")
parser.add_option("-d", "--datatype", type="string", default = 'numucc',
                  help = "datatype")
parser.add_option("--dataset", type="string", default='585',
                  help = "specify dataset: 551 552 553 554 555 556 560 561 562 563 564 565 567 568 571 585")
parser.add_option("-m", "--merged", action="store_true", default = False,
                  help = "datatype")
(options, args) = parser.parse_args()
infile   = options.infile
merged   = options.merged
outpath  = options.outpath
dataset  = options.dataset
datatype = options.datatype

ptype = '12' if 'nue' in datatype else '14' if 'numu' in datatype else '16'
interaction = 1 if 'cc' in datatype else 2

is_dima_set  = True if 'dima' in dataset or dataset in ['600', '601', '603', '604', '605', '606', '608',
                                                        '610', '611', '612', '613', '620', '621', '622', '623', '624'] else False
is_noRDE_set = True if dataset in ['640', '641', '643', '644', '645', '646', '648', '660', '661', '662', '663', '670', '671', '672', '673', '674'] else False
is_nugen = True if dataset in ['11981', '11029', '11477'] else False
is_coin = True if 'coin' in dataset else False
is_SPiceHD = True if dataset in ['544', '545', '546', '547', '548', '549'] else False
is_SPice3X = True if dataset in ['900', '950'] else False

domeff_sets = { '551':'0.85', '552':'0.90', '553':'0.95', '554':'1.05', '555':'1.10', '556':'1.15',
                '601':'0.88', '603':'0.94', '604':'0.97', '605':'1.03', '606':'1.06', '608':'1.12',
                '641':'0.88', '643':'0.94', '644':'0.97', '645':'1.03', '646':'1.06', '648':'1.12',
                'domEff_0.79':'0.79', 'domEff_0.89':'0.89', 'domEff_1':'1', 'domEff_1.09':'1.09', 'domEff_1.19':'1.19',
                'domEff_0.69_dima':'0.69', 'domEff_0.79_dima':'0.79', 'domEff_1_dima':'1'}

for setnum in ['550', '560', '561', '562', '563', '564', '565', '566', '567', '568', '571', '572', '573', '585', '585-550',
               '600', '610', '611', '612', '613', '620', '621', '622', '623', '624', '11981', '11029', '11477',
               '640', '660', '661', '662', '663', '670', '671', '672', '673', '674', 'coincident',
               '544', '545', '546', '547', '548', '549', '900', '950', 
               'holeIce_30', 'holeIce_50', 'holeIce_100', 'holeIce_25_dima', 'holeIce_15_dima', 'holeIce_30_dima', 'forward_m2_dima', 'forward_m4_dima']:
    domeff_sets[setnum] = '1'

holeice_sets = { '560':'100', '561':'30', '562':'0' , '563':'40', '564':'45', '565':'55', '566':'67' ,
                 '571':'33' , '572':'36', '573':'80', #### hadronic sets: '567':'50', '568':'50',                                                                 
                 '610':'15' , '611':'20', '612':'30', '613':'35',
                 '660':'15' , '661':'20', '662':'30', '663':'35',
                 '544':'170', '545':'14', '546':'85', '547':'32', '548':'70', '549':'10',
                 '900':'30' , '950':'30', 
                 'holeIce_30':'30', 'holeIce_50':'50', 'holeIce_100':'100',
                 'holeIce_30_dima':'30', 'holeIce_25_dima':'25', 'holeIce_15_dima':'15'}

forward_sets = {'620':'p2', '621':'n5', '622':'n3', '623':'p1', '624':'n1',
                '670':'p2', '671':'n5', '672':'n3', '673':'p1', '674':'n1',
                '900':'p1', '950':'n1', 
                'forward_m2_dima':'n2', 'forward_m4_dima':'n4'}

for setnum in ['600', '601', '603', '604', '605', '606', '608', '550', '551', '552', '553', '554', '555', '556', '585', '585-550', '11981', '11029', '11477',
               '640', '641', '643', '644', '645', '646', '648', 'coincident',
               'domEff_0.79', 'domEff_0.89', 'domEff_1', 'domEff_1.09', 'domEff_1.19', 'domEff_0.69_dima', 'domEff_0.79_dima', 'domEff_1_dima']:
    holeice_sets[setnum] = '30' if ('dima' in setnum and 'domEff' in setnum) else '25' if is_dima_set or is_noRDE_set or is_coin else '50'
    forward_sets[setnum] = '0'
for setnum in ['610', '611', '612', '613', '660', '661', '662', '663', '544', '545', '546', '547', '548', '549',
               'holeIce_30_dima', 'holeIce_25_dima', 'holeIce_15_dima']:
    forward_sets[setnum] = '0'
for setnum in ['620', '621', '622', '623', '624', '670', '671', '672', '673', '674', 'forward_m2_dima', 'forward_m4_dima']:
    holeice_sets[setnum] = '25' if is_noRDE_set else '30'

domeff  = domeff_sets[dataset]
holeice = holeice_sets[dataset]
forward = forward_sets[dataset] if dataset in forward_sets.keys () else None
ext = '_forward'+forward_sets[dataset] if is_dima_set else \
      '_forward'+forward_sets[dataset]+'_SPice3X' if is_SPice3X else \
      '_forward'+forward_sets[dataset]+'_SPiceHD' if is_SPiceHD else \
      '_forward'+forward_sets[dataset]+'_noRDE' if is_noRDE_set else \
      '_nugen' if is_nugen else '_coincident' if is_coin else \
      '_'+dataset if dataset in ['585', '550'] else ''
ext = ext + '_qIndPegleg'
if merged: ext += '_merged'
outname = datatype + '_domeff' + domeff + '_holeice' + holeice + ext + '.p'
print ('{0} {1} domeff, holeice, forward: {2} {3} {4}'.format (datatype, dataset, domeff, holeice, forward))
outfile = outpath + outname

##################################################################################
##### start repickling
##################################################################################

#### set intereaction type kept for genie sets
events = Map ({})

#### open pickled file
with open (infile, 'rb') as f:
    original = cPickle.load(f)
f.close()
isnugen = original['isNuGen'] if original.has_key ('isNuGen') else None
if not isnugen==None:
    print ('number of isnugen: {0}'.format (np.sum (isnugen.astype (int))))
cut = np.logical_and (original['interaction']==interaction, get_cut (original['nchannel'], original['GeV_per_channel'], original['t_rms'], isnugen=isnugen))
print ('nevents before, after: {0}, {1}'.format (len (original['reco_energy']), len (original['reco_energy'][cut])))

#### for all datatype: reco / hits / w (empty for neutrinos)
events['w'] = []

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

events['mc']   = Map({ 'e'           :original['energy'][cut]        , 
                       'z'           :original['zenith'][cut]        , 
                       'cz'          :np.cos(original['zenith'][cut]), 
                       'pdg'         :original['ptype'][cut]         ,
                       'X'           :original['x'][cut]             ,
                       'Y'           :original['y'][cut]             ,
                       'Z'           :original['z'][cut]             ,
                       'a'           :original['azimuth'][cut]       ,
                       'nprimaries'  :original['nprimaries'][cut]    })

events['misc'] = Map({ 'detector_w':np.ones (len (original['reco_energy'][cut])),
                       'fluxes'    :np.array ([original['weight_e'][cut], original['weight_mu'][cut]]).T,
                       'genie_x'   :original['GENIE_x'][cut],
                       'genie_y'   :original['GENIE_y'][cut],
                       'genie_w'   :[],
                       'scattering':original['scattering'][cut],
                       'ma_qe'     :original['ma_qe'][cut],
                       'ma_res'    :original['ma_res'][cut],
                       'isnugen'   :original['isNuGen'][cut]         })

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

print ('## Save {0} events to {1}'.format(len(events.reco.e), outfile))
##### writing out pickled file
with open(outfile, "wb") as f:
    cPickle.dump( events, f, protocol=2 )
f.close ()

print ('## Done.')
