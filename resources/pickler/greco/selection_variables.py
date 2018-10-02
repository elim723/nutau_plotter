#!/usr/bin/env python

####
#### By Elim Thompson (09/24/2018)
####
#### This script includes functions to pull variables
#### used during the GRECO selections.
####
#######################################################

from __future__ import print_function

#### icecube projects
from icecube import dataclasses, fill_ratio, genie_icetray

### GRECO MC variables for MC info
##  1. ftruth  = dataclasses.get_most_energetic_neutrino (frame['I3MCTree'])
##  2. fweight =  
grecoMC = {'nu': {'energy' :'truth/energy'      ,
                  'zenith' :'truth/dir.zenith'  ,
                  'azimuth':'truth/dir.azimuth' ,
                  'ptype'  :'truth/pdg_encoding',
                  ''}}

############################################################
#### GRECO L4
#### L4 BDT varaibles at level 5 I3 files:
####    -- frame['L4BDT'] (dict)
############################################################
grecoL4 = {'c2qr6'         :'L4BDT/BDTC2QR6'        ,
           'ilf_speed'     :'L4BDT/BDTiLinefitSpeed',
           'vertex_guess_z':'L4BDT/BDTVGZ'          ,
           'bdt_score'     :'L4BDT/BDTScore'        ,
           'qr6'           :'L4BDT/BDTQR6'          ,
           'nabove200'     :'L4BDT/BDTNAbove200'    ,
           'toi_ratio'     :'L4BDT/BDTToIEvalRatio' }

def pull_greco_l4 (frame):
    
    ### holder for greco L4 variables
    L4vars = {}
    for key, path in grecoL4.items ():
        ## define frame obj and variable name
        frame_obj, var_name = path.split ('/')
        ## access via simple dictionary
        L4vars[key] = frame[frame_obj][var_name]
    ### return L4 var dictionary
    return L4vars

############################################################
#### GRECO L5
#### L5 BDT varaibles at level 5 I3 files:
####    -- frame['L5BDT'] (dict)
############################################################
grecoL5 = {'firstHLC_rho'   :'L5BDT/BDTFirstHLCRho'    ,
           'acc_time'       :'L5BDT/BDTAccumulatedTime',
           'vich'           :'L5BDT/BDTVICH'           ,
           'cog_seperation' :'L5BDT/BDTSeparation'     ,
           'ztravel'        :'L5BDT/BDTZTravel'        ,
           'spe11_zenith'   :'L5BDT/BDTZenith'         ,
           'bdt_score'      :'L5BDT/BDTScore'          }

def pull_greco_l5 (frame):

    ### holder for greco L5 variables
    L5vars = {}
    for key, path in grecoL5.items ():
        ## define frame obj and variable name
        frame_obj, var_name = path.split ('/')
        ## access via simple dictionary
        L5vars[key] = frame[frame_obj][var_name]
    ### return L5 var dictionary
    return L5vars

############################################################
#### GRECO L6
#### L6 varaibles at level 6 I3 files
####    -- frame['TauL6_variables']      (dict)
####    -- frame['taul6_bool']           (icetray.I3Bool)
####    -- frame['TauL6_FillRatio']      (fill_ratio.I3FillRatioInfo)
####    -- frame['SRTTWOfflinePulsesDC'] (dataclasses.I3RecoPulseSeriesMapMask)
############################################################
grecoL6 = {'finitereco_rho' :'TauL6_variables/FiniteRecoRho'       ,
           'finitereco_z'   :'TauL6_variables/FiniteRecoZ'         ,
           'corridor_nch'   :'TauL6_variables/CorridorCount'       ,
           ## fill_ratio project is needed
           'fill_ratio'     :'TauL6_FillRatio/fill_ratio_from_mean', 
           ## nchannel is calculated from pulseseries
           'SRT_nch'        :'SRTTWOfflinePulsesDC/'               ,
           'taul6_bool'     :'TauL6_bool/'                         }

def pull_greco_l6 (frame):

    ### holder for greco L6 variables
    L6vars = {}
    for key, path in grecoL6.items ():
        ## define frame obj and variable name
        frame_obj, var_name = path.split ('/')
        ## get the values accordingly
        if key == 'taul6_bool':
            value = frame[frame_obj].value
        elif key == 'fill_ratio':
            value = frame[frame_obj].fill_ratio_from_mean
        elif key == 'SRT_nch':
            hitmap = frame[frame_obj].apply (frame)
            value  = np.sum ([ len(hitmap[dom])>0
                               for dom in hitmap.keys() ])
        else:
            # access TauL6_variables dictionary
            value = frame[frame_obj][var_name]
        ## store the value
        L6vars[key] = value
    ## return L6 var dictionary
    return L6vars

############################################################
#### GRECO L7
#### L7 varaibles at level 7 I3 files 
#### 
#### L7 variables are slightly complicated. An array of
#### variables are defined instead of a dictionary.
############################################################
grecoL7 = ['energy', 'zenith', 'tracklength', 'deltaLLH',
           'z', 'rho', 'gev_per_channel', 't_rms',
           'charge_rms_normalized', 'errata_length']

def get_hit_info (hitmap):
    ### initialize holders
    nch, charge, charge2 = 0., 0., 0.
    times = []
    ### loop through each hit DOM
    for dom in hitmap.keys ():
        ## add number of hit DOMs
        nch += len (hitmap[dom])>0
        ## add charge
        thisq = sum ([ hit.charge for hit in hitmap[dom] ])
        charge += thisq
        charge2 += thisq**2
        ## append to first photon arrival time
        times.append (hitmap[dom][0].time)
    ### return hit info
    return nch, charge, charge2, np.array (times)

def pull_greco_l7 (frame):

    ### initialize info holder
    L7vars = {}

    ###  1. interacting particle 
    reco_cascade = frame['Pegleg_Fit_MN_tol10HDCasc']
    reco_track   = frame['Pegleg_Fit_MN_tol10Track']    
    energy = reco_cascade.energy + reco_track.energy
    zenith = reco_track.dir.zenith
    tracklength = reco_track.length

    L7vars['energy'] = energy
    L7vars['zenith'] = zenith
    L7vars['tracklength'] = tracklength
    
    ### 2. pegleg fit (used by Martin)
    pegleg_llh  = frame['Pegleg_Fit_MN_tol10FitParams'].logl 
    monopod_llh = frame['Monopod_bestFitParams'].logl
    deltaLLH = pegleg_llh - monopod_llh
    
    L7vars['deltaLLH'] = deltaLLH

    ### 3. vertex
    x = reco_cascade.pos.x
    y = reco_cascade.pos.y
    z = reco_cascade.pos.z
    rho = np.sqrt ((x-46.29)**2 + (y+34.88)**2)

    L7vars['z'] = z
    L7vars['rho'] = rho

    ### 4. data / mc agreements
    hitmap = frame['SRTTWOfflinePulsesDC'].apply (frame)
    nch, charge, charge2, times = get_hit_info (hitmap)
    
    t_rms = np.std (times)
    gev_per_channel = energy / nch
    charge_rms_normalized = np.sqrt (charge2)/charge
    errata_length = len(frame["CalibrationErrata"].keys())
    
           'gev_per_channel', 't_rms',
           'charge_rms_normalized', 'errata_length']    
