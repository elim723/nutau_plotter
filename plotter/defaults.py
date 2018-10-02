#!/usr/bin/env python

####
#### By Elim Cheung (07/24/2018)
####
#### This file defines the global variables for
#### various plots in nutau paper:
####
####   1. Common variables
####
####   2. Colors
####     -- data types
####     -- color maps
####
####   3. Hatches for data types
####
####   4. Linestyles for numu 90% contours
####
####   5. Labels
####     -- systematic
####     -- data types (members)
####     -- observables/event variables
####     -- analyses for GRECO and DRAGON
####
###############################################################

####################################
#### Common Variables
####################################
dm21 = 7.49E-5
seconds_per_year = 3600.*24.*365.24
## 2.27 effective years for GRECO
greco_scalar = 2.27 * seconds_per_year
## 2.5  effective years for GRECO
dragon_scalar = 2.5 * seconds_per_year

####################################
#### Colors
####################################
colors = {'numucc' : 'steelblue'   ,
          'nuecc'  : 'seagreen'    ,
          'nutaucc': 'orangered'   ,
          'numunc' : 'navy'        ,
          'nuenc'  : 'darkgreen'   ,
          'nutaunc': 'khaki'       ,
          'numu'   : 'steelblue'   ,
          'nue'    : 'seagreen'    ,
          'nutau'  : 'orangered'   ,
          'nunc'   : 'darkmagenta' ,
          'muon'   : 'lightpink'   ,
          'muongun': 'lightpink'   ,
          'corsika': 'lightpink'   ,
          'iccdata': 'lightpink'   ,
          'noise'  : '#6C4000'     ,
          'data'   : 'black'       ,
          'mc'     : 'red'         ,
          'mcerr'  : 'red'         ,
          'bestfit': 'red'         ,
          'noosc'  : 'blue'        ,

          ## cmaps
          'counts'   : 'Blues'     ,
          'variances': 'Reds'      ,
          'effects'  : 'Spectral'  ,

          ## contours from other exp
          't2k17'    : 'magenta'   ,
          'minos16'  : 'green'     ,
          'sk17'     : 'red'       ,
          'nova17'   : 'blueviolet',
          'ic17'     : 'black'     ,
          'greco'    : 'blue'      }

####################################
#### Hatches for data types only
####################################
hatches = {'numucc' : "\\\\",
           'nuecc'  : "//"  ,
           'nutaucc': ''    ,
           'numunc' : '\\'  ,
           'nuenc'  : '/'   ,
           'nutaunc': '|'   ,
           'numu'   : "\\\\",
           'nue'    : "//"  , 
           'nutau'  : ''    ,
           'nunc'   : '\\'  ,
           'muongun': '+'   ,
           'corsika': '+'   ,
           'muon'   : '+'   ,
           'iccdata': '+'   ,
           'noise'  : '-'   }

####################################
#### Linetyls for numu 90% contours
####################################
linestyles = {'t2k17'  :'--',
              'minos16':'-.',
              'nova17' :'-.',
              'sk17'   :':' ,
              'ic17'   :'-' ,
              'greco'  :'-'  }

####################################
#### Labels
####################################
labels = { 'dm31'            : r'$\Delta$m$^2_{31}$ ($10^{-3}$ eV$^2$)',
           'dm32'            : r'$\Delta$m$^2_{32}$ ($10^{-3}$ eV$^2$)',
           'sin2theta23'     : r'sin$^2$ $\theta_{23}$'        ,
           'theta23'         : r'$\theta_{23}$'                ,
           'theta13'         : r'$\theta_{13}$'                ,
           'muon_flux'       : r'$\gamma_{\mu}$'               ,
           'gamma'           : r'$\gamma_{\nu}$'               ,
           'nue_numu_ratio'  : r'N$_{\nu_e}$ / N$_{\nu_{\mu}}$',
           'barr_uphor_ratio': r'up / horizontal'              ,
           'barr_nubar_ratio': r'$\nu$ / $\bar{\nu}$ '         ,
           'axm_res'         : r'axial res'                    ,
           'axm_qe'          : r'axial qe'                     ,
           'DISa_nu'         : r'DIS $\nu$'                    ,
           'DISa_nubar'      : r'DIS $\bar{\nu}$'              ,
           'spe_corr'        : r'SPE correction'               ,
           'nyears'          : r'global normalization'         ,
           'norm_atmmu'      : r'N$_{\mu}$'                    ,
           'norm_noise'      : r'N$_{\text{noise}}$'           ,
           'norm_numu'       : r'N$_{\nu_{\mu}}$'              ,
           'norm_nutau'      : r'N$_{\nu_{\tau}}$'             ,
           'norm_nc'         : r'N$_{\text{NC}}$'              ,
           'norm_nugen'      : r'N$_{\text{nugen}}$'           ,
           'norm_nugenHE'    : r'N$_{\text{nugenHE}}$'         ,
           'norm_corsika'    : r'N$_{\text{corsika}}$'         ,
           'domeff'          : r'DOM efficiency'               ,
           'holeice'         : r'hole ice'                     ,
           'forward'         : r'hole ice forward'             ,
           'coin'            : r'coincident fraction'          ,
           'absorption'      : r'absorption scaling'           ,
           'scattering'      : r'scattering scaling'           ,

           ## data type members
           'numucc' : r'$\nu_{\mu} + \bar{\nu}_{\mu}$ CC'  ,
           'nuecc'  : r'$\nu_e + \bar{\nu}_e$ CC'          ,
           'nutaucc': r'$\nu_{\tau} + \bar{\nu}_{\tau}$ CC',
           'numunc' : r'$\nu_{\mu} + \bar{\nu}_{\mu}$ NC'  ,
           'nuenc'  : r'$\nu_e + \bar{\nu}_e$ NC'          ,
           'nutaunc': r'$\nu_{\tau} + \bar{\nu}_{\tau}$ NC',
           'numu'   : r'$\nu_{\mu}$ NC $+$ CC'             ,
           'nue'    : r'$\nu_e$ NC $+$ CC'                 ,
           'nutau'  : r'$\nu_{\tau}$ NC $+$ CC'            ,
           'nunc'   : r'$\nu + \bar{\nu}$ NC'              ,
           'muongun': r'atm $\mu$'                         ,
           'corsika': r'atm $\mu$'                         ,
           'muon'   : r'atm $\mu$'                         ,
           'iccdata': r'inverted data'                     ,
           'noise'  : r'noise'                             ,
           'data'   : r'data'                              ,
           'mc'     : r'total MC'                          ,
           'mcerr'  : r'MC Uncertainty'                    ,
           'bestfit': r'Best fit'                          ,
           'noosc'  : r'Null hypothesis'                   ,

           ## analysis observables
           'reco_e'     : r'Reconstructed Energy (GeV)'           ,
           'reco_loge'  : r'log$_{10}$ Reconstructed Energy (GeV)',
           'reco_z'     : r'Reconstructed Zenith (radian)'        ,
           'reco_zdeg'  : r'Reconstructed Zenith (degree)'        ,
           'reco_cz'    : r'Reconstructed Cos(Zenith)'            ,
           'reco_L/E'   : r'Reconstructed L/E (km / GeV)'         ,
           'reco_length': r'Reconstructed track length (m)'       ,
           'reco_dllh'  : r'$\Delta$ log $L_{\text{reco}}$'       ,
           'mc_e'       : r'MC Truth Energy (GeV)'                ,
           'mc_loge'    : r'Log$_{10}$ MC Truth Energy (GeV)'     ,
           'mc_z'       : r'MC Truth Zenith (radian)'             ,
           'mc_zdeg'    : r'MC Truth Zenith (degree)'             ,
           'mc_cz'      : r'MC Truth Cos(Zenith)'                 ,

           ## lower level selection variables
           'nabove200'            : r'Charges (PE) above z = 200m'                              ,
           'c2qr6'                : r'Charge ratio of first 600ns to total without first 2 hits',
           'qr6'                  : r'Charge ratio of first 600ns to total'                     ,
           'vertexguessZ'         : r'Vertical position (m) of the first hit in time'           ,
           'toi_eigenvalueratio'  : r'Tensor of inertia eigenvalue ratio'                       ,
           'ilinefit_Speed'       : r'Speed (c) from improved LineFit'                          ,
           'timeTo75'             : r'Time (ns) to accumulate 75$\%$ of charges'                ,
           'vertexguessrho'       : r'Radial position $\rho_{DOM}$ (m) of the earliest HLC DOM' ,
           'cog_separation'       : r'Distance (m) between CoGs of 1st and 4th quartiles'       ,
           'separation'           : r'Distance (m) between CoGs of 1st and 4th quartiles'       ,
           'ztravel'              : r'Z travel (m)'                                             ,
           'spe_cz'               : r'SPE 11 cos zenith'                                        ,
           'fillratio'            : r'Fill ratio'                                               ,
           'finiterecorho'        : r'Vertex radial position $\rho$ (m) from FiniteReco'        ,
           'finiterecoz'          : r'Vertex Z position (m) from FiniteReco'                    ,
           'corridornch'          : r'Number of channels along corridors'                       ,
           'nchannel'             : r'Number of cleaned hit DOMs'                               ,
           'charge_rms_normalized': r'Normalized RMS of total charges (PE)'                     ,
           'vetocausalhits'       : r'Causually connected veto hits (PE)'                       ,
           't_rms'                : r'RMS event time (ns)'                                      ,
           'dcfiducialpe'         : r'Chargs (PE) in DeepCore fiducial volume'                  ,
           'microcountpe'         : r'MicroCountPE (PE)'                                        ,
           'microcounthits'       : r'MicroCountHits'                                           , 
           'qtotal'               : r'Total charges (PE)'                                       ,
           'total_charge'         : r'Total charges (PE)'                                       ,
           'rtvetoseries250pe'    : r'RTVetoSeries250PE (PE)'                                   ,
           'dcfilterpulses_vetope': r'DCFilterPulses vetoPE (PE)'                               ,
           'geV_per_channel'      : r'Reconstructed energy per channel (GeV)'                   ,
           'bdtscore'             : r'BDT score'                                                ,
           'bdt_score'            : r'BDT score'                                                ,
           'reco_x'               : r'Reconstructed vertex X position (m)'                      ,
           'reco_y'               : r'Reconstructed vertex Y position (m)'                      ,
           'reco_z'               : r'Reconstructed vertex Z position (m)'                      ,
           'reco_rho'             : r'Reconstructed vertex radial $\rho$ position (m)'          ,
           'start_z'              : r'Reconstructed track starting Z position (m)'              ,
           'start_rho'            : r'Reconstructed track starting radial $\rho$ position (m)'  ,
           'stop_z'               : r'Reconstructed track stopping Z position (m)'              ,
           'stop_rho'             : r'Reconstructed track stopping radial $\rho$ position (m)'  ,

           ## contours from other exp
           't2k17'  :r'T2K 2017'     ,
           'minos16':r'MINOS 2016'   ,
           'sk17'   :r'SK 2017'      ,
           'nova17' :r'NO$\nu$A 2017',
           'ic17'   :r'IC 2017'      ,
           
           ## others
           'dragon' :r'Analysis $\mathcal{B}$'     ,
           'greco'  :r'Analysis $\mathcal{A}$'     ,
           'dchi2'  :r'$\Delta\chi^2$'             ,
           'chi2'   :r'$\chi^2$'                   ,
           'modchi2':r'modified $\chi^2$'          ,
           'poisson':r'2 $\times$ Poisson LLH'     ,
           'barlow' :r'2 $\times$ Barlow LLH'      ,
           'trials' :r'number of trials'           ,
           'counts' :r'number of counts in 3 years' }

#######################################
#### Change in systematic effects
#### on analysis histograms
#### 
#### By default, all are 1 sigma up,
#### or 10%.
####
#### For oscillation parameters,
#### defaults are 1 sigma up from
#### nufit 2016:
#### http://www.nu-fit.org/?q=node/139
#######################################
deltas = { 'dm31'            : 2.56e-3, 
           'theta23'         : 0.7523 ,
           'theta13'         : 0.1503 ,
           'muon_flux'       : 1.0    ,
           'gamma'           : 0.1    ,
           'nue_numu_ratio'  : 1.05   ,
           'barr_uphor_ratio': 1.0    ,
           'barr_nubar_ratio': 1.0    ,
           'axm_res'         : 1.0    ,
           'axm_qe'          : 1.0    ,
           'DISa_nu'         : 0.0757 ,
           'DISa_nubar'      : 0.1008 ,
           'spe_corr'        : 0.04   ,
           'nyears'          : 3.0    ,
           'norm_atmmu'      : 1.1    ,
           'norm_noise'      : 1.1    ,
           'norm_numu'       : 1.1    ,
           'norm_nutau'      : 1.1    ,
           'norm_nc'         : 1.2    ,
           'norm_nugen'      : 1.1    ,
           'norm_nugenHE'    : 1.1    ,
           'norm_corsika'    : 1.1    ,
           'domeff'          : 1.1    ,
           'holeice'         : 35     ,
           'forward'         : 1      ,
           'coin'            : 0.1    ,
           'absorption'      : 1.1    ,
           'scattering'      : 1.1     }
