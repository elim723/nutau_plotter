#!/bin/bash

####
#### By Elim Thompson (09/28/2018)
#### 
#### This bash script runs all plotting scripts.
####
#### command:
#### $ bash run_all.sh < out directory >
####
####################################################

### read input argument
outdir=$1
### define where data files are
pdata='/data/user/elims/nutau_pickled_data/'

############################
### make directories
############################
mkdir $outdir/dragon1D/
mkdir $outdir/dragon2D/
mkdir $outdir/greco1D/
mkdir $outdir/greco2D/
mkdir $outdir/syseffects/
mkdir $outdir/classification/
mkdir $outdir/muon_histograms/

############################
### DRAGON
############################
### 1D distribution 
python plot_dragon1D.py --outdir $outdir/dragon1D/ --fitnorm
### 2D distribution 
python plot_dragon2D.py --outdir $outdir/dragon2D/

############################
### GRECO
############################
### 1D distribution plots
python plot_greco1D.py --outdir $outdir/greco1D/ --fitnorm
### 2D distribution plots
python plot_greco2D.py --outdir $outdir/greco2D/
### numu contour
python plot_greco_numu.py --outdir $outdir/
### signal / sqrt (background)
python plot_sig_to_bg.py --outdir $outdir/

############################
### PID
############################
python plot_classification.py --outdir $outdir/classification/ --infile $pdata/pid_track_fractions.p

############################
### PegLeg Resolutions
############################
python plot_resolutions.py --outdir $outdir/resolution/ --infile $pdata/pegleg_resolutions.p 

############################
### Muon Histograms
############################
python plot_muon_histograms.py --outdir $outdir/muon_histograms/
python plot_muon_histograms.py --greco  --outdir $outdir/muon_histograms/
python plot_muon_histograms.py --dragon --outdir $outdir/muon_histograms/

############################
### Systematic Effects
############################
python plot_sys_effects.py --outdir $outdir/syseffects/ --infile $pdata/greco/syseffects_events.p 
