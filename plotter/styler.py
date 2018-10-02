#! /usr/bin/env python

####
#### By Elim Thompson (09/18/2018)
####
#### This script defines the default stylings
#### for all nutau paper plots.
####
#### For each plotting scripts, do
####  `from plotter.styles import *`
####
#### NOTE: This is not the best practice way.
####       But it provides consistency through
####       all nutau paper plots.
####
####################################################

#################################################
### import default variables from defaults
### These variables include colors/hatches
### for each data types, x/y labels etc
#################################################
from defaults import labels, hatches, colors, linestyles
from defaults import dm21, seconds_per_year

#################################################
### matplotlib tools
#################################################
import matplotlib
### no interactive plots
matplotlib.use ('Agg')
### import coloar map
from matplotlib import cm
### import plot-tools
import matplotlib.pyplot as plt
### import gridspec for multi plots
import matplotlib.gridspec as gridspec
### import AnchoredText for text in plots
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

### set up text styles for published plots
plt.rc ('text', usetex=True)
plt.rc ('font', family='sans-serif')
plt.rc ('font', serif='Computer Modern Roman')
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

### set up global axis ticks
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['xtick.major.pad']='5'
plt.rcParams['ytick.major.pad']='5'

### set up hatches color / linewidth
plt.rcParams['hatch.color'] = 'lightgray'
plt.rcParams['hatch.linewidth'] = 1.0
