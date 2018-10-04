#!/usr/bin/env python

####
#### By Elim Thompson (09/28/2018)
####
#### This script made GRECO signal
#### to sqrt (background) per analysis
#### bin.
####
#### Command to run
####  $ python plot_sig_to_bg.py --outdir <plot>
####
###############################################


#### import standard packages
from __future__ import print_function
from optparse import OptionParser
from glob import glob
import numpy as np
import cPickle, os, sys

from plotter.defaults import labels
from plotter.histogram3D import Histogram3D
from plotter.printer import sigtobg_printer_header
from plotter.printer import sigtobg_printer_events
from plotter.printer import sigtobg_printer_ender

#### import default analysis histogram edges
from analyzer.misc import default_edges as edges
from analyzer.misc import greco_nyears, seconds_per_year

scalar = greco_nyears * seconds_per_year

#########################################################
#### define global variables
#########################################################
### directory where preweighted pickled
### data are located
sample = 'greco'
ppath  = '/data/user/elims/nutau_pickled_data/greco/preweighted/'

variables = ['reco_logen', 'reco_coszen', 'tracklength']

#########################################################
#### define plot settings
####
#### NOTE: To make plots of other variabls, add a
####       dictionary to the settings dictionary
####       like below. And make sure you have the
####       same default key names as the ones in
####       resources/plotter/defaults.py.
####
####      cmapname = color map name (google cmap names)
####      htype    = histogram type:
####                 1. numucc = distribution of numu CC
####                             (same for any member)
####                 2. nurate = total neutrino rates
####                 3. mufrac = muon fraction
####                 4. s/b    = nutau / sqrt (not nutau)
####       NOTE: to implement new htype, go to
####             `plotter/histogram2D/_get_content ()`
#########################################################
### plot settings
settings = {'cmapname'    :'Reds'                   ,
            'htype'       :'s/b'                    ,
            'stretch_cscd':False                    ,
            'xedges'      :np.log10 (edges['e'])    ,
            'yedges'      :np.cos (edges['z'])[::-1],
            'zedges'      :edges['p']               ,
            'xticks'      :np.log10 (edges['e'])    ,
            'yticks'      :[-1.0, -0.6, -0.2, 0.2, 0.6, 1.0],
            'xticklabels' :[5.6, 7.5, 10, 13, 18, 24, 32, 53, 56],
            'yticklabels' :[-1.0, -0.6, -0.2, 0.2, 0.6, 1.0] }  

#########################################################
#### set plot format
#########################################################
### global format
format_global     = {'tick_fontsize' :13 ,
                     'label_fontsize':15 ,
                     'title_fontsize':17 ,
                     'grid_alpha'    :0.2,
                     'grid_linestyle':'-',
                     'grid_linewidth':0.5 }

### 2D histograms
format_histograms = {'hist_alpha'         :0.7  ,
                     'hist_print'         :False,
                     'hist_print_rounding':1    ,
                     'hist_print_fontsize':3.0   }

### color bar
format_colorbar   = {'colorbar_alpha'      :0.7            ,
                     'colorbar_format'     :'%10.1f'       ,
                     'colorbar_tick_length':2              ,
                     'colorbar_label'      :r'S / $\sqrt{\text{B}}$',
                     'colorbar_ranges'     :(0.0, 4.0)     ,
                     'colorbar_nbins'      :5               }


### Put all formats into a list
plot_formats = [settings, format_global, 
                format_histograms,
                format_colorbar  ]

#########################################################
### Functions to plot via `Histogram3D`
#########################################################
def plot_sig_to_bg (outdir, events, sample, stretch_cscd=False):

    ''' plot all variables in events

        :param outdir (str) : user's output directory
        :param events (dict): all event variables needed
    '''

    ### define a histogram instance
    ### with all global plot formats 
    histo = Histogram3D ('greco', stretch=stretch_cscd)

    ### set observable keys
    histo.xobs = 'reco_logen'
    histo.yobs = 'reco_coszen'
    histo.zobs = 'tracklength'

    ### set plot styling / formatting
    for formats in plot_formats:
        histo.set_style (**formats)

    ## plot and save !
    histo.plot (events)
    outname = sample+'_signal_to_background.pdf'
    histo.save (outdir+outname)
    return

#########################################################
#### functions to collect events
#########################################################
def get_events (filename):

    ''' get events

        :param filename  (str): path to the pickled file

        :return event (dict): variables / weights needed
    '''

    ### open file
    with open (filename, 'rb') as f:
        ddict = cPickle.load (f)
    f.close ()

    ### apply cuts
    cutkey = 'level7_bins' if 'data' in filename else \
             'level7_bin'  ### for mc 
    cut = ddict[cutkey]

    ### collect events as needed
    ## weight in Hz!
    weight = ddict['weight'].flatten () *scalar
    event = {'weight':weight[cut]}
    for key in ddict.keys ():
        if key.lower () in variables: 
            event[key.lower ()] = ddict[key][cut]
    return event

def collect_events ():

    ''' collect event variables as needed

        :return events (dict): all event variables needed
    '''

    ### holder for all events needed
    events = {}

    ### define filenames
    filenames = sorted (glob (ppath + 'level7*.pckl'))

    ### loop through each filename (i.e. each member)
    for filename in filenames:
        ## define holder for this member
        member = os.path.split (filename)[1].split ('.')[0].split ('_')[-1]
        ## skip data
        if member == 'data': continue
        ## get events
        events[member] = get_events (filename)

    ### print progress
    sigtobg_printer_events (events)
    return events

######################################
#### Everything starts here :)
######################################
if __name__ == "__main__":

    #### parse options/params
    usage = "%prog [--outdir <directory for plots>]"
    parser = OptionParser(usage=usage)
    parser.add_option ("--outdir", type="string", default='~/',
                       help = "out directory for plots")
    (options, args) = parser.parse_args()
    outdir  = options.outdir

    #### print header
    sigtobg_printer_header (outdir, sample.upper ())

    #### get events
    events = collect_events ()

    #### want to stretch cascade bin?
    stretch_cscd = False if sample=='dragon' else \
                   settings['stretch_cscd']

    #### plot variables
    ##   By default, this script plot the 
    ##   S/sqrt(B) from GRECO sample
    plot_sig_to_bg (outdir, events, sample,
                    stretch_cscd = stretch_cscd)

    #### print ender
    sigtobg_printer_ender (outdir)

