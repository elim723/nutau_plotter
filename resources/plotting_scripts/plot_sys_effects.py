#!/usr/bin/env python

####
#### By Elim Thompson (09/28/2018)
####
#### This script made GRECO signal
#### to sqrt (background) per analysis
#### bin.
####
#### Command to run
####  $ python plot_sys_effects.py --outdir <plot>
####
###############################################


#### import standard packages
from __future__ import print_function
from optparse import OptionParser
from glob import glob
import numpy as np
import cPickle, os, sys

from plotter.defaults import labels
from plotter.histogramND import HistogramND
from plotter.printer import syseffects_printer_header
from plotter.printer import syseffects_printer_events
from plotter.printer import syseffects_printer_ender

#### import default analysis histogram edges
from analyzer.misc import default_edges as edges
from analyzer.misc import greco_nyears, seconds_per_year

scalar = greco_nyears * seconds_per_year

#########################################################
#### define global variables
#########################################################
### directory where preweighted pickled
### data are located

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
settings = {'cmapname':'RdBu'                   ,
            'reverse_cmap':True                  ,
            'stretch_cscd':True                 ,
            'norm'    :True                     ,
            'xedges'  :np.log10 (edges['e'])    ,
            'yedges'  :np.cos (edges['z'])[::-1],
            'zedges'  :edges['p']               ,
            'xticks'  :np.log10 ([5.6, 10, 18, 32, 56]),
            'yticks'  :[-1.0, -0.5, 0.0, 0.5, 1.0],
            'xticklabels':[6, 10, 18, 32, 56],
            'yticklabels':[-1.0, -0.5, 0.0, 0.5, 1.0] }  

#########################################################
#### set plot format
#########################################################
### global format
format_global     = {'tick_fontsize' :13 ,
                     'label_fontsize':15 ,
                     'title_fontsize':17 ,
                     'grid_alpha'    :0.2,
                     'grid_linestyle':None,
                     'grid_linewidth':0.  ,
                     'colorbar_tick_length':4 }

### 2D histograms
format_histograms = {'hist_alpha'         :0.9  ,
                     'hist_print'         :False,
                     'hist_print_rounding':1    ,
                     'hist_print_fontsize':3.0   }

### color bar
format_colorbar = {'colorbar_alpha_dm31'      :0.9      ,
                   'colorbar_format_dm31'     :'%10.1f' ,
                   'colorbar_ranges_dm31'     :(-2, 2)  ,
                   'colorbar_nbins_dm31'      :5        ,

                   'colorbar_alpha_axm_res'      :0.9     ,
                   'colorbar_format_axm_res'     :'%10.1f',
                   'colorbar_ranges_axm_res'     :(-4, 4),
                   'colorbar_nbins_axm_res'      :5       ,

                   'colorbar_alpha_barr_nubar_ratio'      :0.9      ,
                   'colorbar_format_barr_nubar_ratio'     :'%10.1f' ,
                   'colorbar_ranges_barr_nubar_ratio'     :(-8, 8),
                   'colorbar_nbins_barr_nubar_ratio'      :5        ,

                   'colorbar_alpha_nue_numu_ratio'      :0.9     ,
                   'colorbar_format_nue_numu_ratio'     :'%10.1f',
                   'colorbar_ranges_nue_numu_ratio'     :(-1, 1) ,
                   'colorbar_nbins_nue_numu_ratio'      :5       ,

                   'colorbar_alpha_forward'      :0.9     ,
                   'colorbar_format_forward'     :'%10.1f',
                   'colorbar_ranges_forward'     :(-6, 6),
                   'colorbar_nbins_forward'      :5        }

### Put all formats into a list
plot_formats = [settings, format_global, 
                format_histograms,
                format_colorbar  ]

#########################################################
### Functions to plot via `HistogramND`
#########################################################
def plot_syseffects (outdir, events, sample, stretch=False):

    ''' plot all variables in events

        :param outdir (str) : user's output directory
        :param events (dict): all event variables needed
    '''

    ### define parameters to be plotted
    ### based on 'param'_weight
    parameters = [ '_'.join (key.split ('_')[:-1])
                   for key in events[events.keys ()[0]].keys ()
                   if 'weight' in key and
                   not 'seeded' in key ]

    ### define a histogram instance
    ### with all global plot formats 
    histo = HistogramND (sample, parameters, stretch=stretch)

    ### set observable keys
    histo.sample = sample
    histo.xobs = 'reco_logen'
    histo.yobs = 'reco_coszen'
    histo.zobs = 'deltallh' if sample=='dragon' else \
                 'tracklength'

    ### set plot styling / formatting
    for formats in plot_formats:
        histo.set_style (**formats)

    ## plot and save !
    histo.plot (events)
    outname = sample+'_syseffects.pdf'
    histo.save (outdir+'/'+outname)
    return

#########################################################
#### functions to know which sample it is
#########################################################
def which_sample (infile):

    ### catch the sample name from infile
    for sam in ['greco', 'dragon']:
        if sam in infile: return sam

    ### if not in infile, ask user
    question = 'May I ask which sample it is ? [greco/dragon]'
    sample = None
    ### check validity
    while sample not in ['greco', 'dragon']:
        ## if not valid, ask again
        sample = raw_input (question)
    
    return sample

#########################################################
#### functions to collect events
#########################################################
def collect_events (infile):

    ''' collect event variables as needed

        :return events (dict): all event variables needed
    '''

    with open (infile, 'rb') as f:
        events = cPickle.load (f)
    f.close ()
    
    ### print progress
    syseffects_printer_events (events)
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
    parser.add_option ("--infile", type="string", default=None,
                       help = "must provide an input file with all events")
    (options, args) = parser.parse_args()
    outdir  = options.outdir
    infile  = options.infile

    #### print header
    syseffects_printer_header (outdir)

    #### get events
    events = collect_events (infile)

    #### get which sample it is
    sample = which_sample (infile)

    #### plot variables
    plot_syseffects (outdir, events, sample,
                     stretch=settings['stretch_cscd'])

    #### print ender
    syseffects_printer_ender (outdir)

