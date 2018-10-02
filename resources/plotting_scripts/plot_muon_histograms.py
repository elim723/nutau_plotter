#!/usr/bin/env python

####
#### By Elim Thompson (09/28/2018)
####
#### This script made the muon histograms
#### from both GRECO and DRAGON samples.
#### By default, 2 x 2 histograms are made
#### from both samples. If a sample is
#### specified, 1 x 2 histograms are made.
####
#### Command to run
####  $ python plot_muon_histograms.py
####          --outdir <plot>
####          (--dragon if DRAGON only OR
####           --greco  if GRECO  only   )
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
from plotter.histogram6D import Histogram6D
from plotter.printer import muhist_printer_header
from plotter.printer import muhist_printer_events
from plotter.printer import muhist_printer_ender

#### import default analysis histogram edges
from analyzer.misc import default_dragon_edges as dragon_edges
from analyzer.misc import default_edges as greco_edges
from analyzer.misc import greco_nyears, dragon_nyears
from analyzer.misc import seconds_per_year

scalars = {'greco' : greco_nyears * seconds_per_year,
           'dragon': dragon_nyears * seconds_per_year}

#########################################################
#### define global variables
#########################################################
### directory where preweighted pickled
### data are located
ppath  = '/data/user/elims/nutau_pickled_data/'
filenames = {'greco' :ppath+'greco/preweighted/level7_muongun.pckl',
             'dragon':ppath+'dragon/preweighted/level6_final_muon.pckl'}

variables = {'greco' :['reco_logen', 'reco_coszen', 'tracklength'],
             'dragon':['reco_logen', 'reco_coszen', 'deltallh']   }

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
settings = {'htype'      :'muon'                             ,
            'xedges'     :np.log10 (greco_edges['e'])        ,
            'xticks'     :np.log10 ([5.6, 10, 18, 32, 56]) ,
            'xticklabels':[6, 10, 18, 32, 56] }

greco_settings = {'cmapname'   :'Oranges'                         ,
                  'stretch_cscd':False                            ,
                  'colorbar_ranges':(0, 200)                      ,
                  'yedges'      :np.cos (greco_edges['z'])[::-1]  ,
                  'yticks'      :[-1.0, -0.5, 0.0, 0.5, 1.0],
                  'yticklabels' :[-1.0, -0.5, 0.0, 0.5, 1.0],
                  'zedges'      :greco_edges['p']                  }

dragon_settings = {'cmapname'   :'Reds'                             ,
                   'stretch_cscd':False                           ,
                   'colorbar_ranges':(0, 80)       ,
                   'yedges'      :np.cos (dragon_edges['z'])[::-1],
                   'yticks'      :[-1.0, -0.5, 0.0, 0.5, 1.0]     ,
                   'yticklabels' :[-1.0, -0.5, 0.0, 0.5, 1.0]     ,
                   'zedges'      :dragon_edges['p']                }

#########################################################
#### set plot format
#########################################################
### global format
format_global     = {'tick_fontsize' :11 ,
                     'label_fontsize':13 ,
                     'title_fontsize':13 ,
                     'grid_alpha'    :0.2,
                     'grid_linestyle':None,
                     'grid_linewidth':0. }

### 2D histograms
format_histograms = {'hist_alpha'         :0.9  ,
                     'hist_print'         :False,
                     'hist_print_rounding':1    ,
                     'hist_print_fontsize':6.0   }

### color bar
format_colorbar   = {'colorbar_alpha'      :0.9            ,
                     'colorbar_format'     :'%10.0f'       ,
                     'colorbar_tick_length':0              ,
                     'colorbar_label'      :r'Number of $\mu$ events',
                     'colorbar_nbins'      :6               }

### Put all formats into a list
plot_formats = [settings, format_global, 
                format_histograms,
                format_colorbar  ]

#########################################################
### Functions to plot via `Histogram3D`
#########################################################
def plot_muhist (outdir, events):

    ### if more than a samples, 2 x 2
    ### histograms are made via Histogram6D
    if len (events) > 1:
        plot_histo2x2 (outdir, events)
        return

    ### only 1 sample, only
    ### Histogram3D is needed
    plot_histo1x2 (outdir, events)
    return

def plot_histo2x2 (outdir, events):

    ### define a histogram instance
    ### with all global plot formats 
    histo = Histogram6D ()

    for sample, event in events.items ():

        ### set greco / dragon settings
        extra_settings = eval (sample + '_settings')
    
        ### set observable keys
        histo.xobs = 'reco_logen'
        histo.yobs = 'reco_coszen'
        histo.zobs = 'deltallh' if sample=='dragon' else \
                     'tracklength'
    
        ### set plot styling / formatting
        for formats in plot_formats + [extra_settings]:
            histo.set_style (**formats)

        ## plotq
        histo.plot (sample, event)

    outname = 'both_muon_histogram.pdf'
    histo.save (outdir+outname)
    return

def plot_histo1x2 (outdir, events):

    ### define sample
    sample = events.keys ()[0]

    ### set greco / dragon settings
    extra_settings = eval (sample + '_settings')
    
    ### if GRECO sample,
    ### want to stretch cascade bin?
    stretch_cscd = False if sample=='dragon' else \
                   extra_settings['stretch_cscd']
    
    ### define a histogram instance
    ### with all global plot formats 
    histo = Histogram3D (sample, stretch=stretch_cscd)
    
    ### set observable keys
    histo.xobs = 'reco_logen'
    histo.yobs = 'reco_coszen'
    histo.zobs = 'deltallh' if sample=='dragon' else \
                 'tracklength'
    
    ### set plot styling / formatting
    for formats in plot_formats + [extra_settings]:
        histo.set_style (**formats)

    ## plot and save !
    histo.plot (events[sample])
    outname = sample+'_muon_histogram.pdf'
    histo.save (outdir+outname)
    return

#########################################################
#### functions to collect events
#########################################################
def get_weight (weight, sample):

    ### For DRAGON only,
    ##  the weight from the pickled files from
    ##  is not quite right. Since muon icc are
    ##  equally weighted, the total rate is
    ##  rescaled to the muon rate in the nutau
    ##  paper.
    if sample=='dragon':
        length = len (weight)
        rate   = 0.0259/1000.
        weight = np.ones (length) * rate / float (length)

    ### rescale to detector livetime
    weight *= scalars[sample]
    
    ### return
    return weight

def get_events (filename, sample):

    ''' get events

        :param filename  (str): path to the pickled file

        :return event (dict): variables / weights needed
    '''

    ### open file
    with open (filename, 'rb') as f:
        ddict = cPickle.load (f)
    f.close ()

    ### apply cuts
    cutkey = 'level7_bins' if sample == 'greco' else \
             'level6_bool'
    ##  additional cut within the histogram
    cut  = ddict['tracklength'] < 1000. if sample == 'greco' else \
           np.logical_and (ddict['deltaLLH'] >= -3 , ddict['deltaLLH'] < 1000.)
    ecut = np.logical_and (ddict['reco_logen']>=0.75, ddict['reco_logen']<1.75)
    cut  = np.logical_and (cut, ecut)
    cut  = np.logical_and (cut, ddict[cutkey].astype (bool))

    ### collect events as needed
    ##  weight in Hz!
    weight = ddict['weight'].flatten ()[cut]
    weight = get_weight (weight, sample)

    ### initialize an event dictionary
    ### with only muon member
    event = {'muon':{'weight':weight}}

    ###  loop through each key 
    for key in ddict.keys ():
        ## only collect variables needed
        if key.lower () in variables[sample]: 
            # keyname must be lower cases to match
            # the dictionaries in `default.py`
            keyname = key.lower ()
            # get the variables
            event['muon'][keyname] = ddict[key][cut]

    ## return all muon events at final levels
    return event

def collect_events (samples, infiles):

    ''' collect event variables as needed

        :return events (dict): all event variables needed
    '''

    ### holder for all events needed
    events = {}

    ### loop through each sample
    for sample in samples:

        ### use input file if available
        infile = infiles[sample]
        if infile:
            with open (infile, 'rb') as f:
                event = cPickle.load (f)
            f.close ()
            ## switch PID key
            events[sample] = switchPID (event)
            ## modify weights
            mweight = get_weight (events[sample]['muon']['weight'], sample)
            events[sample]['muon']['weight'] = mweight
            continue

        ### define muon filename
        filename = filenames[sample]
        ## get muon events
        event = get_events (filename, sample)
        ## switch PID key
        events[sample] = switchPID (event)

    ### print progress
    muhist_printer_events (events)
    return events

def switchPID (ddict):
    
    ### loop through keys
    for key in ddict.keys ():
        ## pick the pid key
        if key in ['deltaLLH', 'tracklength']:
            # create a new pid key
            ddict['pid'] = ddict[key]
    return ddict

######################################
#### Everything starts here :)
######################################
if __name__ == "__main__":

    #### parse options/params
    usage = "%prog [--outdir <directory for plots>]"
    parser = OptionParser(usage=usage)
    parser.add_option ("--outdir", type="string", default='~/',
                       help = "out directory for plots")
    parser.add_option ("--greco_infile", type="string", default=None,
                       help = "in file for greco sample")
    parser.add_option ("--dragon_infile", type="string", default=None,
                       help = "in file for dragon sample")
    parser.add_option ("--dragon", action="store_true", default=False,
                       help = "DRAGON muon histogram only")
    parser.add_option ("--greco", action="store_true", default=False,
                       help = "GRECO muon histogram only")
    (options, args) = parser.parse_args()
    outdir = options.outdir
    dragon = options.dragon
    greco  = options.greco
    infiles = {'greco':options.greco_infile,
               'dragon':options.dragon_infile}

    #### define what samples to be used
    samples = ['greco']  if greco and not dragon else \
              ['dragon'] if dragon and not greco else \
              ['greco', 'dragon']

    #### print header
    muhist_printer_header (outdir, samples)

    #### get events
    events = collect_events (samples, infiles)

    #### plot variables
    plot_muhist (outdir, events)

    #### print ender
    muhist_printer_ender (outdir)

