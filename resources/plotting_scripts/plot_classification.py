#!/usr/bin/env python

####
#### By Elim Thompson (09/21/2018)
####
#### This script plot the PID plot in the reconstruction
#### section of the nutau paper. This plot shows the
#### fraction of track-like events as a function of
#### energy.
#### energy bins. These data points are then passed
#### to the Classification class in `plotter/classification.py`
#### for plotting.
####
#### By default, the left plot is from GRECO (analysis A),
#### and right plot is from DRAGON (analysis B). Only
#### neutrino are included and are separated by
#### interaction types.
####
#### command to run:
####
#### option 1 with an infile
####    $ python plot_classification.py
####             --infile pid_track_fractions.p
####
#### option 2 without an infile
####    $ python plot_classification.py
####             --grecodir pickled_data/greco/final/
####             --dragondir pickled_data/dragon/final/
####
##############################################################

from __future__ import print_function
from optparse import OptionParser
import cPickle, sys, os
from glob import glob
import numpy as np

## custum analysis tool
from analyzer import nuparams, library

## custum functions
from plotter.classification import Classification 
from plotter.printer import classification_printer_header
from plotter.printer import classification_printer_events
from plotter.printer import classification_printer_rates
from plotter.printer import classification_printer_ender

## this append is needed for loading
## files when track fractions need
## to be recalculated
sys.path.append ('../../analyzer/')

## ignore RuntimeWarning
import warnings
warnings.filterwarnings ("ignore")

################################################
#### default variables
################################################
from plotter.styler import labels
from analyzer.misc import default_edges as edges

#### set oscillation setting in case
#### no infile is provided
matter, inverted, oscnc = True, False, True

#### default paths to greco and dragon
#### samples at their final levels
samples   = ['greco', 'dragon']
path      = '/data/user/elims/'
grecodir  = path + 'nutau_pickled_data/greco/final/'
dragondir = path + 'nutau_pickled_data/dragon/final/'
thispath  = os.path.dirname(os.path.abspath( __file__ ))
nufile    = thispath + '/../events_collecter/nufiles/nuparams_template.txt'

#### members = individual member included
#### lines   = the lines to be plotted
####           NC neutrinos will be merged to nunc
members   = ['numucc', 'nuecc', 'nutaucc',
             'numunc', 'nuenc', 'nutaunc']
lines = ['numucc', 'nuecc', 'nutaucc', 'nunc']

#### import a toolbox to merge NC into nunc
from analyzer.misc import Toolbox
toolbox = Toolbox ()

#### this is where the energy bins are pre-defined
loge_edges = np.linspace (0.75, 1.75, 25)

################################################
#### set plotting style
################################################
format_global = {'alpha'         :0.8,
                 'linewidth'     :5.0,
                 'linestyle'     :'-',
                 'label_fontsize':33 ,
                 'tick_fontsize' :27 ,
                 'grid_alpha'    :0.2, 
                 'grid_linewidth':0.5,
                 'grid_linestyle':'-' }

format_axes = {'ylabel'     :'Fraction of events \n classified as track'   ,
               'xlabel'     :labels['mc_e']                                ,
               'xticks'     :np.linspace (loge_edges[0], loge_edges[-1], 9),
               'yticks'     :np.linspace (0.0, 1.0, 6)                     ,
               'xtickslabel':[6, 8, 10, 13, 18, 24, 32, 42, 56]        ,
               'ytickslabel':np.round (np.linspace (0.0, 1.0, 6), 1)       ,
               'xranges'    :(0.7, 1.75)                                   ,
               'yranges'    :(0.0, 1.0 )                                   }

format_legend = {'legend_fontsize':30 , 
                 'legend_loc'     :2  ,
                 'legend_ncols'   :1  ,
                 'legend_alpha'   :1.0 }

format_text = {'text_fontsize'  :33   ,
               'text_frameon'   :False,
               'text_greco_loc' :1    ,
               'text_dragon_loc':1     }

plot_formats = [format_global, format_axes, format_legend, format_text]

################################################
#### functions to plot pid classification
################################################
def plot_pid (outdir, track_fractions):

    ### define canvas for pid plot
    pid = Classification ()
    pid.samples    = samples
    pid.loge_edges = loge_edges
    ### set plot styling / formatting
    for formats in plot_formats:
        pid.set_style (**formats)

    ### plot and save !
    pid.plot (track_fractions)
    outname = 'pid.pdf'
    pid.save (outdir+'/'+outname)
    return

################################################
#### functions to get fraction of track events
################################################
def get_weights (lib):

    ''' calculat weights for all members

        :param lib (`Libraray`): all events and weighters
                                 for all members

        :return dictionary (dict): dictionary with all events
                                   and weights in terms of member
                                   (i.e. individual NC)
    '''
    
    ### define parameter values for weighting
    ### from `nufiles/nuparam_textfile.txt`
    nuparam = nuparams.Nuparams (nufile, isinverted=inverted)
    
    ### set up weighters for baseline
    ### members using injected values
    injected = nuparam.extract_params ('injected')
    lib.set_weighters (injected, matter=matter, oscnc=oscnc)

    ### loop through each member
    ### and get all event variables
    dictionary = {}
    for dtype in members:
        # define member / weighter
        member, weighter = lib (dtype), lib.weighters[dtype]
        event = member ()
        # get weights based on injected values
        event.w = member.get_weights (injected,
                                      weighter=weighter)
        dictionary[dtype] = event

    ### return events in terms of members
    return dictionary

def merge_nc (memevents):
    
    ''' convert events in terms of members to events in
        terms of lines, which are the actual data types
        get plotted. By default, members contain individual
        NC and CC neutrinos (total of 6), whereas lines
        contain 3 CC and 1 merged NC neutrinos (total of 4).
        
        :param memevents (dict): dictionary from 6 individual
                                 members
        
        :return linevents (dict): dictionary from 4 data types
                                  (3 CC and 1 merged NC)
    '''
    
    linevents = {}
    
    ### loop through numucc, nuecc, nutaucc, nunc
    for dtype in lines:
        ## does not change CC
        if 'cc' in dtype.lower ():
            linevents[dtype] = memevents[dtype]
            continue

        ## merge if NC
        linevents[dtype] = toolbox.concat_dicts (memevents['numunc'],
                                                 memevents['nuenc'],
                                                 memevents['nutaunc'])
    return linevents
        
def get_events (grecodir, dragondir):
    
    ''' get all events from both samples calculate their
        weights as well

        :param grecodir  (str): folder of GRECO final level
        :param dragondir (str): folder of DRAGON final level

        :return events (dict): all events
    '''
    
    events = {}
    for sample in samples:
        ## define ppath where pickled data are
        ppath = eval (sample + 'dir')
        ## initialize library instance
        lib = library.Library (members, ppath, edges=edges,
                               isdragon=sample=='dragon',
                               verbose=0)
        ## get weights
        events_in_members = get_weights (lib)
        
        ## merge nc
        events[sample] = merge_nc (events_in_members)

    return events

def classify_events (sample, events):

    ''' classify events into track-like / cascade-like
        of all energy bins in log10 scale to obtain the
        track fraction per energy bin.

        :param sample (str) : greco or dragon
        :param events (dict): all events

        :return fractions (dict): track fractions per member
                                  per energy bin from a sample
    '''
    
    ### define holders for track fractions from all members
    fractions = {}
    ### define pid key
    ## 'pid'  for tracklength in GRECO
    ## 'dllh' for delta LLH in DRAGON
    pkey = 'dllh' if sample=='dragon' else 'pid'

    ### loop through each member to get track fractions
    ### as a function of log10 MC truth energy
    for member, event in events.items ():
        ## define holder to append fractions
        fractions[member] = []
        ## define variables needed
        mcloge, weight = np.log10 (event.mc.e), event.w
        recopid = event.reco[pkey]
        # define boolean of `is track-like events`
        istrack = recopid>=50. if sample=='greco' else recopid>=2.
        
        ## loop through log10 energy bins
        for index, loge in enumerate (loge_edges[:-1]):
            # select events with MC truth energy
            # inside this loge bin
            thisbin = np.logical_and (mcloge>=loge_edges[index],
                                      mcloge<loge_edges[index+1])
            # if no events fall in this bin,
            # append None to fractions
            if len (thisbin[thisbin]) == 0:
                fractions[member].append (None)
                continue
            # if there are events in this bin, calculate
            # ratio of track-like events in this bin to
            # total number of events in this bin
            boolean = np.logical_and (thisbin, istrack)
            trackf = np.sum (weight[boolean]) / np.sum (weight[thisbin])
            # append to fractions
            fractions[member].append (trackf)
        ## numpify array
        fractions[member] = np.array (fractions[member])
    return fractions
    
def collect_track_fractions (infile, outdir, gdir, ddir):

    ''' collect track fractions either from an input
        `infile` if provided or calculate track fractions
        from scratch. If built from scratch, the output
        is stored to `outdir` for future plots.

        :param infile (str): input track fraction dictionary
        :param outdir (str): location of output track fractions
                             will be stored if infile not provided
        :param gdir (str): location of GRECO final level sample
        :param ddir (str): location of DRAGON final level sample
      
        :return track_fracs (dict): track fractions as a function
                                    of MC truth energy for each
                                    member from both samples
    '''
    
    ### if user provides a valid input file,
    ### load it (save time!)
    if infile:
        try: 
            with open (infile, 'rb') as f:
                track_fracs = cPickle.load (f)
            f.close ()
            return track_fracs
        except:
            pass

    ### report to user
    outfile = outdir + '/pid_track_fractions.p'
    classification_printer_events (outfile)
    
    ### get events from scratch
    events = get_events (gdir, ddir)
    classification_printer_rates (events)

    ### classify events for each sample
    track_fracs = {sample:classify_events (sample, events[sample])
                   for sample in samples}

    ### save classigied events for future
    with open (outfile, 'wb') as f:
        cPickle.dump (track_fracs, f, protocol=2)
    f.close ()

    ### return classified events
    return track_fracs
    
#############################################
#### main function
#############################################
if __name__ == "__main__":

    ### parse user's options
    usage = "%prog [--outdir plots/ --infile classified.p]"
    parser = OptionParser(usage=usage)
    parser.add_option( "-o", "--outdir", type="string", default='~/',
                       help = "out directory for plotting")
    parser.add_option( "-i", "--infile", type="string", default=None,
                       help = "pre-calculated track-like fractions")
    parser.add_option( "-d", "--dragondir", type="string", default=dragondir,
                       help = "path to DRAGON final data sets")
    parser.add_option( "-g", "--grecodir", type="string", default=grecodir,
                       help = "path to GRECO final data sets")
    (options, args) = parser.parse_args()
    outdir = options.outdir
    infile = options.infile
    grecodir  = options.grecodir
    dragondir = options.dragondir

    ### print header
    classification_printer_header (outdir)
    
    ### collect events from either infile or from scratch
    ### This is where the fraction of track-like events
    ### in each energy bin is calculated.
    track_fracs = collect_track_fractions (infile, outdir, grecodir, dragondir)

    ### plot PID via Classification
    plot_pid (outdir, track_fracs)

    ### print ender
    classification_printer_ender (outdir)
