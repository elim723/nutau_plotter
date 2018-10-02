#!/usr/bin/env python

####
#### By Elim Thompson (09/23/2018)
####
#### This script plot all six PegLeg resolution
#### plots for nutau paper. The six output plots
#### are:
####       A. Energy Resolution
####          1. log10  reco E vs log10  true E: "log10e_log10e"
####          2. linear reco E vs linear true E: "lineare_lineare"
####       B. Zenith Resolution
####          1. cos zenith   vs log10  true E: "log10e_cosinez"
####          2. zenith (deg) vs log10  true E: "log10e_degreez"
####          3. cos zenith   vs linear true E: "lineare_cosinez"
####          4. zenith (deg) vs iinear true E: "lineare_degreez"
####
#### command to run:
####
#### option 1 with an infile
####    $ python plot_resolutions.py
####             --infile pegleg_resolutions.p
####
#### option 2 without an infile
####    $ python plot_resolutions.py
####             --grecodir pickled_data/greco/final/
####             --dragondir pickled_data/dragon/final/
####
##################################################################

from __future__ import print_function
from optparse import OptionParser
import cPickle, sys, os
from glob import glob
import numpy as np

## custum functions
from plotter.resolution import Resolution
from plotter.printer import resolution_printer_header
from plotter.printer import resolution_printer_events
from plotter.printer import resolution_printer_rates
from plotter.printer import resolution_printer_res
from plotter.printer import resolution_printer_ender

from analyzer.nuparams import Nuparams
from analyzer.library import Library

## ignore RuntimeWarning
import warnings
warnings.filterwarnings("ignore")

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
path  = '/data/user/elims/nutau_pickled_data/'
grecodir  = path + 'greco/final/'
dragondir = path + 'dragon/final/'
thispath  = os.path.dirname(os.path.abspath( __file__ ))
nufile    = thispath + '/../events_collecter/nufiles/nuparams_template.txt'

#### members   = individual member included
#### datatypes = the data types to be plotted
####             NC neutrinos will be merged to nunc
members   = ['numucc', 'nuecc', 'nutaucc',
             'numunc', 'nuenc', 'nutaunc']
datatypes = ['numucc', 'nuecc', 'nutaucc', 'nunc']

#### import a toolbox to merge NC into nunc
from analyzer.misc import Toolbox
toolbox = Toolbox ()

#### this is where the MC truth energy bins
#### (x-axis) are pre-defined
binedges = {'lineare':np.linspace (0, 60, 25),
            'log10e' :np.linspace (0.2, 1.8, 18)}

#### define axis combinations (xaxis_yaxis)
combinations = ['lineare_lineare',
                'log10e_log10e'  ,
                'lineare_cosinez',
                'lineare_degreez',
                'log10e_cosinez' ,
                'log10e_degreez' ]

################################################
#### set plotting style
################################################
#### global styles for all six combinations
format_global = {'alpha'           :0.8    ,
                 'linewidth'       :5.0    ,
                 'median_linestyle':'-'    ,
                 'sigma_linestyle' :':'    ,
                 'greco_color'     :'red'  ,
                 'dragon_color'    :'green',
                 'greco_label'     :r'$\mathcal{A}$',
                 'dragon_label'    :r'$\mathcal{B}$',
                 ### axes fontsizes
                 'label_fontsize'  :38     ,
                 'tick_fontsize'   :30     ,
                 ### legend
                 'legend_fontsize' :22     ,
                 'legend_alpha'    :1.0    ,
                 ### text
                 'text_fontsize'   :40     ,
                 'text_frameon'    :False  ,
                 ### grid 
                 'grid_alpha'      :0.1    ,
                 'grid_linewidth'  :0.3    ,
                 'grid_linestyle'  :'-'    ,
                 ### reference lines
                 'ref_plot'        :True   ,
                 'ref_alpha'       :0.8    ,
                 'ref_color'       :'gray' ,
                 'ref_linewidth'   :1.0    ,
                 'ref_linestyle'   :'-'     }

#### legend styles for each combinations
### "legend_'xaxis'_'yaxis'_plot" = which plot has the legend
### "legend_'xaxis'_'yaxis'_loc"  = location of legend
### "legend_'xaxis'_'yaxis'_ncols"= number of columns of the legend
format_legend = {'legend_lineare_lineare_plot' :'numucc',
                 'legend_lineare_lineare_loc'  :2       ,
                 'legend_lineare_lineare_ncols':1       ,
                 
                 'legend_log10e_log10e_plot'   :'numucc',
                 'legend_log10e_log10e_loc'    :4       ,
                 'legend_log10e_log10e_ncols'  :1       ,

                 'legend_log10e_cosinez_plot'  :'numucc',
                 'legend_log10e_cosinez_loc'   :1       ,
                 'legend_log10e_cosinez_ncols' :2       ,

                 'legend_log10e_degreez_plot'  :'numucc',
                 'legend_log10e_degreez_loc'   :1       ,
                 'legend_log10e_degreez_ncols' :2       ,

                 'legend_lineare_cosinez_plot' :'numucc',
                 'legend_lineare_cosinez_loc'  :1       ,
                 'legend_lineare_cosinez_ncols':2       ,

                 'legend_lineare_degreez_plot' :'numucc',
                 'legend_lineare_degreez_loc'  :1       ,
                 'legend_lineare_degreez_ncols':2        }

#### text styles for each combinations
### "text_'xaxis'_'yaxis'_'datatype'_loc" = location of text
format_text = {'text_lineare_lineare_numucc_loc' :4,
               'text_lineare_lineare_nuecc_loc'  :4,
               'text_lineare_lineare_nutaucc_loc':2,
               'text_lineare_lineare_nunc_loc'   :2,

               'text_log10e_log10e_numucc_loc'   :2,
               'text_log10e_log10e_nuecc_loc'    :2,
               'text_log10e_log10e_nutaucc_loc'  :2,
               'text_log10e_log10e_nunc_loc'     :2,

               'text_log10e_cosinez_numucc_loc'  :4,
               'text_log10e_cosinez_nuecc_loc'   :4,
               'text_log10e_cosinez_nutaucc_loc' :4,
               'text_log10e_cosinez_nunc_loc'    :4,

               'text_log10e_degreez_numucc_loc'  :4,
               'text_log10e_degreez_nuecc_loc'   :4,
               'text_log10e_degreez_nutaucc_loc' :2,
               'text_log10e_degreez_nunc_loc'    :2,

               'text_lineare_cosinez_numucc_loc'  :4,
               'text_lineare_cosinez_nuecc_loc'   :4,
               'text_lineare_cosinez_nutaucc_loc' :4,
               'text_lineare_cosinez_nunc_loc'    :4,

               'text_lineare_degreez_numucc_loc'  :4,
               'text_lineare_degreez_nuecc_loc'   :4,
               'text_lineare_degreez_nutaucc_loc' :2,
               'text_lineare_degreez_nunc_loc'    :2 }

#### axis styles for each combination
### "'xaxis'_xticks"      = x axis ticks
### "'xaxis'_xlabel"      = x axis label
### "'xaxis'_xranges"     = x axis ranges
### "'xaxis'_xtickslabel" = x axis tick labels
### "'yaxis'_yticks"      = y axis ticks
### "'yaxis'_ylabel"      = y axis label
### "'yaxis'_yranges"     = y axis ranges
### "'yaxis'_ytickslabel" = y axis tick labels

lineare_ticks = np.linspace (0, 60, 7)
log10e_ticks  = np.linspace (0.2, 1.8, 5)
cosinez_ticks = np.linspace (-1, 1, 5)
degreez_ticks = np.linspace (-80, 80, 5)

format_axes = {'lineare_xticks'     :lineare_ticks                        ,
               'lineare_xlabel'     :r'MC truth E$_{\nu}$ (GeV)'          ,
               'lineare_xranges'    :(lineare_ticks[0], lineare_ticks[-1]),
               'lineare_xtickslabel':['%.0f' % e for e in lineare_ticks]  ,

               'log10e_xticks'     :log10e_ticks                          ,
               'log10e_xlabel'     :r'MC truth E$_{\nu}$ (GeV)'           ,
               'log10e_xranges'    :(0.1, 1.9)                            ,
               'log10e_xtickslabel':['%.0f' % 10**e for e in log10e_ticks],
               
               'lineare_yticks'     :lineare_ticks                        ,
               'lineare_ylabel'     :r'Reconstructed E$_{\nu}$ (GeV)'     ,
               'lineare_yranges'    :(lineare_ticks[0], lineare_ticks[-1]),
               'lineare_ytickslabel':['%.0f' % e for e in lineare_ticks]  ,

               'log10e_yticks'     :log10e_ticks                          ,
               'log10e_ylabel'     :r'Reconstructed E$_{\nu}$ (GeV)'     ,
               'log10e_yranges'    :(0.1, 1.9)                            ,
               'log10e_ytickslabel':['%.0f' % 10**e for e in log10e_ticks],

               'cosinez_yticks'     :cosinez_ticks                                             ,
               'cosinez_ylabel'     :r'cos $\theta_{\text{reco}}$ - cos $\theta_{\text{true}}$',
               'cosinez_yranges'    :(cosinez_ticks[0], cosinez_ticks[-1])                     ,
               'cosinez_ytickslabel':['%.1f' % z for z in cosinez_ticks]                       ,
               
               'degreez_yticks'     :degreez_ticks                                                  ,
               'degreez_ylabel'     :r'$\theta_{\text{reco}}$ - $\theta_{\text{true}}$ ($^{\circ}$)',
               'degreez_yranges'    :(degreez_ticks[0], degreez_ticks[-1])                          ,
               'degreez_ytickslabel':['%.0f' % z for z in degreez_ticks]                             }
               
plot_formats = [format_global, format_axes, format_legend, format_text]

################################################
#### functions to plot pegleg resolutions
################################################
def plot_resolutions (outdir, resolutions):

    ### define canvas for resolution plot
    res = Resolution ()
    ### set plot styling / formatting
    for formats in plot_formats:
        res.set_style (**formats)

    ### loop through each x/y combinations
    for combin, resolution in resolutions.items ():
        ### set new xaxis / yaxis
        res.xaxis, res.yaxis = combin.split ('_')
        ### plot
        res.plot (resolution)
        ### save
        outname = 'resolutions_'+combin+'.pdf'
        res.save (outdir+'/'+outname)
    return

##################################################
#### functions to get events
##################################################
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
    nuparams = Nuparams (nufile, isinverted=inverted)

    ### set up weighters for baseline
    ### members using injected values
    injected = nuparams.extract_params ('injected')
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
        terms of data types, which are the actual data types
        get plotted. By default, members contain individual
        NC and CC neutrinos (total of 6), whereas datatypes
        contain 3 CC and 1 merged NC neutrinos (total of 4).

        :param memevents (dict): dictionary from 6 individual
                                 members

        :return devents (dict): dictionary from 4 data types
                                (3 CC and 1 merged NC)
    '''

    devents = {}

    ### loop through numucc, nuecc, nutaucc, nunc
    for dtype in datatypes:
        ## does not change CC
        if 'cc' in dtype.lower ():
            devents[dtype] = memevents[dtype]
            continue

        ## merge if NC
        devents[dtype] = toolbox.concat_dicts (memevents['numunc'],
                                               memevents['nuenc'],
                                               memevents['nutaunc'])
    return devents

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
        lib = Library (members, ppath, edges=edges,
                       isdragon=sample=='dragon',
                       verbose=0)
        ## get weights
        events_in_members = get_weights (lib)

        ## merge nc
        events[sample] = merge_nc (events_in_members)

    return events

##################################################
#### functions to get energy / zenith resolutions
##################################################
def statistics (num, w):

    ''' get the reconstructed values at lower 1sigma
        50% middle, and upper 1sigma at a given bin

        :param num (array): reconstructed values at an ebin
        :param w   (array): weights for each values

        :return sigma_lower  (float): xvalue at lower 1sigma
        :return sigma_middle (float): xvalue at 50%
        :return sigma_upper  (float): xvalue at upper 1sigma
    '''
    
    ### if empty array, return 0s
    if len (num) == 0: return 0, 0, 0

    ### sort xvalues
    indicies = np.argsort (num)
    num = np.array ([ num[i] for i in indicies ])
    w = np.array ([ w[i] for i in indicies ])

    ### find the indices at lower 1 sigmas,
    ### middle, and upper 1 sigmas
    cumulative = np.cumsum (w) / np.sum (w)
    # lower part is at 50%-34% = 16%
    lower = (cumulative > 0.16)
    # median
    middle = (cumulative > 0.50)
    #upper is 50%+34% = 84%
    upper = (cumulative > 0.84)

    ### x values at those three points
    sigma_lower  = num[lower][0]
    sigma_middle = num[middle][0]
    sigma_upper  = num[upper][0]
    
    return sigma_lower, sigma_middle, sigma_upper

def get_statistics (yvalues, xvalues, weights, binedge):

    ''' get the lower 1sigma / middle 50% / upper 1sigma
        for each xvalue bins

        :param yvalues (array): reconstructed values
                                either linear / log10 energy
                                or (reco_deg_zenith - true_deg_zenith)
                                or (reco_cos_zenith - true_cos_zenith)
        :param xvalues (array): MC truth energy
        :param weights (array): event weights
        :param binedge (array): energy bin edges
    
        :return rstats (dict): lower / upper 1sigma
                               and middle for each ebin
                               {'xvalues': energy bin centers
                                'lsigma' : lower 1 sigma
                                'median' : median
                                'usigma' : upper 1 sigma}
    '''
    
    stats = {}
    ### loop through each energy bin
    for index, bine in enumerate (binedge[:-1]):
        ## energy range of this bin
        lowere = binedge[index]
        uppere = binedge[index+1]
        bincenter = (lowere + uppere) * 0.5
        
        ## boolean whether an element of xvalues
        ## (MC truth energy) is in this energy range
        thisbin = np.logical_and (xvalues>=lowere, xvalues<uppere)

        ## get lower 1sigma, middle, and upper 1sigma
        if len (thisbin[thisbin]) == 0:
            ## ignore if no elements in this bin,
            continue

        ## if there are events in this bin,
        ## get their yvalues and weights
        yvalue = yvalues[thisbin]
        weight = weights[thisbin]
        ## get the statistics
        lower, middle, upper = statistics (yvalue, weight)
            
        ## store to dictionary
        stats[bincenter] = {'lower':lower, 'middle':middle, 'upper':upper}
        
    ### rearrage stats
    ##  sort bincenters
    bincenters = sorted (stats.keys ())
    rstats = {'xvalues':bincenters,
              'lsigma' :np.array ([ stats[bc]['lower']  for bc in bincenters ]),
              'median' :np.array ([ stats[bc]['middle'] for bc in bincenters ]),
              'usigma' :np.array ([ stats[bc]['upper']  for bc in bincenters ]) }
    return rstats
            
def get_resolutions (event, yaxis=None, xaxis=None):

    ''' calculate energy / zenith resolutions for a
        given data type of a given the xaxis and yaxis
        scaling from both samples

        :param event (dict): events of one data type from one sample
        :param yaxis (str) : 6 optional scalings of y axis
                             reconstructed energy or zenith
                             'lineare', 'log10e', 'cosinez', 'degreez'
        :param xaxis (str) : 2 optional scalings of x axis
                             MC truth energy
                             'lineare', 'log10e'

        :return res (dict): resolutions of energy / zenith
                            as a function of MC truth energy
    '''

    ## define variables needed
    #  1. xvalues = mc truth energy
    xvalues = event.mc.e if xaxis=='lineare' else \
              np.log10 (event.mc.e) ## xaxis=='log10e'
    #  2. yvalues = reconstructed energy or
    #               (reco - true) zenith
    yvalues = event.reco.e if yaxis=='lineare' else \
              np.log10 (event.reco.e) if yaxis=='log10e' else \
              event.reco.cz - event.mc.cz if yaxis=='cosinez' else \
              (event.reco.z - event.mc.z) * 180./np.pi  ## yaxis=='degreez'
    #  3. weights
    weights = event.w
    #  4. energy bin edges
    binedge = binedges[xaxis]

    ## get y values at upper 1sigma / 50% / lower 1sigma
    return get_statistics (yvalues, xvalues, weights, binedge)

def get_all_resolutions (combin, events):

    ''' get resolutions from all four channels (numucc,
        nuecc, nutaucc, nunc) from all six combinations
         -- linear_reco_e vs linear_mc_e
         -- log_reco_e    vs log_mc_e
         -- reco - true cosine_zenith vs linear_mc_e
         -- reco - true cosine_zenith vs log10_mc_e
         -- reco - true degree_zenith vs linear_mc_e
         -- reco - true degree_zenith vs log10_mc_e
       for events from one sample

       :param events (dict): all events from one sample

       :return allres (dict): all resolutions for this sample
    '''

    ### holder to store resolutions for
    ### one of the possible combinations
    ### from all members for both samples
    allres = {}

    ### the x/y scaling
    xaxis, yaxis = combin.split ('_')

    ### loop through each data type
    for dtype in datatypes:
        ## holder for resolutions of a given data type
        allres[dtype] = {}
        ## get and store resolutions from a sample
        for sample in samples:
            event = events[sample][dtype]
            allres[dtype][sample] = get_resolutions (event,
                                                     yaxis=yaxis,
                                                     xaxis=xaxis)
            resolution_printer_res (combin, sample, dtype,
                                    allres[dtype][sample])
    return allres

def collect_resolutions (infile, outdir, gdir, ddir):

    ''' collect resolutions either from an input
        `infile` if provided or calculate resolutions                                                                                         
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
                resolutions = cPickle.load (f)
            f.close ()
            return resolutions
        except:
            pass

    ### report to user
    outfile = outdir + '/pegleg_resolutions.p'
    resolution_printer_events (outfile)

    ### get events from scratch
    events = get_events (gdir, ddir)
    resolution_printer_rates (events)

    ### get resolutions for each sample
    ### and all 6 options
    resolutions = {combin:get_all_resolutions (combin, events)
                   for combin in combinations}

    ### save resolutions for future
    with open (outfile, 'wb') as f:
        cPickle.dump (resolutions, f, protocol=2)
    f.close ()

    return resolutions

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
    resolution_printer_header (outdir)

    ### collect events from either infile or from scratch
    ##  The structure of resolutions:
    ##  {combin: {dtype: {sample: {'xvalues':[],
    ##                             'lsigma' :[],
    ##                             'usigma' :[] }}}}
    resolutions = collect_resolutions (infile, outdir, grecodir, dragondir)
    
    ### plot
    plot_resolutions (outdir, resolutions)

    ### print ender
    resolution_printer_ender (outdir)
