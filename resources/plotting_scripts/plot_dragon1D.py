#!/usr/bin/env python

####
#### By Elim Thompson (09/25/2018)
####
#### This script plots DRAGON 1D selection plots
#### for the combined nutau paper. The variables
#### included in the paper includes
####
#### level 5: total_charge, separation, bdt_score
#### level 6: start_z, start_rho, stop_z, stop_rho
####
#### Command line:
####   python plot_dragon1D.py --outdir <output folder>
####
#### NOTE: These plots are made from Elim's
####       pickled files at lower levels with
####       an extension `.pckl`, located at
####       resources/pickled_data/dragon/preweighted/.
####       These pickled files have weights
####       pre-calculated, which are used for
####       plotting. A pickling script is also
####       provided in case one wants to repickle
####       the lower level i3 files.
####
####       No Level4 variables are plotted because
####       FeiFei's data hdf5 files do not have 
####       pre-level4 information. And, Elim is
####       too lazy to repickle the data files
####       herself. (Sorry.. too much to do now..)
####
#########################################################

#### import standard packages
from __future__ import print_function
from optparse import OptionGroup, OptionParser
from glob import glob
import numpy as np
import cPickle, os

from plotter.defaults import greco_scalar, dragon_scalar
from plotter.histogram1D import Histogram1D
from plotter.fitter1D import Fitter1D
from plotter.printer import selection_printer_header
from plotter.printer import selection_printer_events
from plotter.printer import selection_printer_ender

#########################################################
#### define global variables
#########################################################
### directory where preweighted pickled
### data are located
sample = 'dragon'
scalar = dragon_scalar
ppath  = '/data/user/elims/nutau_pickled_data/dragon/preweighted/'

#########################################################
#### define plot settings
####
#### NOTE: To make plots of other variabls, add a
####       dictionary to the settings dictionary
####       like below. And make sure you have the
####       same default key names as the ones in
####       resources/plotter/defaults.py.
####
#### For each variable,
####   key   : levelX_variable_name
####   values: 
####      xranges = x axis limits
####      yranges = y axis limits on distribution (top) plot
####      rranges = y axis limits on ratio (bottom) plot
####      NOTE: if above ranges are changed, make sure
####            you check their ticks in `get_additional_formats`
#### 
####      nbins   = number of bins in x axis
####      hloc    = location of legend box in histogram plot
####                (1 = upper right , 2 = upper left  ,
####                 3 = lower left  , 4 = lower right ,
####                 8 = lower center, 9 = upper center,
####                 google the rest)
####      rloc    = location of legend box in ratio plot
####      ncols   = number of columns in the legend box
####      cut     = cut value from selection
####                (None for those included in BDTs)
####
#########################################################
settings = {'level5_total_charge'  : {'xranges':(   0,  250)  ,
                                      'yranges':(1e-7, 1e-3)  ,
                                      'rranges':(0.0 ,  2.0)  ,
                                      'nbins'  :50            ,
                                      'hloc'   :1 , 'rloc': 3 ,
                                      'ncols'  :2 , 'cut':None},
            'level5_separation'    : {'xranges':(0   ,   220) ,
                                      'yranges':(1e-5,  0.001),
                                      'rranges':(0   , 2)     ,
                                      'nbins'  :50            , 
                                      'hloc'   :1 , 'rloc': 9 ,
                                      'ncols'  :2 , 'cut':None},
            'level5_bdt_score'     : {'xranges':(-0.4, 0.7)   ,
                                      'yranges':(1e-6, 1e-4)  ,
                                      'rranges':(0.4, 1.6)    ,
                                      'nbins'  :50            ,
                                      'hloc'   :1 , 'rloc': 8 ,
                                      'ncols'  :1 , 'cut':0.2 },
            'level6_start_z'       : {'xranges':(-550, -150)  ,
                                      'yranges':(1e-8, 10e-5) ,
                                      'rranges':(0.5, 1.5)    ,
                                      'nbins'  :50            , 
                                      'hloc'   :1 , 'rloc': 8 ,
                                      'ncols'  :2 , 'cut':None},
            'level6_start_rho'     : {'xranges':(   0,  250)  ,
                                      'yranges':(1e-6, 1e-5)  ,
                                      'rranges':(0.5, 1.5)    ,
                                      'nbins'  :50            ,
                                      'hloc'   :1 , 'rloc': 3 ,
                                      'ncols'  :1 , 'cut':None},
            'level6_stop_z'        : {'xranges':(-600, -100)  ,
                                      'yranges':(1e-7, 0.0001),
                                      'rranges':(0.0, 2.0)    ,
                                      'nbins'  :50            ,
                                      'hloc'   :1 , 'rloc':8  ,
                                      'ncols'  :2 , 'cut':[-500, -200]},
            'level6_stop_rho'      : {'xranges':(   0,  250)  ,
                                      'yranges':(1e-6, 1e-5)  ,
                                      'rranges':(0.0, 2.0)    ,
                                      'nbins'  :50            ,
                                      'hloc'   :1 , 'rloc':3  ,
                                      'ncols'  :1 , 'cut':150 } }

############################################################################
#### define plot formatting
####
#### These variables are to change the plots layouts.
#### Here show all the possible changes one can make.
#### Elim tried to make it as flexible as possible,
#### so you don't have to dig into the code. But,
#### to implement new adjustment, go to `histogram1D`.
####
#### NOTE: hatches linewidth / color are globally
####       defined in `plotter/styler.py` via rcParams.
############################################################################

### Global:
##  - plottotalmc         (bool): if True, plot total MC histogram
##  - hist/ratio_ylabel    (str): label for y axes
##  - hist/ratio_logy     (bool): log y axis if True
##  - hist/ratio_legend_fontsize (float): font sizes of legend for hist / ratio
##  - hist/ratio_legend_alpha    (float): transparency of legend frame (0-1)
##  - tick/label_fontsize        (float): fontsizes for tick / label / legend
##  - grid_alpha/linestyle/linewidth    : grid format

format_global = {'plottotalmc'          :False           ,
                 'hist_ylabel'          :r'Event Rate Hz',
                 'hist_logy'            :True            ,
                 'hist_legend_fontsize' :13              ,
                 'hist_legend_alpha'    :1.0             ,
                 'ratio_ylabel'         :r'Ratio to MC'  ,
                 'ratio_logy'           :False           ,
                 'ratio_legend_fontsize':11              ,
                 'ratio_legend_alpha'   :1.0             ,
                 'tick_fontsize'        :11              ,
                 'label_fontsize'       :13              ,
                 'grid_alpha'           :0.2             ,
                 'grid_linestyle'       :'-'             ,
                 'grid_linewidth'       :0.5             }

### Global format for ratio reference horizontal line
### and cut value vertical line
##  - ratioref/cut_alpha     (float): transparancy of plot between 0 and 1
##  - ratioref/cut_color     (str)  : color of the line
##  - ratioref/cut_linewidth (float): line widths of the lines
##  - ratioref/cut_linestyle (str)  : standard matplotlib linestyle (:/--/-/etc)
format_lines = {'ratioref_alpha'    :0.8    , 'cut_alpha'    :0.8   ,
                'ratioref_color'    :'black', 'cut_color'    :'blue',
                'ratioref_linewidth':2.0    , 'cut_linewidth':2.0   ,
                'ratioref_linestyle':'-'    , 'cut_linestyle':'-'   }

### Total MC / Stacked histograms:
##  - mc/stack_linewidth (float): line widths of the histogram 
##  - mc/stack_linestyle (str)  : standard matplotlib linestyle (:/--/-/etc)
##  - mc/stack_alpha     (float): transparancy of plot between 0 and 1

format_histograms = {'stack_linewidth':0.8 , 'mc_linewidth':2  , 
                     'stack_linestyle':None, 'mc_linestyle':'-', 
                     'stack_alpha'    :0.65, 'mc_alpha'    :0.6 }                     

### Data points / Total MC errorbars / Stacked errorbars /
### Data/MC ratio points / MC uncertainty errorbars :
##  - alpha, linestyle, linewidth
##  - capsize, capthick, elinewidth,
##  - marker, markersize, markeredgewidth, fillstyle

## for data points on upper histogram
format_data_points = {'data_alpha'          :0.7   ,
                      'data_linestyle'      :None  ,
                      'data_linewidth'      :0.0   ,
                      'data_capsize'        :0.0   ,
                      'data_capthick'       :0.0   ,
                      'data_elinewidth'     :3.0   ,
                      'data_marker'         :'o'   ,
                      'data_markersize'     :5.0   ,
                      'data_markeredgewidth':0.0   ,
                      'data_fillstyle'      :'full'}

## for total and stacked MC histograms on upper histogram
format_hist_errors = {'mcerr_alpha'          :0.6 , 'stackerr_alpha'          :0.6 ,
                      'mcerr_linestyle'      :None, 'stackerr_linestyle'      :None,
                      'mcerr_linewidth'      :0.0 , 'stackerr_linewidth'      :0.0 ,
                      'mcerr_capsize'        :0.0 , 'stackerr_capsize'        :0.0 ,
                      'mcerr_capthick'       :0.0 , 'stackerr_capthick'       :0.0 ,
                      'mcerr_elinewidth'     :3.0 , 'stackerr_elinewidth'     :3.0 ,
                      'mcerr_marker'         :'o' , 'stackerr_marker'         :'o' ,
                      'mcerr_markersize'     :0.0 , 'stackerr_markersize'     :0.0 ,
                      'mcerr_markeredgewidth':0.0 , 'stackerr_markeredgewidth':0.0 ,
                      'mcerr_fillstyle'      :None, 'stackerr_fillstyle'      :None}

## for data/mc ratio points ('dataratio') and mc
## uncertainty bands ('mcratio') on lower ratio plot
format_ratios = {'dataratio_alpha'          :0.9   , 'mcratio_alpha'          :0.6 ,
                 'dataratio_linestyle'      :None  , 'mcratio_linestyle'      :None,
                 'dataratio_linewidth'      :0.0   , 'mcratio_linewidth'      :0.0 ,
                 'dataratio_capsize'        :0.0   , 'mcratio_capsize'        :0.0 ,
                 'dataratio_capthick'       :0.0   , 'mcratio_capthick'       :0.0 ,
                 'dataratio_elinewidth'     :3.0   , 'mcratio_elinewidth'     :7.5 ,
                 'dataratio_marker'         :None  , 'mcratio_marker'         :None,
                 'dataratio_markersize'     :0.0   , 'mcratio_markersize'     :0.0 ,
                 'dataratio_markeredgewidth':0.0   , 'mcratio_markeredgewidth':0.0 ,
                 'dataratio_fillstyle'      :'full', 'mcratio_fillstyle'      :None}

### Put all formats into a list
plot_formats = [format_global, format_lines,
                format_histograms , format_data_points,
                format_hist_errors, format_ratios ]

#########################################################
#### functions to collect events
#########################################################
def get_variables ():

    ''' collect variables from the global
        setting dictionary

        :return variables (dict): variable names at each level
                                  {level:[var1, var2]}
    '''

    variables = {}
    ### loop through each key in settings
    for variable in settings:
        level, var = variable.split ('_', 1)
        ## add to variables dictionary
        if level in variables:
            ## if level in variables, append
            variables[level].append (var)
        else:
            ## else create new key / value
            variables[level] = [var]
    return variables

def get_events (level, filename, variable):

    ''' get events

        :param level     (str): at which level it is for cut
        :param filename  (str): path to the pickled file
        :param variable (list): variables needed for plots

        :return event (dict): variables / weights needed
    '''
    
    ### open file
    with open (filename, 'rb') as f:
        ddict = cPickle.load (f)
    f.close ()
    
    ### apply cuts
    cutkey = '_bool'
    cut = np.ones (len (ddict[level + cutkey])).astype (bool)
    ## uncomment this line if actually apply cut
    #cut = ddict[level + cutkey] 
    
    ### collect events as needed
    ## weight in Hz!
    event = {'weight':ddict['weight'].flatten ()[cut]}
    for key in ddict.keys ():
        if key.lower () in variable:
            event[key.lower ()] = ddict[key][cut]
    return event

def collect_events ():

    ''' collect event variables as needed

        :return events (dict): all event variables needed
    '''
    
    ### holder for all events needed
    events = {}
    ### get {level:[var1, var2]} dictionary
    variables = get_variables ()

    ### loop through each level
    for level, variable in variables.items ():
        ## holder for events at this level
        events [level] = {}
        ## filenames/members at this level
        filenames = sorted (glob (ppath + level + '*.pckl'))
        ## temperary fix in filenames
        if level=='level6':
            filenames = sorted (glob (ppath + 'level6_postreco_*.pckl'))
        ## loop through each filename (i.e. each member)
        for filename in filenames:
            member = os.path.split (filename)[1].split ('.')[0].split ('_')[-1]
            events[level][member] = get_events (level, filename, variable)
    ### print progress 
    selection_printer_events (events)
    return events

###########################################
### Functions to plot via `Histogram1D`
###########################################
def get_additional_formats (variable):
    
    formats = settings[variable]
    ## add prefix to x label to show level
    formats['xlabel_prefix'] = 'Level ' + variable.split ('_', 1)[0][-1]
    ## set specific ratio y ticks
    formats['ratio_yticks'] = np.linspace (formats['rranges'][0],
                                           formats['rranges'][1], 5)
    ## set specific histogram y ticks
    #formats['hist_yticks'] = np.linspace (formats['yranges'][0],
    #                                      formats['yranges'][1], 5)
    ## set specific histogram x ticks
    ## (automatically applied to ratio x ticks
    #formats['hist_xticks'] = np.linspace (formats['xranges'][0],
    #                                      formats['xranges'][1], 6)
    return formats

def plot_variables (outdir, events, fitnorm=False):

    ''' plot all variables in events
        
        :param outdir (str) : user's output directory
        :param events (dict): all event variables needed
        :param fitnorm (bool): If True, fit events
    '''
    
    ### define a histogram instance
    ### with all global plot formats
    histo = Histogram1D ()
    ### set scalar for MC error bars
    histo.scalar = scalar
    ### set plot styling / formatting
    for formats in plot_formats:
        histo.set_style (**formats)
        
    ### create a fitter instance
    fitter = Fitter1D (events, scalar)

    ## loop through each variable
    for variable in settings:
        ## define level and variable
        level, var = variable.split ('_', 1)
        
        ## set any format for this variable plot
        formats = get_additional_formats (variable)
        histo.set_style (**formats)
            
        ## add ratio refence line
        #refratio = {'ratioref_yvalue':1.0}
        #histo.set_style (**refratio)

        event = events[level]
        ## fit data / mc via norms if asked
        if fitnorm:
            fitter.xranges = formats['xranges']
            fitter.nbins   = formats['nbins']
            fitter.dtypes  = sorted (event.keys ())
            event = fitter.fit_events (level, var)
        
        ## plot and save !
        histo.plot (event, var)
        outname = sample+'_'+level+'_'+var+'.pdf' 
        histo.save (outdir+'/'+outname)
    return

######################################
#### Everything starts here :)
######################################
if __name__ == "__main__":

    #### parse options/params 
    usage = "%prog [--outdir <directory for plots>]"
    parser = OptionParser(usage=usage)
    parser.add_option( "-o", "--outdir", type="string", default='~/',
                       help = "out directory for plots")
    parser.add_option ("--fitnorm", action="store_true", default=False,
                       help = "If True, fit normalizations.")
    (options, args) = parser.parse_args()
    outdir  = options.outdir
    fitnorm = options.fitnorm

    #### print header
    selection_printer_header (outdir, sample.upper ())

    #### get events
    events = collect_events ()

    #### plot variables
    plot_variables (outdir, events, fitnorm=fitnorm)

    #### print ender
    selection_printer_ender (outdir)
