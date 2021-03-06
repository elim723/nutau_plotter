#!/usr/bin/env python

####
#### By Elim Thompson (09/28/2018)
####
#### This script made DRAGON 2D selection
#### plots for the containment cuts.
####
#### 1. Starting position at Level 6
####     xaxis = PegLeg Rho
####     yaxis = PegLeg Z
####     color = summed neutrino rates
####
#### Command to run
####  $ python plot_dragon2D.py --outdir <plot>
####
###############################################


#### import standard packages
from __future__ import print_function
from optparse import OptionParser
from glob import glob
import numpy as np
import cPickle, os

from plotter.defaults import labels
from plotter.histogram2D import Histogram2D
from plotter.printer import selection_printer_header
from plotter.printer import selection_printer_events
from plotter.printer import selection_printer_ender

#########################################################
#### define global variables
#########################################################
### directory where preweighted pickled
### data are located
sample = 'dragon'
ppath  = '/data/user/elims/nutau_pickled_data/dragon/preweighted/'

#########################################################
#### define variables related to PegLeg cuts
#### vline / hline / line are lists of dictionaries
####
#### DO NOT TOUCH THIS BLOCK
#### unless you know what you are doing
#########################################################
pegleg_title   = r'$\mu$ fraction at Level 6'

## xvalues for the cut line
xvalues = np.linspace (0, 100, 1000)
yvalues = np.sqrt (125**2 - xvalues**2) - 350-75

## diagonal line
pegleg_line  = [{'xvalues':xvalues,
                 'yvalues':yvalues}]

## vertical cut line
pegleg_vline = [{'xvalue' :100,
                 'yranges':(60/700., 350/700.)}]

## horizontal cut lines
pegleg_hline = [{'yvalue' :-475,
                 'xranges':(0, 100/200.)}]

#########################################################
#### define plot settings
####
#### NOTE: To make plots of other variabls, add a
####       dictionary to the settings dictionary
####       like below. And make sure you have the
####       same default key names as the ones in
####       resources/plotter/defaults.py.
####
#### For each 2D plot,
####   key   : levelX/xvariable/yvariable
####   values: 
####      xranges  = x axis limits
####      yranges  = y axis limits
####      colorbar_ranges  = z axis (colorbar) limits
####      cmapname = color map name (google cmap names)
####      xnbins   = number of bins in x axis
####      ynbins   = number of bins in y axis
####      colorbar_nbins   = number of ticks in color bar
####      title    = title of the plot
####      htype    = histogram type:
####                 1. numucc = distribution of numu CC
####                             (same for any member)
####                 2. nurate = total neutrino rates
####                 3. mufrac = muon fraction
####                 4. s/b    = nutau / sqrt (not nutau)
####       NOTE: to implement new htype, go to
####             `plotter/histogram2D/_get_content ()`
####      cut_vline = any vertical lines to be plotted
####                  set to None if no cut
####                  type = list of dictionaries
####      cut_hline = any horizontal lines to be plotted
####                  set to None if no cut
####      cut_line  = any lines to be plotted
####                  set to None if no cut
#########################################################
settings = {'level6/start_rho/start_z':{'xranges'        :(   0,  200),
                                        'yranges'        :(-500, -200),
                                        'colorbar_ranges':(0, 0.8)    ,
                                        'cmapname'       :'Blues'     ,
                                        'xnbins'         :150         ,
                                        'ynbins'         :150         ,
                                        'colorbar_nbins' :9           ,
                                        'title'          :pegleg_title,
                                        'htype'          :'mufrac'    ,
                                        'cut_vline'      :pegleg_vline,
                                        'cut_hline'      :pegleg_hline,
                                        'cut_line'       :pegleg_line  } }

#########################################################
#### set plot styling / formating
#########################################################
### Global 
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
format_colorbar   = {'colorbar_alpha'      :0.7              ,
                     'colorbar_format'     :'%10.2f'         ,
                     'colorbar_tick_length':2                ,
                     'colorbar_label'      :r'$\mu$ fraction' }

### format cut lines
format_lines      = {'cut_alpha'    :0.9  ,
                     'cut_color'    :'red',
                     'cut_linewidth':3.0  ,
                     'cut_linestyle':'-'   }

### Put all formats into a list
plot_formats = [format_global, format_lines,
                format_histograms, format_colorbar]

#########################################################
### Functions to plot via `Histogram2D`
#########################################################
def get_additional_formats (variable):

    formats = settings[variable]

    ## set x / ylabels
    level, xvar, yvar = variable.split ('/', 2)
    formats['hist_xlabel'] = labels[xvar]
    formats['hist_ylabel'] = labels[yvar]

    ## set specific histogram x ticks
    formats['hist_xticks'] = np.linspace (formats['xranges'][0],
                                          formats['xranges'][1], 11)
    ## set specific histogram y ticks
    formats['hist_yticks'] = np.linspace (formats['yranges'][0],
                                          formats['yranges'][1], 7)

    return formats

def plot_variables (outdir, events):

    ''' plot all variables in events

        :param outdir (str) : user's output directory
        :param events (dict): all event variables needed
    '''

    ### define a histogram instance
    ### with all global plot formats 
    histo = Histogram2D ()

    ### set plot styling / formatting
    for formats in plot_formats:
        histo.set_style (**formats)

    ## loop through each variable
    for variable in settings:
        ## define level and variable
        level, xvar, yvar = variable.split ('/', 2)

        ## set any format for this variable plot
        formats = get_additional_formats (variable)
        histo.set_style (**formats)

        ## plot and save !
        histo.plot (events[level], xvar, yvar)
        outname = sample+'_'+level+'_'+xvar+'_'+yvar+'.pdf'
        histo.save (outdir+'/'+outname)
    return

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
        level, xvar, yvar = variable.split ('/', 2)
        ## create new key if not has one
        if not level in variables:
            variables[level] = []
        ## add to variables dictionary
        variables[level].append (xvar)
        variables[level].append (yvar)
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

######################################
#### Everything starts here :)
######################################
if __name__ == "__main__":

    #### parse options/params
    usage = "%prog [--outdir <directory for plots>]"
    parser = OptionParser(usage=usage)
    parser.add_option( "-o", "--outdir", type="string", default='~/',
                       help = "out directory for plots")
    (options, args) = parser.parse_args()
    outdir  = options.outdir

    #### print header
    selection_printer_header (outdir, sample.upper ())

    #### get events
    events = collect_events ()

    #### plot variables
    plot_variables (outdir, events)

    #### print ender
    selection_printer_ender (outdir)

