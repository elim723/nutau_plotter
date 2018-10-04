#!/usr/bin/env python

####
#### By Elim Thompson (09/23/2018)
####
#### This script plot the 90% contours of standard
#### oscillation parameter measurements done by
#### various experiments, including greco and dragon.
####
#### command to run:
#### $ python plot_greco_numu.py --oudir <output folder>
####
#########################################################

from __future__ import print_function
from optparse import OptionParser
import cPickle, sys, os
from glob import glob
import numpy as np

## import printers
from plotter.printer import contour_printer_header
from plotter.printer import contour_printer_ender

## custum classes / functions
from plotter.contour import Contour

################################################
#### default variables
################################################
#### filename to all numu 90% contours
numu90s_filename = '/data/user/elims/nutau_pickled_data/numu_contours.p'

################################################
#### plot formating / styling
################################################
#### all colors / linestyles / labels are
#### defined in `plotter/defaults.py`

format_global = {'alpha'         :0.55,
                 'linewidth'     :3.0,
                 'tick_fontsize' :17 ,
                 'label_fontsize':20 ,
                 'grid_alpha'    :0.2,
                 'grid_linewidth':0.5,
                 'grid_linestyle':'-' }

format_marker = {'marker'          :'x',
                 'marker_size'     :40 ,
                 'marker_alpha'    :1.0,
                 'marker_linewidth':3.0 }

format_legend = {'legend_fontsize':17 ,
                 'legend_loc'     :1  ,
                 'legend_ncols'   :3  ,
                 'legend_alpha'   :1.0 }

format_axes = {'dm_ticks'    :np.linspace (2, 3.2, 7),
               'dm_ranges'   :(1.9, 3.2),
               'theta_ticks' :np.linspace (0.3, 0.7, 5),
               'theta_ranges':(0.25, 0.75),
               'chi2_ticks'  :np.linspace (0, 5, 6),
               'chi2_ranges' :(0, 5)                }

plot_formats = [format_global, format_marker, format_legend, format_axes]

########################################
#### function to plot
########################################
def plot_contours (contours, bestfits, dm_scans, theta_scans, outdir):

    ### define canvas for contour plot
    contour = Contour ()
    ### set plot styling / formatting
    for formats in plot_formats:
        contour.set_style (**formats)

    ### plot and save !
    contour.plot (contours, bestfits,
                  dm_scans, theta_scans)
    outname = 'greco_numu.pdf'
    contour.save (outdir+'/'+outname)
    return

########################################
#### function to collect data points
########################################
def collect_contours ():

    #### contour data points are stored
    #### in a pickled dictionary
    with open (numu90s_filename, 'rb') as f:
        numu90s = cPickle.load (f)
    f.close ()

    #### separate data points
    ### 1. contours    = {exp:{'x':sin2theta  , 'y':dm  }}
    ### 2. dm_scans    = {exp:{'x':chi2       , 'y':dm  }}
    ### 3. theta_scans = {exp:{'x':sin2theta  , 'y':chi2}}
    ### 4. bestfit     = {exp:{'x':bfsin2theta, 'y': bfdm}
    contours, bestfits    = {}, {}
    dm_scans, theta_scans = {}, {}
    
    ### loop through each experiment
    for exp, datapoints in numu90s.items ():
        ## collect 90% contour
        contours[exp] = numu90s[exp]['contour']
        
        ## collect scans if available
        if 'dm23' in numu90s[exp]:
            dm_scans[exp] = numu90s[exp]['dm23']
        if 'theta23' in numu90s[exp]:
            theta_scans[exp] = numu90s[exp]['theta23']
        
        ## collect best fit points if available
        if 'bestfit' in numu90s[exp]:
            bfpoints = numu90s[exp]['bestfit']
            # change to dictionary if its not
            if type (bfpoints) in [tuple, list, np.array]:
                bfpoints = {'x':bfpoints[0], 'y':bfpoints[1]}
            bestfits[exp] = bfpoints

    ### return all dictionaries
    return contours, bestfits, dm_scans, theta_scans

#############################################
#### main function
#############################################
if __name__ == "__main__":
7
    ### parse user's options
    usage = "%prog [--outdir plots/]"
    parser = OptionParser (usage=usage)
    parser.add_option( "-o", "--outdir", type = "string", default = '~/',
                       help = "out put path that the contour plot will be saved")
    (options, args) = parser.parse_args()
    outdir = options.outdir

    ### print header
    contour_printer_header (outdir)
    
    ### collect data points
    contours, bestfits, dm_scans, theta_scans = collect_contours ()
    
    ### plot contours
    plot_contours (contours, bestfits, dm_scans, theta_scans, outdir)

    ### print ender
    contour_printer_ender (outdir)
