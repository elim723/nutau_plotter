#!/usr/bin/env python

####
#### By Elim Cheung (07/24/2018)
####
#### This script contains the plotter class that makes
#### plotting easier for a numu disappearance analysis.
####
#### The actual plotting scripts to run are located in
#### plot_scripts folder.
####
####################################################################

from __future__ import print_function
from copy import deepcopy
import numpy as np

######################################
#### plot settings
######################################
from styler import colors, hatches, labels
from styler import plt, gridspec
from styler import matplotlib

from analyzer.misc import default_edges as edges

######################################
#### standard plotting canvas class
######################################
class Canvas (object):

    ''' a class that define canvas for various nutau plots '''

    def __init__ (self, figsize, dimensions, 
                  width_ratios=None,
                  height_ratios=None,
                  wspace=0.1, hspace=0.1,
                  top=None, bottom=None,
                  right=None, left=None):

        ''' initialize canvas by its dimensions / spacing / ratios

            :param figsize (tuple): canvas width / height in inches
            :param dimensions (tuple): number of subplots in row and column
            :param width_ratios (list): ratio of the width of subplots
            :param height_ratios (list): ratio of the height of subplots
            :param wspace (float): spacing in width between plots
            :param hspace (float): spacing in height between plots
            :param top/bottom/right/left (float): extent of the subplots as a fraction
                                                  of figure width or height
                                                  (left cannot be larger than right,
                                                   bottom cannot be larger than top)
        '''
        
        self.h  = plt.figure (figsize=figsize)
        self.gs = gridspec.GridSpec (dimensions[0], dimensions[1],
                                     width_ratios=width_ratios,
                                     height_ratios=height_ratios)
        self.gs.update (wspace=wspace, hspace=hspace, top=top,
                        bottom=bottom, right=right, left=left)

    def set_style (self, **kwargs):

        for key, value in kwargs.iteritems ():
            setattr (self, key, value)

    @staticmethod
    def get_range (*contents):

        ''' obtain the range from at least one content

            :param content (list): content of a 1D / 2D histogram
        '''

        xmin, xmax = np.inf, -np.inf
        for content in contents:
            array = np.nan_to_num (content).flatten ()
            if xmin > np.min (array):
                xmin = np.min (array)
            if np.max (array) > xmax:
                xmax = np.max (array)
        return xmin, xmax

    @staticmethod
    def _get_bin_center (edge):

        ''' return bin center from edge

            :param edge (np.array): bin edges
            :return bin_center (np.array): center of each bin
        '''

        return np.array ([ edge[i] + (edge[i+1]-edge[i])/2.
                           for i in range (len (edge)-1) ])
    
    def get_power10s (self, power10x=False, power10y=False):

        ''' print tick labels as power of 10 

            :param power10x (bool): if True, print 10**self.edges['x'] instead
            :param power10y (bool): if True, print 10**self.edges['y'] instead

            :return eticks (dict): tick labels to be printed
                                   {'x':[], 'y':[]}
        '''
        
        eticks = deepcopy (self.edges)
        for a in ['x', 'y']:
            edge = np.power (10, self.edges[a]) if eval ('power10'+a) else self.edges[a]
            fmt = '%.1f' if eval ('power10'+a) else '%.2f'
            eticks[a] = [fmt % x for x in edge]
        return eticks
    
    def save (self, name, title=None):

        ''' save current canvas to name
        
            :param name  (str): full path (/address/name.png or .pdf)
            :param title (str): super title of the canvas
        '''

        if title: plt.suptitle (title, fontsize=self.labelsize)        
        self.h.savefig (name)
        plt.close('all')
        return

##########################################
#### a canvas class specific for 2D plots
##########################################
class Canvas2D (Canvas):

    def __init__ (self, figsize, dimensions,
                  width_ratios=None,
                  height_ratios=None,
                  wspace=0.1, hspace=0.1,
                  top=None, bottom=None,
                  right=None, left=None,
                  **kwargs):

        ### default 2D plot arragements
        super (Canvas2D, self).__init__ (figsize, dimensions,
                                         width_ratios=width_ratios,
                                         height_ratios=height_ratios,
                                         wspace=wspace, hspace=hspace,
                                         top=top, bottom=bottom,
                                         right=right, left=left)

        ### set style (e.g. labelsize/ticksize/etc)
        self.set_style (**kwargs)

    def _set_edges (self):
            
        ### define color map
        self.cmap = plt.get_cmap (self.cmapname)

        ### define x / y / z edges
        for axis in ['x', 'y', 'z']:
            
            ## only has edges if `ranges` is defined
            ## by user
            ranges = getattr (self, axis+'ranges', -1)
            if ranges==-1: continue

            ## define nbins for an axis
            ## default = 9 bins
            nbins = getattr (self, axis+'nbins', 9)

            ## define edges
            edge  = np.linspace (ranges[0], ranges[1], nbins)

            ## store it to `self`
            setattr (self, axis+'edges', edge)

        return
        
    def _plot2D (self, axis, content, xedges=[], yedges=[]):

        ''' plot 2D content

            :param axis (`matplotlib.axes`): current axis for plotting
            :param content (np.array): content per bin
        '''
        
        ### transpose between numpy and matplotlib
        mH = np.ma.masked_array (content).T
        mH.mask = np.logical_not (mH)

        ## define x/y edges
        xedges = xedges if len (xedges) > 0 else self.xedges
        yedges = yedges if len (yedges) > 0 else self.yedges

        axis.pcolormesh (xedges, yedges, mH,
                         vmin=self.colorbar_ranges[0],
                         vmax=self.colorbar_ranges[1],
                         alpha=self.hist_alpha,
                         cmap=self.cmap)

        ### print bin content if asked
        if self.hist_print: self._print2D (axis, content)

    def _print2D (self, axis, content):

        ''' print 2D content per bin

            :param content (np.array): content per bin
        '''

        rounding = self.hist_print_rounding
        fontsize = self.hist_print_fontsize

        ### loop through each bin
        for xdex, xedge in enumerate (self.xedges[:-1]):
            for ydex, yedge in enumerate (self.yedges[:-1]):
                ## rounding the value in this bin
                value = round (content[xdex][ydex], rounding)
                ## set text and its color
                text = '-' if not np.isfinite (value) else str (value)
                darker = value > self.colorbar_ranges[1]*0.85 or \
                         value < self.colorbar_ranges[0]*0.85
                tcolor = 'white' if text=='-' or darker else 'black'
                ## set the coordinate of the printout
                x = self.xedges[xdex] + (self.xedges[xdex+1]-self.xedges[xdex])/8.
                y = self.yedges[ydex] + (self.yedges[ydex+1]-self.yedges[ydex])/3.
                ## actually print
                axis.annotate (text, xy=(x, y),
                               xycoords='data',
                               color=tcolor   ,
                               fontsize=fontsize)

    def _format2D (self, axis, xticks, yticks, xlabel, ylabel,
                   xticklabels=[], yticklabels=[]):

        ''' format 2D axis

            :param axis (`Matplotlib.Axes`): axis to be formatted
            :param title (str): title of the subplot
        '''

        #### set x, y limits
        axis.set_xlim (xticks[0], xticks[-1])
        axis.set_ylim (yticks[0], yticks[-1])

        #### set x, y ticks
        axis.set_xticks (xticks)
        axis.set_yticks (yticks)
        axis.tick_params (axis='x', labelsize=self.tick_fontsize)
        axis.tick_params (axis='y', labelsize=self.tick_fontsize)

        #### set x, y ticklabels
        xticklabels = xticklabels if len (xticklabels) > 0 else \
                      None if xticklabels == [] else \
                      getattr (self, 'hist_xticklabels', xticks)
        yticklabels = yticklabels if len (yticklabels) > 0 else \
                      None if yticklabels == [] else \
                      getattr (self, 'hist_yticklabels', yticks)

        if xticklabels:
            axis.set_xticklabels (xticklabels)
        else: 
            axis.get_xaxis ().set_ticks ([])

        if yticklabels:
            axis.set_yticklabels (yticklabels)
        else: 
            axis.get_yaxis ().set_ticks ([])

        #### set grid lines
        for xtick in xticks:
            axis.axvline (x=xtick,
                          color='grey',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)
        for ytick in yticks:
            axis.axhline (y=ytick,
                          color='grey',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)

        #### set x, y labels and title
        if xlabel: axis.set_xlabel (xlabel, fontsize=self.label_fontsize)
        if ylabel: axis.set_ylabel (ylabel, fontsize=self.label_fontsize)
        return
                
    def _plot_colorbar (self, axis):

        ''' plot colorbar

            :param axis (`matplotlib.axes`): the axis things to be plotted
        '''
        
        ### normalize colorbar
        gradients = np.linspace (self.colorbar_ranges[0],
                                 self.colorbar_ranges[1], 1000)
        norm = matplotlib.colors.BoundaryNorm (gradients, self.cmap.N)

        ### define and plot colorbar !
        ticks = np.linspace (self.colorbar_ranges[0],
                             self.colorbar_ranges[1],
                             self.colorbar_nbins     )
        cb = matplotlib.colorbar.ColorbarBase (axis,
                                               norm=norm,
                                               ticks=ticks,
                                               cmap=self.cmap, 
                                               alpha=self.colorbar_alpha,
                                               spacing='uniform',
                                               orientation='vertical',
                                               format=self.colorbar_format)

        cb.ax.tick_params (labelsize=self.tick_fontsize,
                           length=self.colorbar_tick_length)
        cb.set_label (self.colorbar_label,
                      fontsize=self.label_fontsize)
        
