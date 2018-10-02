#!/usr/bin/env python

####
#### By Elim Thompson (09/26/2018)
####
#### This script contains a histogram2D class
#### which deals with all 2D distribution plots.
####
################################################

#### import packages
from canvases import Canvas2D
from styler import plt
import numpy as np
import sys

#### import plot settings
from defaults import colors, hatches, labels

################################################
#### histogram2D inherited from Canvas2D
################################################
class Histogram2D (Canvas2D):

    def __init__ (self, **kwargs):

        ### default 2D distribution plots (python 2.7)
        super (Histogram2D, self).__init__ ((7.5, 5.5), (1, 2), 
                                            width_ratios=[28, 1],
                                            wspace=0.07, hspace=0,
                                            top=0.94, bottom=0.12,
                                            right=0.9, left=0.13)
        
        ### set any input style
        self.set_style (**kwargs)

        ### set subplot axes
        ## left = 2D histogram
        self.axhist     = self.h.add_subplot (self.gs[0])
        ## right = colorbar
        self.axcolorbar = self.h.add_subplot (self.gs[1])

    def _clean_up (self):

        self.axhist.cla ()
        self.axcolorbar.cla ()
        return

    def _format_hist (self, xvar, yvar):

        axis = self.axhist
        ### set ticks
        xticks = self.hist_xticks if hasattr (self, 'hist_xticks') else \
                 axis.xaxis.get_majorticklocs ()
        yticks = self.hist_yticks if hasattr (self, 'hist_yticks') else \
                 axis.yaxis.get_majorticklocs ()
        ### set label
        xlabel = self.hist_xlabel if hasattr (self, 'hist_xlabel') else \
                 labels[xvar]
        ylabel = self.hist_ylabel if hasattr (self, 'hist_ylabel') else \
                 labels[yvar]
        ### set format
        self._format2D (axis, xticks, yticks, xlabel, ylabel)
        return
        
    def _pull_member (self, histos, member):

        ''' get the muon distribution '''

        ### included members
        members = np.array (histos.keys ())

        ### just return it if member is already in members
        if member in histos:
            return histos[member]
            
        ### if it is muon related, look for a possible histogram
        if member == 'muon':
            # possible muon key names
            mutypes = ['muon', 'muongun', 'corsika',
                       'icc', 'icc_data', 'iccdata']
            haskey  = np.in1d (members, mutypes)
            if haskey.any ():
                keyname = members[haskey][0]
                return histos[keyname]

        ### if not the above, print warning and return None
        print ('You are trying to get {0}, which'.format (member))
        print ('is not in your histogram dictionary ...')
        print ('You will see an empty plot.')
        return None

    def _get_content (self, histos):

        ''' determine content to be plotter
        
            4 options
            ---------
               1. numucc = distribution of numu CC (same for any member)
               2. nurate = total neutrino rates
                           (default to show micro Hz)
               3. mufrac = muon fraction
               4. s/b    = signal / sqrt (background)
        '''

        ### option 1: histogram of a specific member
        if self.htype in histos.keys ():
            return histos[self.htype]

        ### option 2: summed neutrino histograms
        if self.htype == 'nurate':
            ## only add the neutrino in dtype 
            print ([ dtype for dtype in histos.keys () ])
            return sum ([ histo for dtype, histo in histos.items ()
                          if 'nu' in dtype ])

        ### option 3: muon fraction
        if self.htype == 'mufrac':
            totalmc  = sum ([ histo for dtype, histo in histos.items () ])
            muon     = self._pull_member (histos, 'muon')
            return np.divide (muon, totalmc)

        ### option 4: signal / sqrt (background)
        if self.htype == 's/b':
            ## by default, signal     = nutau (cc + nc)
            ##             background = everything else
            background = sum ([ histo for dtype, histo in histos.items ()
                                if not 'nutau' in dtype ])
            signal     = sum ([ histo for dtype, histo in histos.items ()
                                if 'nutau' in dtype ])
            return np.divide (signal, np.sqrt (background))

        ### if not the above, print warning and return None
        print ('You are trying to get {0}, which'.format (self.htype))
        print ('is not implemented ...')
        print ('You will see an empty plot.')
        return None

    def _get_a_histogram (self, event, xvar, yvar):

        ''' get a 2D histogram for a data type. '''

        ### define values and weights
        xvalues = event[xvar]
        yvalues = event[yvar]
        weights = event['weight']

        ### make sure they are finite
        finite = np.logical_and (np.isfinite (xvalues), np.isfinite (yvalues))
        finite = np.logical_and (np.isfinite (weights), finite)

        ### define 2D histogram edges / ranges
        edges  = (self.xedges, self.yedges)
        ranges = [ [self.xedges[0], self.xedges[-1]],
                   [self.yedges[0], self.yedges[-1]] ]

        ### build a 2D histogram
        data = [xvalues[finite], yvalues[finite]]
        H    = np.histogramdd (np.array(data).T,
                               edges           ,
                               range=ranges    ,
                               weights=weights[finite])[0]
        ### return 2D histogram
        return H

    def _plot_histograms (self, events, xvar, yvar):

        ### build a 2D histogram for each data type
        histograms = {}
        for member, event in events.items ():
            histograms[member] = self._get_a_histogram (event, xvar, yvar)

        ### get content
        content = np.nan_to_num (self._get_content (histograms))
        
        ### plot content
        self._plot2D (self.axhist, content)

        return 

    def _plot_cut_vline (self):

        ### nothing to be plotted
        ### if cut_vline is not defined
        if not hasattr (self, 'cut_vline') or \
           not self.cut_vline: return

        ### must be hist axis
        axis = self.axhist

        ### loop through each available cut lines
        for cut in self.cut_vline:

            ## define values
            xvalue = cut['xvalue']
            ymin   = cut['yranges'][0]
            ymax   = cut['yranges'][1]

            ## plot vertical line
            axis.axvline (x=xvalue,
                          ymin=ymin,
                          ymax=ymax,
                          alpha=self.cut_alpha,
                          color=self.cut_color,
                          linestyle=self.cut_linestyle,
                          linewidth=self.cut_linewidth)
        return

    def _plot_cut_hline (self):

        ### nothing to be plotted
        ### if cut_hline is not defined
        if not hasattr (self, 'cut_hline') or \
           not self.cut_hline: return

        ### must be hist axis
        axis = self.axhist

        ### loop through each available cut lines
        for cut in self.cut_hline:

            ## define values
            yvalue = cut['yvalue']
            xmin   = cut['xranges'][0]
            xmax   = cut['xranges'][1]

            ## plot vertical line
            axis.axhline (y=yvalue,
                          xmin=xmin,
                          xmax=xmax,
                          alpha=self.cut_alpha,
                          color=self.cut_color,
                          linestyle=self.cut_linestyle,
                          linewidth=self.cut_linewidth)
        return

    def _plot_cut_line (self):

        ### nothing to be plotted
        ### if cut_line is not defined
        if not hasattr (self, 'cut_line') or \
           not self.cut_line: return

        ### must be hist axis
        axis = self.axhist

        ### loop through each available cut lines
        for cut in self.cut_line:

            ## define values
            xvalues = cut['xvalues']
            yvalues = cut['yvalues']

            ## plot line
            axis.plot (xvalues, yvalues,
                       alpha=self.cut_alpha,
                       color=self.cut_color,
                       linestyle=self.cut_linestyle,
                       linewidth=self.cut_linewidth)
        return
        
    def plot (self, events, xvariable, yvariable):

        self._clean_up ()

        ### set edges
        self._set_edges ()
        
        ### left (self.axhist) - distribution plot
        self._plot_histograms (events, xvariable, yvariable)
        ## plot the cut lines
        self._plot_cut_vline ()
        self._plot_cut_hline ()
        self._plot_cut_line  ()
        ## format self.axhist
        self._format_hist (xvariable, yvariable)

        ### right (self.axcolorbar) - colorbar
        self._plot_colorbar (self.axcolorbar)
        
        return 


