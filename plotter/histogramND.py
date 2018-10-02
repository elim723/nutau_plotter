#!/usr/bin/env python

####
#### By Elim Thompson (10/01/2018)
####
#### This script contains a histogramND class
#### which plots systematic effect plots.
####
################################################

#### import packages
from canvases import Canvas2D
from copy import deepcopy
from styler import plt
import numpy as np
import sys

#### import plot settings
from defaults import colors, hatches, labels

################################################
#### histogramND inherited from Canvas2D
################################################
class HistogramND (Canvas2D):

    def __init__ (self, sample, parameters,
                  stretch=False, **kwargs):

        ### define subplot settings
        nrows = len (parameters) 
        height = 3*len (parameters)
        width_ratios = [30, 22, 1] if stretch else \
                       [30, 30, 1]

        ### default ND distribution plots (python 2.7)
        super (HistogramND, self).__init__ ((7.5, height), (nrows, 3),
                                            width_ratios=width_ratios,
                                            wspace=0.15, hspace=0.2,
                                            top=0.97, bottom=0.04,
                                            right=0.9, left=0.11)
        
        ### set any input style
        self.set_style (**kwargs)

        ### define observable
        self._xobs = None
        self._yobs = None
        self._zobs = None
        self._sample = sample.lower ()
        self._parameters = sorted (parameters)

        ### define axes
        self._set_axes (stretch)

    def _set_axes (self, stretch):

        ### set subplot axes
        for index, param in enumerate (self.parameters):

            ## left = cascade-like histogram
            setattr (self, 'axcascade_'+param,
                     self.h.add_subplot (self.gs[3*index+0]))

            ## middle = track-like histogram
            if self._sample=='greco' and not stretch:
                width = (0.915-0.1)*30/(30+30+1.)
                height = 0.91-0.20
                left = width + (width - width*0.67) + 0.063
                bottom = 0.20
                axtrack = self.h.add_axes ([left, bottom, width*0.67, height])
            else:
                axtrack = self.h.add_subplot (self.gs[3*index+1])
            setattr (self, 'axtrack_'+param, axtrack)

            ## right = colorbar
            setattr (self, 'axcolorbar_'+param,
                     self.h.add_subplot (self.gs[3*index+2]))

        return

    @property
    def sample (self):
        return self._sample
        
    @sample.setter
    def sample (self, sample):
        self._sample = sample

    @property
    def parameters (self):
        return self._parameters
        
    @parameters.setter
    def parameters (self, parameters):
        self._parameters = parameters

    @property
    def xobs (self):
        return self._xobs
        
    @xobs.setter
    def xobs (self, xobs):
        self._xobs = xobs

    @property
    def yobs (self):
        return self._yobs
        
    @yobs.setter
    def yobs (self, yobs):
        self._yobs = yobs

    @property
    def zobs (self):
        return self._zobs
        
    @zobs.setter
    def zobs (self, zobs):
        self._zobs = zobs

    def _clean_up (self):

        ### set subplot axes
        for param in self.parameters:

            ## left = cascade-like histogram
            axis = getattr (self, 'axcascade_'+param)
            axis.cla ()

            ## middle = track-like histogram
            axis = getattr (self, 'axtrack_'+param)
            axis.cla ()

            ## right = colorbar
            axis = getattr (self, 'axcolorbar_'+param)
            axis.cla ()

        return
    
    def _format (self, axis, xticks, yticks,
                 xlabel, ylabel, right_title,
                 left_title, xticklabels, yticklabels):

        ### set format
        self._format2D (axis, xticks, yticks, xlabel, ylabel,
                        xticklabels=xticklabels,
                        yticklabels=yticklabels )

        ### set title
        if right_title: axis.set_title (right_title, loc='right',
                                        fontsize=self.title_fontsize)
        if left_title: axis.set_title (left_title, loc='left',
                                       fontsize=self.title_fontsize)
        return

    def _format_cascade (self, axis, param, index):

        is_last_cascade = index == len (self.parameters)-1

        ### set x/y label and title
        ylabel = labels['reco_cz'] 
        xlabel = labels['reco_e'] if is_last_cascade else None
        right_title = r'Cascade-like' if index==0 else None
        left_title  = labels[param]

        ### x/y ticklabels
        xticklabels = self.xticklabels if is_last_cascade else []

        ### format
        self._format (axis, self.xticks, self.yticks,
                      xlabel, ylabel, right_title,
                      left_title, xticklabels, self.yticklabels )
        return

    def _format_track (self, axis, index):

        is_last_track = index == len (self.parameters)-1

        ### set x/y label and title
        ylabel = None
        xlabel = labels['reco_e'] if is_last_track else None
        right_title  = r'Track-like' if index==0 else None
        left_title = None

        ### set initial ticks
        xticks = self.xticks
        xticklabels = getattr (self, 'xticklabels', xticks)
        ##  for analysis histogram, ticks = edges
        ##  if greco and track, skip the first two tick
        if self.sample.lower ()=='greco':
            xticks = xticks[1:]
            xticklabels = xticklabels[1:]

        ### reset ticklabels
        yticklabels = []
        xticklabels = xticklabels if is_last_track else []

        ### format
        self._format (axis, xticks, self.yticks,
                      xlabel, ylabel, right_title,
                      left_title, xticklabels, yticklabels)
        return

    def _set_colorbar (self, param):

        ''' plot colorbar (local function instead of
            the one from Canvas2D

            :param axis (`matplotlib.axes`): the axis things to be plotted
        '''

        self.colorbar_ranges = getattr (self, 'colorbar_ranges_' + param)
        self.colorbar_nbins  = getattr (self, 'colorbar_nbins_'  + param)
        self.colorbar_alpha  = getattr (self, 'colorbar_alpha_'  + param)
        self.colorbar_format = getattr (self, 'colorbar_format_' + param)
        self.colorbar_label  = r'Percentage change (\%)'
        return
        
    def _get_content (self, base, changed):

        ''' determine content to be plotter

            syseffects = % change in counts per bin compared
                         to seeded (baseline) weight
        '''

        total_base = sum ([ base[dtype] for dtype in base.keys () ])
        total_changed = sum ([ changed[dtype] for dtype in changed.keys () ])
        if self.norm:
            total_changed *= np.sum (total_base) / np.sum (total_changed)

        diff = total_changed - total_base
        return np.divide (diff, total_base) * 100.

    def _get_a_histogram (self, event, weightkey):

        ''' get a 2D histogram for a data type. '''

        ### define values and weights
        xvalues = event[self.xobs]
        yvalues = event[self.yobs]
        zvalues = event[self.zobs]
        weights = event[weightkey] 
        
        ### make sure they are finite
        finite = np.logical_and (np.isfinite (xvalues), np.isfinite (yvalues))
        finite = np.logical_and (np.isfinite (zvalues), finite)
        finite = np.logical_and (np.isfinite (weights), finite)

        ### define 3D histogram edges / ranges
        edges  = (self.xedges, self.yedges, self.zedges)
        ranges = [ [self.xedges[0], self.xedges[-1]],
                   [self.yedges[0], self.yedges[-1]],
                   [self.zedges[0], self.zedges[-1]] ]

        ### build a 3D histogram
        data = [xvalues[finite], yvalues[finite], zvalues[finite]]
        H    = np.histogramdd (np.array(data).T,
                               edges           ,
                               range=ranges    ,
                               weights=weights[finite])[0]
        ### return 3D histogram
        return H

    def _get_histograms (self, events):

        ### define how many histograms to
        ### be made (seeded + n_parameters)
        weightkeys = [ key for key in
                       events[events.keys ()[0]].keys ()
                       if 'weight' in key ]

        ### initialize histogram container
        histograms = {}

        ### loop through each weight key
        ### and build a histogram for it
        for weightkey in weightkeys:

            ## define param name for this weightkey
            param = '_'.join (weightkey.split ('_')[:-1])

            ## build a 3D histogram for each data type
            histograms[param] = {'cascade':{},
                                 'track'  :{}}
            
            ## loop through each member
            for member, event in events.items ():
                # build a 3D histogram for this member
                histos = self._get_a_histogram (event, weightkey)
                # define cascade / track histograms
                cascade = histos[:,:,0]
                track   = histos[:,:,1]
                # ignore first two track bins if greco
                if self.sample.lower () == 'greco':
                    track = track[1:]
                # store histograms
                histograms[param]['cascade'][member] = cascade
                histograms[param]['track'][member]   = track
        
        return histograms

    def _plot (self, axis, base, changed, ptype=None):

        ### return None if no valid ptype
        if not ptype: return None

        ### get content
        content = self._get_content (base[ptype], changed[ptype])
        content = np.nan_to_num (content)

        ### set edges
        xticks = deepcopy (self.xedges)
        if self.sample.lower ()=='greco' and \
           ptype=='track': xticks = xticks[1:]

        ### plot content
        self._plot2D (axis, content, xedges=xticks)

        return 

    def plot (self, events):

        self._clean_up ()

        ### set edges
        self._set_edges ()

        ### build all 3D histograms
        histograms = self._get_histograms (events)

        ### loop through each parameter
        for index, param in enumerate (self.parameters):

            ### set colorbar formatting for this param
            self._set_colorbar (param)

            ### define histograms needed
            ### for this param
            base = histograms['seeded']
            changed = histograms[param]
            
            ### left (self.axcascade) - cascade distribution
            axcascade = getattr (self, 'axcascade_'+param)
            self._plot (axcascade, base, changed, ptype='cascade')
            self._format_cascade (axcascade, param, index)
            
            ### middle (self.axtrack) - track distribution
            axtrack = getattr (self, 'axtrack_'+param)
            self._plot (axtrack, base, changed, ptype='track')
            self._format_track (axtrack, index)
            
            ### right (self.axcolorbar) - colorbar
            axcolorbar = getattr (self, 'axcolorbar_'+param)
            self._plot_colorbar (axcolorbar)
        
        return 


