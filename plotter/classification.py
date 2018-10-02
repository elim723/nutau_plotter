#!/usr/bin/env python

####
#### By Elim Thompson (09/19/2018)
####
#### This script contains a Classification
#### class which plot the fraction of
#### track-like events per MC truth energy
#### slice.
####
#### By default, left plot is greco,
#### right plot is dragon. Only neutrinos
#### are included separated by interaction
#### types.
################################################

#### import packages
from canvases import Canvas
from styler import plt, AnchoredText
import numpy as np

#### import plot settings
from defaults import colors, hatches, labels

################################################
#### Classification inherited from Canvas
################################################
class Classification (Canvas):

    def __init__ (self, samples=[],
                  loge_edges=None,
                  **kwargs):

        ### default 1D distribution plots (python 2.7)
        super (Classification, self).__init__ ((17, 6), (1, 2), 
                                               width_ratios=[1, 1],
                                               wspace=0.15, hspace=0.4,
                                               top=0.95, bottom=0.17,
                                               right=0.98, left=0.1)

        self._samples = samples
        self._loge_edges = loge_edges
        ### set any input style
        self.set_style (**kwargs)

        ### set subplot axes
        ## left = track fraction from greco (analysis A)
        self.axgreco  = self.h.add_subplot (self.gs[0])
        ## right = track fraction from dragon (analysis B)
        self.axdragon = self.h.add_subplot (self.gs[1])

    @property
    def samples (self):
        return self._samples

    @samples.setter
    def samples (self, samples):
        self._samples = samples

    @property
    def loge_edges (self):
        return self._loge_edges

    @loge_edges.setter
    def loge_edges (self, loge):
        self._loge_edges = loge

    def _clean_up (self):

        self.axgreco.cla ()
        self.axdragon.cla ()
        return
    
    def _format (self, sample):

        ### define axis (left = greco, right = dragon)
        axis = self.axgreco if sample == 'greco' else \
               self.axdragon

        ### set x axis
        axis.set_xlim (self.xranges)
        axis.set_xticks (self.xticks)
        axis.set_xticklabels (self.xtickslabel)
        axis.tick_params (axis='x', labelsize=self.tick_fontsize)
        axis.set_xlabel (self.xlabel, fontsize=self.label_fontsize)

        ### draw grids
        for ytick in self.yticks:
            axis.axhline (y=ytick, color='gray',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)
        for xtick in self.xticks:
            axis.axvline (x=xtick, color='gray',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)

        ## set y label / tick params
        axis.set_yticks (self.yticks)
        axis.set_ylim (self.yranges)
            
        if sample == 'greco':
            axis.set_ylabel (self.ylabel,
                             multialignment='center',
                             fontsize=self.label_fontsize)
            axis.set_yticklabels (self.ytickslabel)
            axis.tick_params (axis='y', labelsize=self.tick_fontsize)
            ## set legend
            axis.legend (loc=self.legend_loc,
                         ncol=self.legend_ncols,
                         framealpha=self.legend_alpha,
                         prop={'size':self.legend_fontsize})

        if sample == 'dragon':
            ## hide y axis if right dragon plot
            axis.get_yaxis().set_ticks ([])
        return
        
    def _plot_a_track_frac (self, sample, track_frac):

        ### define axis (left = greco / right = dragon)
        axis = self.axgreco if sample == 'greco' else \
               self.axdragon
        ### xvalues = bin center of loge edges
        xvalues = self._get_bin_center (self.loge_edges)

        ### loop through each member / data type 
        for member, trackf in track_frac.items ():
            yvalues = track_frac[member]
            color   = colors[member]
            label   = labels[member]
            axis.plot (xvalues, yvalues,
                       color=color,
                       label=label,
                       alpha=self.alpha,
                       linewidth=self.linewidth,
                       linestyle=self.linestyle)

        ### write analysis A / B
        text_loc = self.text_greco_loc if sample == 'greco' else \
                   self.text_dragon_loc
        at = AnchoredText (labels[sample],
                           prop=dict (size=self.text_fontsize),
                           frameon=self.text_frameon,
                           loc=text_loc)
        at.patch.set_boxstyle ("round,pad=0.,rounding_size=0.5")
        axis.add_artist(at)
        return
    
    def plot (self, track_fracs):

        ### start with cleaned subplots
        self._clean_up ()

        ### loop through samples 
        for sample in self.samples:
            ## plot track fractions of a sample
            self._plot_a_track_frac (sample, track_fracs[sample])
            ## format the subplot
            self._format (sample)
        
        return 


