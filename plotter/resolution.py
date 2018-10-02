#!/usr/bin/env python

####
#### By Elim Thompson (09/23/2018)
####
#### This script contains a Resolutiion class
#### which plot the energy / zenith resolution
#### (i.e. median / 1 sigma deviations) for
#### both dragon and greco samples and for
#### 3 CC and 1 merged NC interaction types.
####
################################################

#### import packages
from canvases import Canvas
import numpy as np

#### import plot settings
from styler import AnchoredText
from defaults import colors, labels

################################################
#### Resolution inherited from Canvas
################################################
class Resolution (Canvas):

    def __init__ (self, **kwargs):

        ### default 1D distribution plots (python 2.7)
        super (Resolution, self).__init__ ((14, 11), (2, 2),
                                           height_ratios=[1, 1],
                                           width_ratios=[1, 1],
                                           wspace=0.1, hspace=0.15,
                                           top=0.98, bottom=0.10,
                                           right=0.98, left=0.11)

        ### xaxis and yaxis as properties
        self._xaxis = None
        self._yaxis = None
        
        ### set any input style
        self.set_style (**kwargs)

        ### set subplot axes
        ## upper left  = numucc 
        self.axnumucc  = self.h.add_subplot (self.gs[0])
        ## upper right = nuecc
        self.axnuecc   = self.h.add_subplot (self.gs[1])
        ## lower left  = nutaucc
        self.axnutaucc = self.h.add_subplot (self.gs[2])
        ## lower right = nunc
        self.axnunc    = self.h.add_subplot (self.gs[3])

    @property
    def xaxis (self):
        return self._xaxis

    @xaxis.setter
    def xaxis (self, xaxis):
        self._xaxis = xaxis

    @property
    def yaxis (self):
        return self._yaxis

    @yaxis.setter
    def yaxis (self, yaxis):
        self._yaxis = yaxis
        
    def _clean_up (self):

        self.axnumucc.cla ()
        self.axnuecc.cla ()
        self.axnutaucc.cla ()
        self.axnunc.cla ()
        return
    
    def _format_ticks (self, axis, xticks, yticks, xranges,
                       yranges, xtickslabel, ytickslabel):
        
        axis.set_ylim (yranges)
        axis.set_xlim (xranges)

        axis.set_xticks (xticks)
        axis.set_yticks (yticks)

        axis.set_xticklabels (xtickslabel)
        axis.set_yticklabels (ytickslabel)
        
        axis.tick_params (axis='x', labelsize=self.tick_fontsize)
        axis.tick_params (axis='y', labelsize=self.tick_fontsize)
        
        for xtick in xticks:
            axis.axvline (x=xtick, color='gray',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)

        for ytick in yticks:
            axis.axhline (y=ytick, color='gray',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)
        
        return

    def _format_legend (self, axis, member):

        ### exit if plot legend in this
        ### member resolution plot
        legend_plot_key = 'self.legend_'+self.xaxis+'_'+self.yaxis+'_plot'
        if not eval (legend_plot_key) == member: return

        ### define location and ncolumns
        loc   = eval ('self.legend_'+self.xaxis+'_'+self.yaxis+'_loc')
        ncols = eval ('self.legend_'+self.xaxis+'_'+self.yaxis+'_ncols')
        
        axis.legend (loc=loc,
                     ncol=ncols,
                     framealpha=self.legend_alpha,
                     prop={'size':self.legend_fontsize})
        return
    
    def _format_label (self, axis, member, xlabel, ylabel):
        
        ### set only y label if numucc or nutaucc
        if member in ['numucc', 'nutaucc']:
            axis.set_ylabel (ylabel, fontsize=self.label_fontsize)
        ### set only x label if nutaucc or nunc
        if member in ['nutaucc', 'nunc']:
            axis.set_xlabel (xlabel, fontsize=self.label_fontsize)
        
        ### tick labels
        ### hide x axis tick labels if numucc or nuecc
        if member in ['numucc', 'nuecc']:
            axis.get_xaxis ().set_ticklabels ([])
        ### hide y axis tick labels if nuecc or nunc
        if member in ['nuecc', 'nunc']:
            axis.get_yaxis ().set_ticklabels ([])        
        return

    def _format_text (self, axis, member):

        ### define text location
        loc_key = 'self.text_'+self.xaxis+'_'+self.yaxis+'_'+member+'_loc'
        loc = eval (loc_key)

        ### print data type
        at = AnchoredText (labels[member],
                           prop=dict (size=self.text_fontsize),
                           frameon=self.text_frameon,
                           loc=loc)
        at.patch.set_boxstyle ("round,pad=0.,rounding_size=0.5")
        axis.add_artist(at)
        return
            
    def _format (self, member):

        try:
            ### define axis
            ## upper left  = numucc
            ## upper right = nuecc
            ## lower left  = nutaucc
            ## lower right = nunc
            axis = eval ('self.ax' + member)
        except:
            ### if ptype is not 'numucc' or 'nuecc' or
            ### 'nutaucc' or 'nunc', nothing is formatted
            return None

        ### set ticks
        xticks  = eval ('self.' + self.xaxis + '_xticks')
        yticks  = eval ('self.' + self.yaxis + '_yticks')
        xranges = eval ('self.' + self.xaxis + '_xranges')
        yranges = eval ('self.' + self.yaxis + '_yranges')
        xtickslabel = eval ('self.' + self.xaxis + '_xtickslabel')
        ytickslabel = eval ('self.' + self.yaxis + '_ytickslabel')
        self._format_ticks (axis,
                            xticks , yticks ,
                            xranges, yranges,
                            xtickslabel, ytickslabel)

        ### set x / y labels
        xlabel  = eval ('self.' + self.xaxis + '_xlabel')
        ylabel  = eval ('self.' + self.yaxis + '_ylabel')
        self._format_label (axis, member, xlabel, ylabel)

        ### set legend
        self._format_legend (axis, member)

        ### set text
        self._format_text (axis, member)
        return

    def _plot_ref_line (self, axis, xvalues):

        ### does not plot reference line
        ### if not asked to 
        if not self.ref_plot: return

        ### define yvalues
        yvalues = xvalues if self.yaxis in ['lineare', 'log10e'] else \
                  np.zeros (len (xvalues))

        ### plot the reference line
        axis.plot (xvalues, yvalues,
                   color=self.ref_color,
                   alpha=self.ref_alpha,
                   linestyle=self.ref_linestyle,
                   linewidth=self.ref_linewidth)
        return
    
    def _plot (self, sample, member, res):

        try:
            ### define axis
            ## upper left  = numucc
            ## upper right = nuecc
            ## lower left  = nutaucc
            ## lower right = nunc
            axis = eval ('self.ax' + member)
        except:
            ### if ptype is not 'numucc' or 'nuecc' or
            ### 'nutaucc' or 'nunc', nothing is plotted
            return None

        ### define global color / label 
        color = eval ('self.'+sample+'_color') 
        label = eval ('self.'+sample+'_label') 

        ### plot median and sigma lines
        for sigma in ['lsigma', 'median', 'usigma']:
            ## define label and linestyle
            lab = label + r' 50\%' if sigma == 'median' else \
                  label + r' 68\%' if sigma == 'lsigma' else None
            linestyle = self.median_linestyle if sigma == 'median' else \
                        self.sigma_linestyle
            ## plot line
            axis.plot (res['xvalues'], res[sigma],
                       color=color,
                       label=lab,
                       alpha=self.alpha,
                       linestyle=linestyle,
                       linewidth=self.linewidth)

        ### plot reference line
        self._plot_ref_line (axis, res['xvalues'])
        return

    def plot (self, resolutions):

        ### start with cleaned subplots
        self._clean_up ()

        ### loop though each data type
        for member, resolution in resolutions.items ():
            ## loop through each sample
            for sample, res in resolution.items ():
                ## plot
                self._plot (sample, member, res)
                ## format
                self._format (member)

        return 


