#!/usr/bin/env python

####
#### By Elim Thompson (09/23/2018)
####
#### This script contains a Contour class
#### which plot a 2D contour for a standard
#### atmospheric oscillation parameter
#### measurement.
####
################################################

#### import packages
from canvases import Canvas
import numpy as np

#### import plot settings
from defaults import colors, labels, linestyles

################################################
#### Contour inherited from Canvas
################################################
class Contour (Canvas):

    def __init__ (self, **kwargs):

        ### default 1D distribution plots (python 2.7)
        super (Contour, self).__init__ ((11, 7.5), (2, 2),
                                        height_ratios=[1, 3],
                                        width_ratios=[4, 1],
                                        wspace=0, hspace=0,
                                        top=0.97, bottom=0.08,
                                        right=0.98, left=0.08)

        ### set any input style
        self.set_style (**kwargs)

        ### set subplot axes
        ## upper left  = dchi2 vs sin2 theta23
        self.axtheta   = self.h.add_subplot (self.gs[0])
        ## upper right = empty
        ## lower left  = contour
        self.axcontour = self.h.add_subplot (self.gs[2])
        ## lower right = dm23 vs dchi2
        self.axdm      = self.h.add_subplot (self.gs[3])

    def _clean_up (self):

        self.axtheta.cla ()
        self.axcontour.cla ()
        self.axdm.cla ()
        return
    
    def _format (self, axis, xticks, yticks, xranges, yranges):
        
        axis.set_ylim (yranges)
        axis.set_xlim (xranges)

        axis.set_xticks (xticks)
        axis.set_yticks (yticks)
        
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

    def _format_contour (self):

        axis = self.axcontour

        xticks  = self.theta_ticks
        yticks  = self.dm_ticks
        xranges = self.theta_ranges
        yranges = self.dm_ranges
        self._format (axis,
                      xticks , yticks ,
                      xranges, yranges)

        ## set x/y labels
        xlabel  = labels['sin2theta23']
        ylabel  = labels['dm32']
        axis.set_xlabel (xlabel, fontsize=self.label_fontsize)
        axis.set_ylabel (ylabel, fontsize=self.label_fontsize)
        
        ## set legend
        axis.legend (loc=self.legend_loc,
                     ncol=self.legend_ncols,
                     framealpha=self.legend_alpha,
                     prop={'size':self.legend_fontsize})
        return

    def _format_dm (self):

        axis = self.axdm

        xticks  = self.chi2_ticks
        yticks  = self.dm_ticks
        xranges = self.chi2_ranges
        yranges = self.dm_ranges
        self._format (axis,
                      xticks , yticks ,
                      xranges, yranges)

        ## set only x label
        xlabel  = labels['dchi2']
        axis.set_xlabel (xlabel, fontsize=self.label_fontsize)
        
        ## hide y axis tick labels
        axis.get_yaxis ().set_ticklabels ([])
        ## hide the first tick on x axis
        xticks = axis.xaxis.get_major_ticks ()
        xticks[0].set_visible (False)
        ## put x labdl on top axis
        axis.xaxis.set_label_coords (0.5, 1.13)
        ## put x ticks on top axis
        axis.xaxis.set_ticks_position ("top")
        return

    def _format_theta (self):

        axis = self.axtheta

        xticks  = self.theta_ticks
        yticks  = self.chi2_ticks
        xranges = self.theta_ranges
        yranges = self.chi2_ranges
        self._format (axis,
                      xticks , yticks ,
                      xranges, yranges)

        ## set only y label
        ylabel  = labels['dchi2']
        axis.set_ylabel (ylabel, fontsize=self.label_fontsize)
        
        ## hide x axis tick labels
        axis.get_xaxis ().set_ticklabels ([])
        ## hide the first tick on y axis
        yticks = axis.yaxis.get_major_ticks ()
        yticks[0].set_visible (False)
        return
    
    def _plot (self, values, ptype=None):

        try:
            ### define axis
            ## lower left  = contour
            ## upper left  = chi2 vs theta
            ## lower right = dm vs chi2
            axis = eval ('self.ax' + ptype)
        except:
            ### if ptype is not 'contour' or
            ### 'dm' or 'theta', nothing is plotted
            return None

        ### linestyles / colors / labels are
        ### defined in `plotter/defaults.py`
        for name, measurement in values.items ():
            axis.plot (measurement['x'],
                       measurement['y'],
                       color=colors[name],
                       label=labels[name],
                       linestyle=linestyles[name],
                       alpha=self.alpha,
                       linewidth=self.linewidth)
        return

    def _plot_bestfits (self, bestfits):

        axis = self.axcontour

        for name, measurement in bestfits.items ():
            axis.scatter (measurement['x'],
                          measurement['y'],
                          color=colors[name],
                          facecolors=colors[name],
                          edgecolors=colors[name],
                          marker=self.marker,
                          s=self.marker_size,
                          alpha=self.marker_alpha,
                          linewidths=self.marker_linewidth)
        return
    
    def plot (self, contours, bestfits, dm_scans, theta_scans):

        ### start with cleaned subplots
        self._clean_up ()

        ### plot and format contour
        self._plot (contours, ptype='contour')
        self._plot_bestfits (bestfits)
        self._format_contour ()

        ## plot and format dm31 scan
        #  NOTE: y values must be dm
        #        x values must be chi2 
        self._plot (dm_scans, ptype='dm')
        self._format_dm ()
        
        ## plot and format theta23 scan
        #  NOTE: y values must be chi2
        #        x values must be sin2theta23 
        self._plot (theta_scans, ptype='theta')
        self._format_theta ()        
        
        return 


