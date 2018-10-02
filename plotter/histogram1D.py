#!/usr/bin/env python

####
#### By Elim Thompson (09/19/2018)
####
#### This script contains a histogram1D class
#### which deals with all 1D distribution plots
#### given the variables/weights/plot settings.
####
#### By default, distribution is on the upper
#### plot, and MC/data ratio is on the lower
#### plot.
################################################

#### import packages
from canvases import Canvas
from styler import plt
import numpy as np
import sys

#### import toolbox
from analyzer.misc import Toolbox
toolbox = Toolbox ()

#### import plot settings
from defaults import colors, hatches, labels

errorbar_kwargs = ['alpha', 'capsize', 'capthick',
                   'elinewidth', 'linestyle',
                   'linewidth', 'marker',
                   'markersize', 'fillstyle',
                   'markeredgewidth']

################################################
#### histogram1D inherited from Canvas
################################################
class Histogram1D (Canvas):

    def __init__ (self, **kwargs):

        ### default 1D distribution plots (python 2.7)
        super (Histogram1D, self).__init__ ((7.5, 5.5), (2, 1), 
                                            height_ratios=[3.5, 1],
                                            wspace=0, hspace=0.1,
                                            top=0.96, bottom=0.12,
                                            right=0.97, left=0.1)

        self._scalar = 1.0
        ### set any input style
        self.set_style (**kwargs)

        ### set subplot axes
        ## upper = histogram
        self.axhist  = self.h.add_subplot (self.gs[0])
        ## lower = ratio
        self.axratio = self.h.add_subplot (self.gs[1])

    @property
    def scalar (self):
        return self._scalar

    @scalar.setter
    def scalar (self, scalar):
        self._scalar = scalar

    @staticmethod
    def _order_stack (events):

        ''' order the MC member by rates from lowest
            to highest so that the contribution from
            the lowest shows up in a logy plot.

            :param events (dict): all event variables
        
            :return members (list): ordered members
        '''
        
        rates = {member:np.sum (event['weight'])
                 for member, event in events.items ()
                 if not 'data' in member}
        rates = sorted (rates.items (),
                        key = lambda t: t[1])
        ## rates here is a list of tuple (member, rate)
        return [ tup[0] for tup in rates ]

    def _scale_to_counts (self, array):

        return np.array (array) * self._scalar

    def _scale_to_hertz (self, array):

        return np.array (array) / self._scalar
    
    def _clean_up (self):

        self.axhist.cla ()
        self.axratio.cla ()
        return
    
    def _set_errkwargs (self, plottype):

        ### modify errorbar_kwargs above to
        ### add new plotting format parameters
        ### related to data points / error bars
        
        kwargs = {}
        for errkwarg in errorbar_kwargs:
            kwargs[errkwarg] = eval ('self.'+plottype+'_'+errkwarg) \
                               if hasattr (self, plottype+'_'+errkwarg) else \
                                  None
        return kwargs

    def _format_plot (self, axis, xticks, yticks, ylabel):

        ### set y axis
        axis.set_ylabel (ylabel, fontsize=self.label_fontsize)
        axis.tick_params (axis='y', labelsize=self.tick_fontsize)

        ### set ticks
        axis.set_yticks (yticks)
        axis.set_xticks (xticks)
        
        ### set x axis
        axis.set_xlim (self.xranges[0], self.xranges[1])

        ### draw grids
        for ytick in yticks:
            axis.axhline (y=ytick, color='gray',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)
        for xtick in xticks:
            axis.axvline (x=xtick, color='gray',
                          alpha=self.grid_alpha,
                          linestyle=self.grid_linestyle,
                          linewidth=self.grid_linewidth)
        return

    def _format_hist (self):

        axis = self.axhist
        ### set log if asked
        if self.hist_logy: axis.set_yscale ("log", nonposy='clip')
        ### set y limits
        axis.set_ylim (self.yranges[0], self.yranges[1])
        ### set ticks
        xticks = self.hist_xticks if hasattr (self, 'hist_xticks') else \
                 axis.xaxis.get_majorticklocs ()
        yticks = self.hist_yticks if hasattr (self, 'hist_yticks') else \
                 axis.yaxis.get_majorticklocs ()
        ### set format
        self._format_plot (axis, xticks, yticks, self.hist_ylabel)
        ### set legend
        axis.legend (loc=self.hloc, ncol=self.ncols,
                     framealpha=self.hist_legend_alpha,
                     prop={'size':self.hist_legend_fontsize})
        ### hide x tick parameters 
        axis.get_xaxis ().set_ticks ([])        
        return
        
    def _format_ratio (self, xLabel):

        axis = self.axratio
        ### set log if asked
        if self.ratio_logy: axis.set_yscale ("log", nonposy='clip')
        ### set y limits
        axis.set_ylim (self.rranges[0], self.rranges[1])
        ### set ticks
        xticks = self.ratio_xticks if hasattr (self, 'ratio_xticks') else \
                 axis.xaxis.get_majorticklocs ()
        yticks = self.ratio_yticks if hasattr (self, 'ratio_yticks') else \
                 axis.yaxis.get_majorticklocs ()
        ### set format
        self._format_plot (axis, xticks, yticks, self.ratio_ylabel)
        ### set legend
        ##  default two ncols for ratio plots
        axis.legend (loc=self.rloc, ncol=2,
                     framealpha=self.ratio_legend_alpha,
                     prop={'size':self.ratio_legend_fontsize})
        ### set x axis
        xlabel = self.xlabel_prefix + ' ' + xLabel
        axis.set_xlabel (xlabel, fontsize=self.label_fontsize)
        axis.tick_params (axis='x', labelsize=self.tick_fontsize)
        ### show reference line if available
        if hasattr (self, 'ratioref_yvalue'):
            axis.axhline (y=self.ratioref_yvalue,
                          alpha=self.ratioref_alpha,
                          color=self.ratioref_color,
                          linestyle=self.ratioref_linestyle,
                          linewidth=self.ratioref_linewidth)
        return
        
    def _get_stack_err (self, stackedHs, xvalues, weights, Colors):

        ### center bin for plotting errorbars
        bin_center = self._get_bin_center (self.edge)
        
        ### sum err from all MC per bin for total MC err
        totalsqrtH = np.zeros (len (bin_center))
        
        ### loop through each row (i.e. sorted member)
        for H, xvalue, weight, color in zip (stackedHs, xvalues, weights, Colors):
            ## build a histogram with weight**2
            count = self._scale_to_counts (weight)
            H2 = np.histogram (xvalue, bins=self.nbins,
                               range=self.xranges,
                               weights=(count)**2)[0]
            H2 = np.nan_to_num (H2)
            yerror = self._scale_to_hertz (np.sqrt (H2))
            totalsqrtH += yerror

            ## set errorbar kwargs
            kwargs = self._set_errkwargs ('stackerr')
            self.axhist.errorbar (bin_center, H, xerr=0,
                                  yerr=yerror,
                                  ecolor=color, color=color,
                                  markeredgecolor=color,
                                  markerfacecolor=color,
                                  **kwargs)
        ### return total MC uncertainties
        return totalsqrtH
        
    def _get_stack_hist (self, events, variable):

        ### holders for stacked histograms
        xvalues, weights = [], []
        Hatches, Labels, Colors = [], [], []
        Hs, H2s, edges = [], [], []

        ### loop through each MC member
        ### sorted from lowest rate to highest
        for member in self._order_stack (events):
            ## append info as a giant array
            xvalues.append (events[member][variable])
            weights.append (events[member]['weight'])
            ## color / label / hatch set in `default.py`
            Hatches.append (hatches[member])
            Colors.append  (colors[member])
            Labels.append  (labels[member])

        ### plot stacked histogram
        stackedHs, edge, patches = self.axhist.hist (xvalues, weights=weights,
                                                     bins =self.nbins,
                                                     range=self.xranges,
                                                     color=Colors, label=Labels,
                                                     histtype='stepfilled',
                                                     stacked=True, fill=True,
                                                     alpha=self.stack_alpha,
                                                     linewidth=self.stack_linewidth,
                                                     linestyle=self.stack_linestyle)
        
        ### add patches
        for patch, hatch in zip(patches, Hatches):
            plt.setp (patch, hatch=hatch)

        ### store edge for the rest of the plots
        self.edge = edge
            
        ### plot stacked errorbars
        totalsqrtH = self._get_stack_err (stackedHs, xvalues, weights, Colors)
        
        ### return total mc
        return {'H':stackedHs[-1], 'sqrtH':totalsqrtH}

    def _plot_mc_histogram (self, events, variable):

        ### plot MC stacked histogram
        totalmc = self._get_stack_hist (events, variable)

        ### plot total MC histo if asked
        if self.plottotalmc:
            ## color / label set in `default.py`
            color, label = colors['mc'], labels['mc']
            bin_center = self._get_bin_center (self.edge)
            self.axhist.plot (bin_center, totalmc['H'],
                              drawstyle='steps-pre',
                              color=color, label=label,
                              alpha=self.mc_alpha,
                              linewidth=self.mc_linewidth,
                              linestyle=self.mc_linestyle)
            ### plot total MC errorbar
            kwargs = self._set_errkwargs ('mcerr')
            self.axhist.errorbar (bin_center, totalmc['H'],
                                  yerr=totalmc['sqrtH'],
                                  ecolor=color, color=color,
                                  markeredgecolor=color,
                                  markerfacecolor=color,
                                  **kwargs)

        ### plot cut line
        self._plot_cut_line (self.axhist)
        return totalmc

    def _plot_data_histogram (self, events, variable, mcrate):

        ### build data histogram
        xvalue = events['data'][variable]
        ## NOTE: data weight is scaled to total mcrate
        weight = np.ones (len (xvalue)).astype (float)
        weight *= mcrate / np.sum (weight)
        dataH = np.histogram (xvalue, bins=self.nbins,
                              range=self.xranges,
                              weights=weight)[0]
        ## scale up to counts to get error bar
        sqrtH = np.sqrt (self._scale_to_counts (dataH))
        sqrtH = self._scale_to_hertz (sqrtH)

        ### plot data points
        ## color / label  set in `default.py`
        color, label = colors['data'], labels['data']
        bin_center = self._get_bin_center (self.edge)
        kwargs = self._set_errkwargs ('data')
        
        self.axhist.errorbar (bin_center, dataH,
                              yerr=sqrtH,
                              label=label, color=color,
                              ecolor=color,
                              markeredgecolor=color,
                              markerfacecolor=color,
                              **kwargs)
        return {'H':dataH, 'sqrtH':sqrtH}

    def _plot_cut_line (self, axis):

        ### leave if no cut lines need to be plotted
        if not hasattr (self, 'cut'): return
        if not self.cut: return

        ### change self.cut into list
        ### this feature is needed in case
        ### multiple cut lines are drawn
        if not toolbox.is_array (self.cut):
            self.cut = [self.cut]

        ### loop through each cut lines
        for cut in self.cut:
            ### show cut value line 
            axis.axvline (x=cut,
                          alpha=self.cut_alpha,
                          color=self.cut_color,
                          linestyle=self.cut_linestyle,
                          linewidth=self.cut_linewidth)
        return
    
    def _plot_a_ratio (self, ratiotype, ratio, error):

        ## color / label  set in `default.py`
        key = 'mcerr' if ratiotype == 'mc' else ratiotype
        color = colors[key]
        label = labels[key]
        bin_center = self._get_bin_center (self.edge)
        kwargs = self._set_errkwargs (ratiotype+'ratio')
        
        self.axratio.errorbar (bin_center, ratio,
                               yerr=error,
                               label=label, color=color,
                               ecolor=color,
                               markeredgecolor=color,
                               markerfacecolor=color,
                               **kwargs)
        return

    def _plot_ratios (self, mc, data):

        ### plot two ratio plots
        ### 1. data / mc in points
        ### 2. mc / mc in bands for MC uncertainties

        ### define denominator (i.e. mc)
        dH, sqrtdH = mc['H'], mc['sqrtH']
        for ratiotype in ['mc', 'data']:
            ## define the numerators
            numerator   = eval (ratiotype)            
            nH, sqrtnH = numerator['H'], numerator['sqrtH']
            ## calculate ratio / its errors
            ratio = np.divide (nH, dH)
            error = ratio * np.sqrt ((sqrtnH/nH)**2 + (sqrtdH/dH)**2)
            ## plot ratio
            self._plot_a_ratio (ratiotype, ratio, error)

        ## plot cut line
        self._plot_cut_line (self.axratio)
        return
    
    def plot (self, events, variable):

        self._clean_up ()
        
        ### upper (self.axhist) - distribution plot
        ## plot stacked MC and total MC histograms
        totalmc = self._plot_mc_histogram (events, variable)
        ## plot data points
        data = self._plot_data_histogram (events, variable,
                                          np.sum (totalmc['H']))
        ## format self.axhist
        self._format_hist ()

        ### lower (self.axratio) - ratio plot
        self._plot_ratios (totalmc, data)
        ## format self.axratio
        self._format_ratio (xLabel=labels[variable])
        
        return 


