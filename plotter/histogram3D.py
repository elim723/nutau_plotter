#!/usr/bin/env python

####
#### By Elim Thompson (09/26/2018)
####
#### This script contains a histogram2D class
#### which deals with all 3D distribution plots.
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
#### histogram3D inherited from Canvas2D
################################################
class Histogram3D (Canvas2D):

    def __init__ (self, sample, stretch=False, **kwargs):

        width_ratios = [30, 22, 1] if stretch else \
                       [30, 30, 1]

        ### default 3D distribution plots (python 2.7)
        super (Histogram3D, self).__init__ ((7.5, 3), (1, 3), 
                                            width_ratios=width_ratios,
                                            wspace=0.15, hspace=0,
                                            top=0.91, bottom=0.20,
                                            right=0.915, left=0.1)
        
        ### set any input style
        self.set_style (**kwargs)

        ### define observable
        self._xobs = None
        self._yobs = None
        self._zobs = None
        self._sample = sample.lower ()

        ### set subplot axes
        ## left = cascade-like histogram
        self.axcascade  = self.h.add_subplot (self.gs[0])
        ## middle = track-like histogram
        if self._sample=='greco' and not stretch:
            width = (0.915-0.1)*30/(30+30+1.)
            height = 0.91-0.20
            left = width + (width - width*0.67) + 0.063
            bottom = 0.20
            self.axtrack = self.h.add_axes ([left, bottom, width*0.67, height])
        else:
            self.axtrack = self.h.add_subplot (self.gs[1])
        ## right = colorbar
        self.axcolorbar = self.h.add_subplot (self.gs[2])

    @property
    def sample (self):
        return self._sample
        
    @sample.setter
    def sample (self, sample):
        self._sample = sample

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

        self.axcascade.cla ()
        self.axtrack.cla ()
        self.axcolorbar.cla ()
        return
    
    def _format (self, axis, xticks, yticks,
                 xlabel, ylabel, title,
                 xticklabels, yticklabels):

        ### set format
        self._format2D (axis, xticks, yticks, xlabel, ylabel,
                        xticklabels=xticklabels,
                        yticklabels=yticklabels )

        ### set title
        axis.set_title (title, fontsize=self.title_fontsize)
        return

    def _format_cascade (self):

        ### define axis
        axis = self.axcascade

        ### set x/y label and title
        xlabel = labels['reco_e']
        ylabel = labels['reco_cz'] 
        title  = r'Cascade-like'

        ### format
        self._format (axis, self.xticks, self.yticks    ,
                      xlabel, ylabel, title             ,
                      self.xticklabels, self.yticklabels )
        return

    def _format_track (self):

        ### define axis
        axis = self.axtrack

        ### set ticks
        xticks = self.xticks
        xticklabels = getattr (self, 'xticklabels', xticks)
        ##  for analysis histogram, ticks = edges
        ##  if greco and track, skip the first two tick
        if self.sample.lower ()=='greco':
            xticks = xticks[1:]
            xticklabels = xticklabels[1:]

        ### set x/y label and title
        xlabel = labels['reco_e']
        ylabel = None
        title  = r'Track-like'

        ### format
        self._format (axis, xticks, self.yticks,
                      xlabel, ylabel, title,
                      xticklabels, self.yticklabels)

        ### hide y axis ticks
        axis.get_yaxis ().set_ticks ([])

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

    def _get_content (self, histos, ptype=None):

        ''' determine content to be plotter

            :param ptype (str): either cascade or track
        
            4 options
            ---------
               1. numucc = distribution of numu CC (same for any member)
               2. nurate = total neutrino rates
               3. mufrac = muon fraction
               4. s/b    = signal / sqrt (background)
        '''

        ### option 1: histogram of a specific member
        if self.htype in histos.keys ():
            return histos[self.htype]

        ### option 2: summed neutrino histograms
        if self.htype == 'nurate':
            ## only add the neutrino in dtype 
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

    def _get_a_histogram (self, event):

        ''' get a 2D histogram for a data type. '''

        ### define values and weights
        xvalues = event[self.xobs]
        yvalues = event[self.yobs]
        zvalues = event[self.zobs]
        weights = event['weight']
        
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

        ### build a 3D histogram for each data type
        histograms = {'cascade':{},
                      'track'  :{}}

        ### loop through each member
        for member, event in events.items ():
            ## build a 3D histogram for this member
            histos = self._get_a_histogram (event)
            ## define cascade / track histograms
            cascade = histos[:,:,0]
            track   = histos[:,:,1]
            ## ignore first two track bins if greco
            if self.sample.lower () == 'greco':
                track = track[1:]
            ## store histograms
            histograms['cascade'][member] = cascade
            histograms['track'][member]   = track
        
        return histograms

    def _plot (self, histograms, ptype=None):

        ### return None if no valid ptype
        if not ptype or not hasattr (self, 'ax'+ptype):
            return None

        ### define axis
        axis = getattr (self, 'ax'+ptype)

        ### get content
        content = self._get_content (histograms[ptype])
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

        ### build 3D histograms
        histograms = self._get_histograms (events)
        
        ### left (self.axcascade) - cascade distribution
        self._plot (histograms, ptype='cascade')
        self._format_cascade ()

        ### middle (self.axtrack) - track distribution
        self._plot (histograms, ptype='track')
        self._format_track ()

        ### right (self.axcolorbar) - colorbar
        self._plot_colorbar (self.axcolorbar)
        
        return 


