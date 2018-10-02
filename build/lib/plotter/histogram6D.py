#!/usr/bin/env python

####
#### By Elim Thompson (09/26/2018)
####
#### This script contains a histogram6D class
#### which deals with all two 3D distribution
#### plots.
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
#### histogram6D inherited from Canvas2D
################################################
class Histogram6D (Canvas2D):

    def __init__ (self, **kwargs):

        ### default two 3D distribution plots (python 2.7)
        super (Histogram6D, self).__init__ ((7.5, 6), (2, 3), 
                                            width_ratios=[30, 30, 1],
                                            wspace=0.15, hspace=0.10,
                                            top=0.96, bottom=0.085,
                                            right=0.915, left=0.1)
        
        ### set any input style
        self.set_style (**kwargs)

        ### define observable
        self._xobs = None
        self._yobs = None
        self._zobs = None

        ###############################
        ##  GRECO 
        ###############################
        ##  upper left   = cascade-like
        self.axGcascade = self.h.add_subplot (self.gs[0])
        ##  upper middle = track-like
        width = (0.915-0.1)*30/(30+30+1.)
        height = (0.96-0.085)/2.*0.953
        left = width + (width - width*0.67) + 0.06
        bottom = 0.20*2.72
        self.axGtrack = self.h.add_axes ([left, bottom,
                                          width*0.68, height])

        ###############################
        ##  DRAGON
        ###############################
        ##  lower left  = cascade-like
        self.axDcascade = self.h.add_subplot (self.gs[3])
        ##  lower right = track-like
        self.axDtrack = self.h.add_subplot (self.gs[4])

        ###############################
        ##  Colorbars
        ###############################
        ##  upper = GRECO Colorbar
        self.axGcb   = self.h.add_subplot (self.gs[2])
        ##  lower = DRAGON Colorbar
        self.axDcb   = self.h.add_subplot (self.gs[5])

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

        self.axGcascade.cla ()
        self.axGtrack.cla ()
        self.axDcascade.cla ()
        self.axDtrack.cla ()
        self.axGcb.cla ()
        self.axDcb.cla ()
        return
    
    def _format (self, axis, xticks, yticks,
                 xlabel, ylabel, title,
                 xticklabels, yticklabels):

        ### set format
        self._format2D (axis, xticks, yticks, xlabel, ylabel,
                        xticklabels=xticklabels,
                        yticklabels=yticklabels )
 
        ### set title
        if title: axis.set_title (title, fontsize=self.title_fontsize)
        return

    def _format_cascade (self, sample):

        ### define axis
        axis = getattr (self, 'ax'+sample[0].upper () + 'cascade')

        ### set x/y label and title
        xlabel = labels['reco_e'] if sample=='dragon' else None
        ylabel = labels['reco_cz'] 
        title  = r'Cascade-like' if sample=='greco' else None

        ### set ticks
        xticks = self.xticks if sample=='dragon' else self.xedges
        xticklabels = getattr (self, 'xticklabels', xticks)
        if sample=='greco': xticklabels = []

        ### format
        self._format (axis, xticks, self.yticks    ,
                      xlabel, ylabel, title        ,
                      xticklabels, self.yticklabels )
        return

    def _format_track (self, sample):

        ### define axis
        axis = getattr (self, 'ax'+sample[0].upper () + 'track')

        ### set ticks
        xticks = self.xticks if sample=='dragon' else self.xedges
        xticklabels = getattr (self, 'xticklabels', [])
        ##  for analysis histogram, ticks = edges
        ##  if greco and track, skip the first two tick
        if sample.lower ()=='greco':
            xticks = xticks[2:]
            xticklabels = []

        ### set x/y label and title
        xlabel = labels['reco_e'] if sample=='dragon' else None
        ylabel = None
        title  = r'Track-like' if sample=='greco' else None

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

    def _get_histograms (self, events, sample):

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
            if sample.lower () == 'greco':
                track = track[2:]
            ## store histograms
            histograms['cascade'][member] = cascade
            histograms['track'][member]   = track
        
        return histograms

    def _plot (self, histograms, sample, ptype=None):

        ### return None if no valid ptype
        if not ptype: return None

        ### get axis
        axis = getattr (self, 'ax'+sample[0].upper ()+ptype, None)

        ### if None axis, return
        if not axis: return None

        ### get content
        content = self._get_content (histograms[ptype])
        content = np.nan_to_num (content)

        ### set edges
        xticks = deepcopy (self.xedges)
        if sample.lower ()=='greco' and \
           ptype=='track': xticks = xticks[2:]

        ### plot content
        self._plot2D (axis, content, xedges=xticks)

        return 

    def plot (self, sample, events):

        ### set edges
        self._set_edges ()

        ### build 3D histograms
        histograms = self._get_histograms (events, sample)
        
        ### left - cascade distribution
        self._plot (histograms, sample, ptype='cascade')
        self._format_cascade (sample)

        ### middle - track distribution
        self._plot (histograms, sample, ptype='track')
        self._format_track (sample)

        ### right (self.axcolorbar) - colorbar
        cbaxis = getattr (self, 'ax' + sample[0].upper ()+'cb')
        self._plot_colorbar (cbaxis)
        
        return 


