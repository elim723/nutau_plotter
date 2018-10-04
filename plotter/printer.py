#!/usr/bin/env python

####
#### By Elim Thompson (09/19/2018)
####
#### This script collects all printers which print info
#### as a plotting script is running.
####
#########################################################

#### import standard packages
from __future__ import print_function
import numpy as np

#########################################################
#### printers for collect events (collect_events.py)
#########################################################
def collectevents_print_header (outfile):

    ''' print header

        :param outfile  (str) : location of output pickled dictionary
    '''

    print ('################# collect baseline events ##################')
    print ('#### Welcome! I will collect the events your want in')
    print ('#### no time. The output dictionary with the correct')
    print ('#### weights will be located at ...')
    print ('####       {0}'.format (outfile))
    print ('####')
    return

def collectevents_print_params (params, do_print=False, forsyseffects=False):

    ''' print input values of parameters to be used

        :param params   (dict): injected values to be used for weighting
        :param forsyseffects (bool): If true, params are baseline values
                                              for systematic effects on
                                              baseline histograms
        :param do_print (bool): if True, do print
    '''

    if do_print:
        text = 'baseline ' if forsyseffects else ''
        print ('#### Your input {0}values of parameters are ...'.format (text))
        line = '####       {0:<16}: {1:4f}'
        for param in params:
            print (line.format (param, round (params[param], 4)))
        print ('####')
    return

def collectevents_print_library (members, do_print=False):

    ''' print library progress

        :param members  (list): data types included in the library
        :param do_print (bool): if True, do print
    '''

    if do_print:
        statement = '#### Collecting {0} members in a library ...'
        print (statement.format (len (members)))
        print ('#### NOTE: This may take some time. Be patient!')
        print ('####')
    return

def collectevents_print_dtype (dtype, unweighted, prerate, postrate, do_print=False):

    ''' print information for each member / data type

        :param dtype      (str): name of the member / data type
        :param unweighted (int): number of unweighted simulated events
        :param prerate  (float): total rates before hyperplane is applied
        :param postrate (float): total rates after hyperplane is applied
        :param do_print  (bool): if True, do print
    '''

    if do_print:
        print ('#### +-----------------------------------------------')
        print ('#### | working on {0} ...'.format (dtype))
        print ('#### |   unweighted counts : {0} events'.format (unweighted))
        print ('#### |   pre-hplaned  rates: {0} mHz'.format (round (prerate , 4)))
        print ('#### |   post-hplaned rates: {0} mHz'.format (round (postrate, 4)))
        print ('#### +-----------------------------------------------')
        print ('####')
    return

def collectevents_print_ender (dtime, outfile):

    ''' print ender

        :param dtime   (float): time taken for this script to run
        :param outfile   (str): location of output pickled dictionary
    '''

    print ('#### My job is done, and your dictionary is at ...')
    print ('####       {0}'.format (outfile))
    print ('#### Nice working with you! Bye :)')
    print ('############################################################')
    print ('#### It took {0} minutes ...'.format (dtime))
    return

#########################################################
#### printers for collect sysevents (collect_sysevents.py)
#########################################################
def collectsys_print_header (outfile):

    ''' print header

        :param outfile  (str) : location of output pickled dictionary
    '''

    print ('############### collect sys effect events ##################')
    print ('#### Welcome! I will collect the events your want in')
    print ('#### no time. The output dictionary with the correct')
    print ('#### weights will be located at ...')
    print ('####       {0}'.format (outfile))
    print ('####')
    return

def collectsys_print_params (params, parameters, do_print=False):

    ''' print input values of parameters to be used

        :param params   (dict): injected values to be used for weighting
        :param forsyseffects (bool): If true, params are baseline values
                                              for systematic effects on
                                              baseline histograms
        :param do_print (bool): if True, do print
    '''
    
    collectevents_print_params (params, do_print=do_print,
                                forsyseffects=True)
    ### report what parameters to be included
    nparams = len (parameters)
    print ('#### Output will contains {0} sets of weights.'.format (nparams))
    print ('#### Each of them corresponds to a change')
    print ('#### in the following parameters:')
    for param, delta in parameters.items ():
        print ('####     -- {0:<16}: {1:4f}'.format (param,
                                                     np.round (delta, 4)))
    print ('####')
    return

def collectsys_print_library (members, do_print=False):

    ''' print library progress

        :param members  (list): data types included in the library
        :param do_print (bool): if True, do print
    '''

    if do_print:
        statement = '#### Collecting {0} members in a library ...'
        print (statement.format (len (members)))
        print ('#### NOTE: This may take some time. Be patient!')
        print ('####')
    return

def collectsys_print_dtype (dtype, events, do_print=False):

    ''' print information for each member / data type

        :param dtype     (str): name of the member / data type
        :param events   (dict): dictionary of all event info
        :param do_print (bool): if True, do print
    '''

    if do_print:
        print ('#### +-----------------------------------------------')
        print ('#### | working on {0} ...'.format (dtype))

        ## obtain rates for each weight key
        for key, value in events.items ():
            if 'weight' in key:
                weightkey = ' '.join (key.split ('_')[:-1])
                weight = round (np.sum (value)*1000., 4)
                print ('#### |   {0:<20} rate: {1} mHz'.format (weightkey, weight))
        print ('#### +-----------------------------------------------')
        print ('####')
    return

def collectsys_print_ender (dtime, outfile):

    collectevents_print_ender (dtime, outfile)
    return

#########################################################
#### printers for selection plots (plot_greco1D.py)
#########################################################
def selection_printer_header (outdir, sample):
    
    ''' print header for selection plots in plot_greco1D.py

        :param outdir (str): folder to output plots
        :param sample (str): GRECO or DRAGON?
    '''

    head = '################# plot {0:6} 1D selection  #################'
    print (head.format (sample.upper ()))
    print ('#### Welcome! This script produces 1D distribution plots')
    print ('#### from GRECO Level 4 - 7. Your plots will be found in')
    print ('####       {0}'.format (outdir))
    print ('####')
    return

def selection_printer_events (events):

    ''' print rates for selection plots in plot_greco1D.py '''

    print ('#### NOTE: These rates are higher than the ones quoted')
    print ('####       in the paper because these rates are pre-cut')
    print ('####       We show distribution plots *before* the cut')
    print ('####       is applied.')
    print ('####       If you want to apply cut, see `get_events ()`')
    print ('####')

    line = '#### |      {0:<7}: {1} mHz'
    for level, members in events.items ():
        print ('#### +----------------------------------')
        print ('#### | At {0} ...'.format (level))
        ## print rates from each MC member
        mc = 0.
        for member, event in members.items ():
            # skip data for now
            if member == 'data': continue
            rate = np.sum (event['weight'])
            rate = round (rate*1000., 4)
            # add up to total mc rate
            mc += rate
            print (line.format (member, rate))
        print ('#### |      ============================')
        ## print mc total rate
        print (line.format ('MC', mc))
        ## print data rate
        data = np.sum (members['data']['weight'])
        data = round (data*1000., 4)
        print (line.format ('data', data))
        ## print data / mc
        ratio = round (data / mc, 4)
        print ('#### |      data/MC: {0}'.format (ratio))                                                                                          
        print ('#### |      ============================')
        print ('#### +----------------------------------')
        print ('####')
        
def selection_printer_ender (outdir):

    ''' print ender for selection plots in plot_greco1D.py

        :param outdir (str): folder to output plots
    '''
    
    print ('#### My job is done. Your plots can be found in')
    print ('####       {0}'.format (outdir))
    print ('#############################################################')
    print ('####')
    return

#########################################################
#### printers for classification plots (plot_pid.py)
#########################################################
def classification_printer_header (outdir):
    
    ''' print header for classification plots in plot_classification.py

        :param outdir (str): folder to output plots
    '''

    print ('################# plot PID classification #################')
    print ('#### Welcome! This script produces PID classification')
    print ('#### plots from both samples at their final levels.')
    print ('#### The output plots show the ratio of track events')
    print ('#### in an energy slice to total number of events in')
    print ('#### the same energy slice. Your plots will be found in')
    print ('####       {0}'.format (outdir))
    print ('####')
    return

def classification_printer_events (outfile):

    ''' print event info for classification plots in
        plot_classification.py
    
        :param outfile (str): location of output file
                              with track fractions
    '''
    
    print ('#### You did not provide an input file. Getting events')
    print ('#### from scratch may take a while.')
    print ('#### Your patience is appreciated :)')
    print ('####')
    print ('#### The track fractions for this PID classification')
    print ('#### plot is stored as a pickled dictionary at')
    print ('####       {0}'.format (outfile))
    print ('####')
    print ('#### Next time, feed the file as `infile` to save time!')
    print ('####')
    return

def classification_printer_rates (events):

    ''' print rates for classification plots in
        plot_classification.py

        :param events (dict): event rates
    '''
    
    line = '#### |      {0:<7}: {1} mHz'
    for sample in events.keys ():
        print ('#### +----------------------------------')
        print ('#### | sample {0} rate info'.format (sample))
        print ('#### |      ===================')
        for member, event in events[sample].items ():
            ## print rates from each member
            rate = np.sum (event.w)
            rate = round (rate*1000., 4)
            print (line.format (member, rate))
        print ('#### +----------------------------------')
        print ('####')
    return

def classification_printer_ender (outdir):

    ''' let user know where the plots are for classification
        plots in plot_classification.py
    '''
    selection_printer_ender (outdir)
    return

#############################################################
#### printers for resolution plots (plot_resolutions.py)
#############################################################
def resolution_printer_header (outdir):
    
    ''' print header for resolution plots in plot_resolution.py

        :param outdir (str): folder to output plots
    '''

    print ('##################### plot resolution #####################')
    print ('#### Welcome! This script produces pegleg resolution')
    print ('#### plots from both samples at their final levels.')
    print ('#### A total of six resolution plots will be created:')
    print ('####    -- linear reco energy vs linear mc energy')
    print ('####    -- log10  reco energy vs log10  mc energy')
    print ('####    -- reoc - true cos zenith vs linear mc energy')
    print ('####    -- reoc - true cos zenith vs log10 mc energy')
    print ('####    -- reoc - true zenith (deg) vs linear mc energy')
    print ('####    -- reoc - true zenith (deg) vs log10 mc energy')
    print ('#### All outputs will be found at ')
    print ('####       {0}'.format (outdir))
    print ('####')
    return

def resolution_printer_events (outfile):

    ''' print event info for resolution plots in
        plot_resolution.py
    
        :param outfile (str): location of output file
                              pegleg resolutions
    '''
    
    print ('#### You did not provide an input file. Getting events')
    print ('#### from scratch may take a while.')
    print ('#### Your patience is appreciated :)')
    print ('####')
    print ('#### The pegleg resolutions for this resolution')
    print ('#### plot are stored as a pickled dictionary at')
    print ('####       {0}'.format (outfile))
    print ('####')
    print ('#### Next time, feed the file as `infile` to save time!')
    print ('####')
    return

def resolution_printer_rates (events):

    ''' print rates for resolution plots in
        plot_resolution.py

        :param events (dict): event rates
    '''
    
    classification_printer_rates (events)
    return

def resolution_printer_res (combination, sample,
                            member, res, do_print=False):

    ''' print the resolution info of a given
        x/y combination, a given data type, a
        given sample

        :param combination  (str): 'xaxis'_'yaxis'
        :param sample       (str): dragon / greco
        :param member       (str): data type 
        :param res         (dict): resolution info
    '''

    if do_print:
        print ('+----------------------------------------------------')
        print ('| {0} {1} ...'.format (sample, member))
        print ('|     xaxis: {0}'.format (combination.split ('_')[0]))
        print ('|     yaxis: {0}'.format (combination.split ('_')[1]))
        for key, value in res.items ():
            print ('|        {0:7}: {1}'.format (key, value))
    return
    
def resolution_printer_ender (outdir):

    ''' let user know where the plots are for resolution
        plots in plot_resolution.py
    '''
    selection_printer_ender (outdir)
    return

#############################################################
#### printers for sig / sqrt (bg) plots (plot_sig_to_by.py)
#############################################################
def sigtobg_printer_header (outdir, sample):
    
    ''' print header for sig / sqrt bg plots in plot_sig_to_bg.py

        :param outdir (str): folder to output plots
        :param sample (str): either GRECO or DRAGON
    '''

    print ('################# plot sig / sqrt (bg) ####################')
    print ('#### Welcome! This script produces a sig / sqrt (bg)')
    print ('#### plots from {0} sample. By default, signal is'.format (sample.upper ()))
    print ('#### oscillated nutau and background is everything else.')
    print ('#### All outputs will be found at ')
    print ('####       {0}'.format (outdir))
    print ('####')
    return

def sigtobg_printer_events (events):

    ''' print rates for sig / sqrt (bg) plots in
        plot_sig_to_bg.py

        :param events (dict): event rates
    '''

    line = '#### |      {0:<8}: {1}'
    #### start printing
    print ('#### +----------------------------------')
    print ('#### | At final level within bins ...')
    mc = 0.
    for member, event in events.items ():
        ## print rates from each member
        rate = round (np.sum (event['weight']), 0)
        mc += rate
        print (line.format (member, rate))
    print ('#### |      ============================')
    print (line.format ('total MC', mc))
    print ('#### +----------------------------------')
    print ('####')
    return

def sigtobg_printer_ender (outdir):

    ''' let user know where the plots are for sig / sqrt (bg)
        plots in plot_sig_to_bg.py
    '''
    selection_printer_ender (outdir)
    return

#################################################################
#### printers for muon histogram plots (plot_muon_histograms.py)
#################################################################
def muhist_printer_header (outdir, samples):
    
    ''' print header for muon histograms in plot_muon_histograms.py

        :param outdir (str): folder to output plots
        :param sample (str): either GRECO or DRAGON
    '''

    text = ''
    for index, sample in enumerate (samples):
        text += sample.upper ()
        ext = ' and ' if index==len (samples)-2 else \
              ''     if index==len (samples)-1 else \
              ', ' 
        text += ext

    print ('################# plot muon histogram #####################')
    print ('#### Welcome! This script produces muon histogram(s)')
    print ('#### from {0} sample.'.format (text))
    print ('#### All outputs will be found at ')
    print ('####       {0}'.format (outdir))
    print ('####')
    return

def muhist_printer_events (events):

    ''' print rates for muon histograms in
        plot_muon_histograms.py

        :param events (dict): event rates
    '''

    line = '#### |      {0:<8}: {1}'
    #### start printing
    for sample, event in events.items ():
        print ('#### +----------------------------------')
        print ('#### | At final level of {0} ...'.format (sample))
        ## print rates from each member
        rate = round (np.sum (event['muon']['weight']), 0)
        print (line.format ('muon', rate))
    print ('#### +----------------------------------')
    print ('####')
    return

def muhist_printer_ender (outdir):

    ''' let user know where the plots are for muon histograms
        plots in plot_muon_histograms.py
    '''
    selection_printer_ender (outdir)
    return

#################################################################                                                           
#### printers for syseffect plots (plot_sys_effects.py)                                                                     
#################################################################                                                           
def syseffects_printer_header (outdir):

    ''' print header for systematic effecs in plot_sys_effects.py

        :param outdir (str): folder to output plots
        :param sample (str): either GRECO or DRAGON
    '''

    print ('################# plot systematic effects #################')
    print ('#### Welcome! This script produces the systematic effecst')
    print ('#### on histograms in percentage change.')
    print ('#### All outputs will be found at ')
    print ('####       {0}'.format (outdir))
    print ('####')
    return

def syseffects_printer_events (events):

    ''' print rates for muon histograms in plot_sys_effects.py

        :param events (dict): event rates
    '''

    for dtype, event in events.items ():
        collectsys_print_dtype (dtype, event, do_print=True)

    return

def syseffects_printer_ender (outdir):

    ''' let user know where the plots are for systematic
        effects on histograms in plot_sys_effects.py                                                                                        
    '''
    selection_printer_ender (outdir)
    return

#########################################################
#### printers for greco numu plots (plot_greco_numu.py)
#########################################################
def contour_printer_header (outdir):

    ''' print header for contour in plot_greco_numu.py

        :param outdir (str): folder to output plots
    '''

    print ('#################### plot 90% contour #####################')
    print ('#### Welcome! This script produces a 90% contour plots')
    print ('#### of standard atmospheric oscillation measurements')
    print ('#### from various experiments. Your plot will be in')
    print ('####       {0}'.format (outdir))
    print ('####')
    return

def contour_printer_ender (outdir):

    ''' let user know where the plots are located
        from plot_greco_numu.py
    '''
    selection_printer_ender (outdir)
    return

###########################################
#### printers for create_nufile.py
###########################################
def createnu_print_header (infile, outfile, do_print=False):

    ''' print header before modifying text file if verbose

        :param infile    (str): nuparam textfile template
        :param outfile   (str): output file with users settings
        :param do_print (bool): If True, do print.
    '''

    if do_print:
        print ('################ modification starts ################')
        print ('#### The nuparam textfile template is ...')
        print ('####      {0}'.format (infile))
        print ('#### Your output nuparam textfile will be ...')
        print ('####      {0}'.format (outfile))
        print ('####')
    return

def createnu_print_line (param, line,do_print=False):

    ''' print updated line of a given parameter if verbose

        :param param     (str): name of the parameter to be changed
        :param line      (str): the updated line after changes
        :param do_print (bool): If True, do print.
    '''

    if do_print:
        print ('#### +---------------------------------------------')
        print ('#### | {0} is updated to be ...'.format (param))
        print ('#### | {0}'.format (' '.join (line.strip ().split ())))
        print ('#### +---------------------------------------------')
    return

def creatnu_print_ender (outfile, do_print=False):

    ''' print end of work if verbose

        :param outfile   (str): output file with users settings
        :param do_print (bool): If True, do print.
    '''

    if do_print:
        print ('####')
        print ('#### My job is done.')
        print ('#### Your output nuparam textfile is now at ...')
        print ('####      {0}'.format (outfile))
        print ('#### Nice working with you! Bye :)')
        print ('#####################################################')
        print ('####')
    return
