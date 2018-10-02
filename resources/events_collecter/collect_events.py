#!/usr/bin/env python

####
#### By Elim Cheung (09/18/2018)
####
#### This script creates a dictionary with all baseline events
#### and their corresponding weights calcuated based on a given
#### `nufile`. Only baseline events are stored in the outputs.
#### Their weights are calculated from the `injected` (2nd
#### column) values of the input `nufile`. A hyperplane object
#### is also created to modify the weights.
####
#### Command line to run:
#### $ python step2_events.py
####            --outfile <address/to/outputfile.p>
####		--nufile <nuparam file> --verbose 1
####            (optional: --oscnc)
####
#### NOTE: 1. other constants (e.g. edges) are imported from
####          `../analysistools/misc.py` by default. Change the
####          histogram bin edges in this script if needed.
####       2. `inverted` is set to be `False` so normal
####          mass heirarchy is assumed.
####       3. `matter` is set to be `True` so matter effect
####          by Earth is included. 
####
#### This script does the following procedures:
####
#### 0. Define `nuparams` instance based on input `nufile`
#### 1. a. Set up a `Library` instance which includes all
####       baseline events with the corresponding weighters
####    b. A `Hyperplane` instance is also created with all
####       fits per bin per data type (or member)
#### 2. Collect baseline events as a dictionary
####    Their weights are also calculated at this step
#### 3. Save the dictionary as a pickled file
####
###############################################################

from __future__ import print_function
from optparse import OptionParser
import time, os, cPickle, sys
from copy import deepcopy
import numpy as np

## customized classes
from analyzer.library import Library
from analyzer.nuparams import Nuparams
from analyzer.misc import eedges, zedges, pedges
from analyzer.misc import Toolbox, Map

## import printers
from plotter.printer import collectevents_print_header
from plotter.printer import collectevents_print_params
from plotter.printer import collectevents_print_library
from plotter.printer import collectevents_print_dtype
from plotter.printer import collectevents_print_ender

## ignore RuntimeWarning
import warnings
warnings.filterwarnings("ignore")

###########################################
#### default variables
###########################################
sample = 'greco'
## these edges are imported from `misc.py`
## change here if needed
edges   = {'e':eedges, 'z':zedges, 'p':pedges}
## these members are included for GRECO 
members = ['numucc', 'nuecc', 'nutaucc', 'numunc',
           'nuenc', 'nutaunc', 'noise', 'muon']
## default including matter effect assuming
## normal mass hierarchy
matter  = True
inverted = False
## default path
ppath   = os.path.dirname(os.path.abspath( __file__ )) + \
          '/../pickled_data/greco/final/'
## keys in new dictionary
keys = ['reco_logen', 'reco_coszen',
        'logen', 'coszen', 'weight']
if sample=='greco':
    keys += ['tracklength']
if sample=='dragon':
    keys += ['deltaLLH']

###########################################
#### functions to deal with users inputs
###########################################
def check_args (nufile, ppath, outfile):

    ''' check if paths/files exist before work

        :param nufile  (str): location of input nu textfile
        :param ppath   (str): path to all pickled data files
        :param outfile (str): location of output pickled object
    '''
    
    toolbox = Toolbox ()
    ## does nufile exist?
    toolbox.check_file (nufile)
    ## is ppath valid?
    toolbox.check_path (ppath)
    ## is outdir valid?
    toolbox.check_path (os.path.split (outfile)[0])
    return

def define_args (options):

    ''' define variables for this script from users

        :param options (dict): options parsed from users
    '''

    ppath = options.ppath
    nufile, outfile = options.nufile, options.outfile
    check_args (nufile, ppath, outfile)

    oscnc   = options.oscnc
    verbose = options.verbose
    return ppath, nufile, outfile, oscnc, verbose

###########################################
#### functions to get weights
###########################################
def get_library (members, ppath, edges, nuparams,
                 matter, oscnc, do_print=False):

    ''' set up library, its weighters, and a hyperplane

        :param members  (list): data types included
        :param ppath    (str) : path to the pickled data files
        :param edges    (dict): histogram bin edges
        :param nuparams (`nuparams`): a `Nuparams` instance
        :param matter   (bool): If True, matter effect included
        :param oscnc    (bool): If True, oscillate NC
        :param do_print (int) : If > 0 , print progress in `library`

        :return lib        (`Library`): a library instance
        :return hplanes (`Hyperplane`): a hyperplane instance
    '''

    ## report to user
    collectevents_print_library (members, do_print=do_print)
    
    ## initialize library instance
    ## NOTE: set verbose = 2 for `Library`
    ##       to monitor progress
    lib = Library (members, ppath,
                   edges  =edges,
                   verbose=0)
    ## set up weighters for baseline
    ## members using injected values
    injected = nuparams.extract_params ('injected')
    lib.set_weighters (injected,
                       matter=matter,
                       oscnc =oscnc)
    
    ## initialize hyperplane instance
    ## NOTE: set verbose = 1 or 2 for `Hyperplane`
    ##       to monitor progress
    hplanes = lib.get_hplanes (nuparams,
                               verbose=0)
    return lib, hplanes

def obtain_factors (dtype, hplane, nuparams):

    ''' apply a hyperplane of a given data type (or
        member) to get the factors from the injected
        discrete systematic parameters

        :param dtype           (str): name of this member
        :param hplane (`Hyperplane`): the `Hyperplane` instance
                                      for this member
        :param nuparams (`Nuparams`): a `Nuparam` instance from
                                      the input nufile

        :return factors   (np.array): factors for each bin in
                                      the same shape as the
                                      histogram
    '''
    
    ## store injected values of discrete parameters as a dict
    params = {sys:nuparams (sys).injected
              for sys in nuparams.get_hplaned_dparams (dtype)}
    ## obtain factors for all bins by applying the hplane
    return hplane.apply (params)

def apply_factors (dtype, hplane, nuparams, ddict, weights):

    ''' modify weights of a given member using factors
        obtained by evaluate its hyperplane at the injected
        values of the discrete parameters

        :param dtype           (str): name of this member
        :param hplane (`Hyperplane`): the `Hyperplane` instance
                                      for this member
        :param nuparams (`Nuparams`): a `Nuparam` instance from
                                      the input nufile
        :param ddict          (dict): observables/variables from
                                      all events of this member
        :weights          (np.array): weights before modification
                                      from hyperplane

        :weights (np.array): weights after modification by hplanes
    '''
    
    ## return weights if no hyperplane object
    ## e.g. noise
    if not hplane: return weights
    
    ## get factors based on nuparams and hplane
    factors = obtain_factors (dtype, hplane, nuparams)

    ## rebuild baseline histogram 
    e = np.array (ddict.reco.e)
    z = np.array (ddict.reco.z)
    p = np.array (ddict.reco.pid)
    w = np.array (weights)
    finite = np.logical_and (np.isfinite (e), np.isfinite (z))
    finite = np.logical_and (finite, np.isfinite (p))
    finite = np.logical_and (finite, np.isfinite (w))
    events = np.array ([e[finite], z[finite], p[finite]]).T

    ## define modified weights
    modweights = deepcopy (weights)
    
    ## loop through each event
    for j, event in enumerate (events):
        ## which histogram bin does this event fall into ?
        indices = [ np.digitize ([event[i]], edges[n])[0]-1
                    for i, n in enumerate (['e', 'z', 'p']) ]
        ## multiply this weight by the factor of this bin
        modweights[j] *= factors[indices[0]][indices[1]][indices[2]]
        
    return modweights

def get_events (lib, hplanes, nuparams, do_print=False):

    ''' get event dictionary with weights calculated
        based on injected values and modified by hplane

        :param lib      (`Library`)   : a `Library` instance with
                                        all events
        :param hplanes  (`Hyperplane`): a `Hyperplane` instance with
                                        fits per bin per member
        :param nuparams (`Nuparams`)  : a `Nuparams` instance with
                                        user's input as `injected`
        :param do_print (bool)        : if True, print progress

        :return dictionary (dict): observables / variables from
                                   all members with weights
    '''
    
    injected = nuparams.extract_params ('injected')
    dictionary = Map ({})

    ### loop through each data type / member
    for dtype in lib.members:
        ## define member / weighter
        member = lib (dtype)
        weighter = lib.weighters [dtype]
        ## get weights based on injected values
        weights = member.get_weights (injected, weighter=weighter)
        ## modify weights from hyperplane of this member
        hplane = hplanes['nunc'] if 'nc' in dtype else \
                 hplanes[dtype] if dtype in hplanes else None
        modweights = apply_factors (dtype, hplane, nuparams,
                                    member (), weights)
        ## set weights for this member
        member ().w = modweights
        ## print progress in mHz for rates
        collectevents_print_dtype (dtype, len (weights),
                                   np.sum (weights)*1000.,
                                   np.sum (modweights)*1000.,
                                   do_print=do_print)
        ## store in dictionary
        dictionary [dtype] = convert_dictionary (member ())
    return dictionary

def convert_dictionary (ddict):

    ### initialize new dictionary
    ndict = {}

    ### loop through keys to be pickled
    for key in keys:

        value = []

        ## true info
        if key in ['logen', 'coszen']:
            # move on if no true info
            if 'mc' not in ddict: continue
            # determine value
            value = np.log10 (ddict.mc.e)  \
                    if 'logen' in key else \
                       np.cos (ddict.mc.z)
        elif key in ['tracklength', 'deltaLLH']:
            ## PID info
            value = ddict.reco.pid 

        elif 'reco' in key:
            ## reco info
            value = np.log10 (ddict.reco.e) \
                    if 'logen' in key else  \
                    np.cos (ddict.reco.z) 

        elif key=='weight':
            ## weight
            value = ddict.w

        ## pickle
        if len (value)>0:
            ndict[key] = np.array (value).flatten ()

    return ndict

###########################################
#### main function
###########################################
if __name__ == "__main__":

    ### parse user's options
    usage = "%prog [--outfile library.p --nufile nufiles/nuparams_template.txt --verbose 1 (--oscnc)]"
    parser = OptionParser(usage=usage)
    parser.add_option ("--ppath", type="string", default=ppath,
                       help = "directory of pickled data")
    parser.add_option ("--outfile", type="string", default='library.p',
                       help = "address and name for out put file (with extension .p)")
    parser.add_option ("--nufile", type="string", 
                       default='~/nuisance_parameters_settings.txt',
                       help = "a textfile with a table of nuisance parameter settings")
    parser.add_option ("--verbose", type="int", default=1,
                       help = "0 for no print out; 2 for detailed print out")
    parser.add_option ("--oscnc", action="store_true", default=False,
                       help = "oscillate NC neutrinos")
    (options, args) = parser.parse_args()

    ### print header
    ppath, nufile, outfile, oscnc, verbose = define_args (options)
    collectevents_print_header (outfile)
    start_time = time.time()

    ### set up parameter values
    nuparams = Nuparams (nufile, isinverted=inverted)
    collectevents_print_params (nuparams.extract_params ('injected'),
                                do_print=verbose)

    ### set up library and hyperplanes
    lib, hplanes = get_library (members, ppath, edges, nuparams,
                                matter, oscnc, do_print=verbose)

    ### get dictionary of baseline events with weights
    bevents = get_events (lib, hplanes, nuparams, do_print=verbose)

    ### save dictionary
    with open (outfile, 'wb') as f:
        cPickle.dump (bevents, f, protocol=2)
    f.close ()

    ### print ender
    dtime = round ((time.time() - start_time)/60., 2) # in minute
    collectevents_print_ender (dtime, outfile)
