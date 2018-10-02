#!/usr/bin/env python

####
#### This script merge the pickled files
#### from different years into one
####
##############################################
import numpy as np
import cPickle

is_array = lambda var: isinstance (var, list) or isinstance (var, np.ndarray)

def merge (arrays):
    
    carray = np.array (arrays[0])
    if len (arrays[0]) > 1:
        for i in np.arange (len (arrays)-1):
            if is_array (arrays[i+1]):
                carray = np.concatenate ((carray, np.array (arrays[i+1])))
    return carray

def concat_dicts (d1, d2, d3):

    cdict = {}
    for key in d1.keys():
        if is_array (d1[key]):
            cdict[key] = merge ([ d[key] for d in [d1, d2, d3] ])
            continue
        if isinstance (d1[key], dict):
            cdict[key] ={}
            for skey in d1[key].keys():
                cdict[key][skey] = merge ([ d[key][skey] for d in [d1, d2, d3] ])
    return cdict


folder = '/data/user/elims/for_nutau_paper/resources/pickled_data/dragon/preweighted/'

for ptype in ['data', 'muon']:
    
    outfile = 'level6_postreco_'+ptype+'.pckl'
    mdict = {}
    for year in ['2012', '2013', '2014']:
        filename = 'level6_postreco_'+ptype+'_'+year+'.pckl'
        with open (folder+filename, 'rb') as f:
            mdict[year] = cPickle.load (f)
        f.close ()
    
    
    cdict = concat_dicts (mdict['2012'], mdict['2013'], mdict['2014'])
    with open (folder + outfile, 'wb') as f:
        cPickle.dump (cdict, f, protocol=2)
    f.close ()
