#!/usr/bin/env python

####
#### By Elim Cheung (09/17/2018)
####
#### This script modifies the nuparam_textfile for
#### reweighting when making distribution plots.
####
#### For each parameter, only the injected values
#### are changed (if asked). Distribution plots
#### will be generated by injected values. The
#### seeded values will be used for making hyperplane
#### and should be the baseline values for detector
#### systematics.
####
####
#### Command line to run:
#### $ python step1_nuparams.py
####          --outdir <path/to/out/folder>
####         (--verbose)
####
#### Note: user can also change the textfile
#### template itself directly without this script.
####
#### The script does the following steps
####
#### 0. ask user's inputs
####     For each parameter,
####      Q: Do you want to change the default value?
####         If yes, what is your new value?
####
#### 1. create a new nufile with the users values.
####
#########################################################

from __future__ import print_function
from optparse import OptionGroup, OptionParser
import numpy as np
import sys, os

### import printers
from plotter.printer import createnu_print_header
from plotter.printer import createnu_print_line
from plotter.printer import createnu_print_ender

from analyzer.Nuparams import nuparams

###########################################
#### define variables
###########################################
indir = os.path.dirname (os.path.realpath(__file__))
nufile = indir + '/../nufiles/nuparams_template.txt'

### define new injected values for non-osc parameters
nuParams = nuparams.Nuparams (nufile)
### all available parameters
params = sorted (nuParams.get_all_params ())

###########################################
#### function to ask for input values
###########################################
def ask_param (param):

    ''' ask for user input for a given param 

        :param param (str): name of the given parameter
        
        :return value (str): value input by user
    '''
    
    ### report default param value to user
    statement = '#### The current value of {0} is {1}.'
    print (statement.format (param, nuParams (param).value))
    
    ### ask if user wants to change the current value
    question = '#### Do you want to change the default value (y/n) ? '
    warning  = '#### Your input {0} is invalid. Please answer y or n!'
    while True:
        change = raw_input (question)
        ## only accept `y` or `n`
        if change in ['y', 'n']: break
        print (warning.format (change))

    ### return None if user doesn't change the value
    if change == 'n': return None
    
    ### print warning for dm31
    if param == 'dm31':
        print ('#### NOTE: dm31 is in unit of 10^-3 eV^2')
        print ('####       e.g. if dm31 = 2.526e-3, enter `2.526`.')
    
    ### ask for a valid number
    question = '#### What would the new value for {0} be ? '
    warning  = '#### Your input {0} is invalid. Please enter a number!'
    while True:
        value = raw_input (question.format (param))
        ## accept the input if `n` or numbers
        if value.replace ('.', '', 1).replace ('-', '', 1).isdigit () or \
           value == 'n': break
        print (warning.format (value))
    
    ### return None if user change his/her mind
    if value=='n':
        out = '#### You changed your mind.. {0} remains as the default value.'
        print (out.format (param))
        return None
    
    ### report input value
    print ('#### {0} input value is {1}.'.format (param, value))
    return value

def ask_inputs (outfile):

    ''' ask user for his/her values of all parameters

        :param outfile (str): output file with users settings
    '''
    
    ### holder for user inputs for
    ### parameter setting
    users = {}
    
    ### print header
    print ('################# generate nufiles  #################')
    print ('#### This script generate a new nuparams textfile')
    print ('#### based on your following inputs.')
    print ('####')
    print ('#### NOTE: For plotting purposes, you can only')
    print ('####       change the values; Whether or not the')
    print ('####       the parameter is included as nuisance')
    print ('####       parameter or other minimizer / prior')
    print ('####       setting will be untouched.')
    print ('####')
    ### start interactive inputs
    for param in params:
        ## ask for user's input
        value = ask_param (param)
        ## update to `users` dictionary only
        ## if value is changed (keep it as str)
        if value: users [param] = value
        print ('####')
    print ('#### This is it.')
    print ('#### Thanks for your inputs :)')
    print ('#### Your inputs will be seen in the output file ...')
    print ('####      {0}'.format (outfile))
    print ('#####################################################')
    print ('####')
    return users

###########################################
#### function to split out file
###########################################
def generate_outfile (infile, outfile, users, do_print=False):

    ''' generate a nuparam textfile given users setting

        :param infile    (str): nuparam textfile template
        :param outfile   (str): output file with users settings
        :param users    (dict): new values of parameters based
                                on user's inputs
                                {param:value}
        :param do_print (bool): If True, print progress.
    '''
    
    ### open template nufile
    intxt = open (infile, 'r') 
    ### open output file
    outtxt = open (outfile, 'w')
         
    ### change parameters line by line
    for index, line in enumerate (intxt):
        if line=='' or line[0] in ['#', ' ', '\n'] or \
           line.strip ().split ()[0] not in users:
            ## rewrite any empty lines or
            ## when the line is split, its first element is not one of
            ## the parameters user wants to change
            outtxt.write (line)
        else:
            ## if the first element of the line is one of the
            ## parameters that user wants to change,
            ## then change the initial value in template to
            ## the user's value for the first three occurances
            param = line.strip ().split ()[0]
            initial = '2.526' if 'dm31' in param else \
                      str (int (nuParams (param).injected)) \
                      if 'holeice' in param else \
                      str (nuParams (param).injected)
            oline = line.replace (initial, users[param], 2).replace (users[param], initial, 1)
            ## print updated line if verbose
            createnu_print_line (param, oline, do_print=do_print)
            ## write to the output file
            outtxt.write(oline)

    ### close all files
    intxt.close()
    outtxt.close()
    return
        
###########################################
#### main function
###########################################
if __name__ == "__main__":
    
    ### parse user's options    
    usage = "%prog --outfile test.txt --verbose"
    parser = OptionParser(usage=usage)
    parser.add_option ("-o", "--outfile", type="string", default='~/test.txt',
                       help = "output filename with `.txt` extension")
    parser.add_option ("-v", "--verbose", action="store_true", default=False,
                       help = "print progress of this script as it runs")
    (options, args) = parser.parse_args()
    
    outfile = options.outfile
    verbose = options.verbose

    ### first check if the path for outfile exists
    if not os.path.exists (os.path.split (outfile)[0]+'/'):
        print ('Sorry.. your outfile path does not exist.')
        print ('Check that, and we will get something done then!')
        sys.exit ()
    
    ### ask for user's inputs on values
    ### of available parameters
    users = ask_inputs (outfile)

    ### report input and output files
    createnu_print_header (nufile, outfile, do_print=verbose)

    ### generate the nuparam textfile
    ### based on user's inputs
    generate_outfile (nufile, outfile, users, do_print=verbose)
    
    ### remind user where the outfile is if verbose
    createnu_print_ender (outfile, do_print=verbose)
