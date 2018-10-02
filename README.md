Welcome! This project is specifically designed to make plots that are used in the IceCube nutau paper. The plotting tools are located in `resources/plotting_scripts/`. Please read through the READMEs to see how the software is used for its purposes.

Overall Goal
============

This plotter contains plotting scripts, as well as the tools for pickling data, calculating weights, building a hyperplane, etc. *No data is attached to this tool.* Data and specific software tools for weight calculations (e.g. prob3, NuFlux, and weighting) can be found in IceCube Wisconsin machines.

Getting Started
---------------

1. You must have an IceCube account to get the data and run these codes. 
2. All plotting scripts can be run using only Python > 2.7. So, they can be run without IceCube software.
3. To install this project, one should first clone it to a wisconsin machine.
> $ mkdir nutau_plotter

> $ cd nutau_plotter/

> $ git clone https://github.com/elim723/nutau_plotter.git .

Then, install the pacakge locally.

> $ python setup.py install --user

Running Tests
-------------

Given the time limit I have, I did not write a test for this tool. However, one can test if the tool works via the followings.
> $ cd nutau_plotter/resources/plotting_scripts/

> $ python plot_dragon1D.py --outdir plots/

Now you should see a few 1D distribution plots from Analysis B (DRAGON sample) in `plots/`.

Getting Deeper
--------------

Most plotting scripts only require Python because I have prepared the input files with some pre-calculated values (weights or track fractions or resolutions). When those values need to be re-calculated, extra software is needed (e.g. prob3, NuFlux, weighting). I have provided an IceCube enviornment for all the tools needed just in case. When needed, one should set up the enviornment as below.
> $ eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/setup.sh`

> $ bash /data/user/elims/nutau_software/debugv3/env-shell.sh

Please note that these are the latest software at the time (10/02/2018). Updated software may create conflicts.

Final Notes
===========

I do not have time to write a webpage-style documentations. But please do read through each README and comments in scripts to understand the purpose of each script, how to run it, and how to change the plotting styles.