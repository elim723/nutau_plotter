Welcome! This project is specifically designed to make plots that are used in the IceCube nutau paper. The plotting tools are located in `resources/plotting_scripts/`. Please read through the READMEs to see how the software is used for its purposes.

Overall Goal
============

This plotter contains plotting scripts, as well as the tools for pickling data, calculating weights, building a hyperplane, etc. *No data is attached to this tool.* Data and specific software tools for weight calculations (e.g. prob3, NuFlux, and weighting) can be found in IceCube Wisconsin machines.

Getting Started
---------------

1. You must have an IceCube account to get the data and run these codes. So, log into one of the cobalt machines.
2. All plotting scripts can be run using only Python. Please provide a clean enviornment (remove any `cvmfs` lines in your `/.bashrc`) and start with `py2-v3.0.1`, the latest stable `cvmfs` available on 10/02/2018.
```
$ eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/setup.sh`
```
3. To install this project, one should first clone it.
```
$ mkdir nutau_plotter
$ cd nutau_plotter/
$ git clone https://github.com/elim723/nutau_plotter.git .
```
Then, install the pacakge locally.
```
$ python setup.py install --user
```
Running Tests
-------------

Given the time limit I have, I did not write a test for this tool. However, one can test if the tool works via the followings.
```
$ cd nutau_plotter/resources/plotting_scripts/
$ python plot_dragon1D.py --outdir plots/
```
Now you should see a few 1D distribution plots from Analysis B (DRAGON sample) in `plots/`.

Once the above is working, you can make all the plots provided by this code via
```
$ cd nutau_plotter/resources/plotting_scripts/
$ bash run_all.sh plots/
```

Plotting Scripts Overview
-------------------------

If you are just changing the plotting styles, `resources/plotting_scripts/` are the plotting scripts. I tried to make the scripts as generalized as possible, such that user can change any plotting options via dictionaries. For example, in `plot_dragon1D.py`, there is a `format_global` variable:

```
format_global = {'plottotalmc'          :False           ,
                 'hist_ylabel'          :r'Event Rate Hz',
                 'hist_logy'            :True            ,
                 'hist_legend_fontsize' :13              ,
                 'hist_legend_alpha'    :1.0             ,
                 'ratio_ylabel'         :r'Ratio to MC'  ,
                 'ratio_logy'           :False           ,
                 'ratio_legend_fontsize':11              ,
                 'ratio_legend_alpha'   :1.0             ,
                 'tick_fontsize'        :11              ,
                 'label_fontsize'       :13              ,
                 'grid_alpha'           :0.2             ,
                 'grid_linestyle'       :'-'             ,
                 'grid_linewidth'       :0.5             }
```

Explanations to the style-keys can be found in the script, and users can simply change the style-values to change the plotting styles.

If you cannot find the style-key you want, you may look into the relevant classes in `plotter/`. For example, if you want to change the canvas settings of plots made by `plot_dragon1D.py`, look into the `__init__ ()` in `plotter/histogram1D.py`. Note that, stylings in those plotter classes are global, meaning that it will also affect the plots from e.g. `plot_greco1D.py`.

Getting Deeper
--------------

Most plotting scripts only require Python because I have prepared the input files with some pre-calculated values (weights or track fractions or resolutions). When those values need to be re-calculated, extra software is needed (e.g. prob3, NuFlux, weighting). I have provided an IceCube enviornment for all the tools needed just in case. When needed, one should set up the enviornment as below.
```
$ eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/setup.sh`
$ bash /data/user/elims/nutau_software/debugv3/env-shell.sh
```
Please note that these are the latest software at the time (10/02/2018). Updated software may create conflicts.

Final Notes
===========

I do not have time to write a webpage-style documentations. But please do read through each README and comments in scripts to understand the purpose of each script, how to run it, and how to change the plotting styles.