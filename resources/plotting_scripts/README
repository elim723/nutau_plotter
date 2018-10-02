This folder contains the plotting scripts to generate plots in the nutau paper. The current latest version (10/02/2018) is available is located in:

https://wiki.icecube.wisc.edu/images/f/fa/Nutau_v2.0.pdf

The following listed out which figure in the pdf is made by which python script and how. Note that all plotting scripts require one input argument `$outdir`, which is the address where the output pdfs are dumped to.

plot_greco1D.py
---------------

This script generates 1D histograms from the GRECO sample, which includes

- Figure 2 : GRECO BDT distribution at L4
- Figure 3 : GRECO BDT distribution at L5
- Figure 19: GRECO ToI distribution at L4
- Figure 20: GRECO VICH distribution at L5
- Figure 21: GRECO Fill-Ratio distribution at L6
- Figure 22: GRECO NChannel distribution at L6

To run this script:
```
$ python plot_greco1D.py --outdir $outdir (--fitnorm)
```

The script itself does two things. First, it defines all the plot stylings in forms of dictionaries. Second, it extracts the information from GRECO datasets located at `/data/user/elims/nutau_pickled_data/greco/preweighted/`. Note that, when plotting additional 1D histograms, make sure the variables are in those pickled data.

The actual plot commands are done by `histogram1D` class defined in `plotter/Histogram1D.py`. All the plot stylings / formattings are passed to this class as dictionaries. Multiple styling dictionaries are defined in `plot_greco1D.py`. For example, `settings` variable contains many dictionaries, each of which define specific stylings for one variable; these stylings include the cut values / lcoation of legends / the x and y ranges for the distribution and ratio plots / etc. Global stylings, such as grid lines / fontsizes / ratio reference lines / cut lines / etc, are also defined in `plot_greco1D.py`. Those input stylings are then passed to the `histogram1D` instance when it is created. Then, within the class, the stylings are applied to the plotting commands accordingly.

The `fitnorm` is an optional argument. If not provided, each histogram from a data type (numu / nue / nutau / noise / muon) has a normalization of 1.0. If `fitnorm` is turned on, the normalizations are fitted to best match the data points. This `fitnorm` feature is added to better match MC and data at lower selection levels. The fitting is done by the `fitter1D` class in `plotter/Fitter1D.py` via scipy.minimizer.

plot_greco2D.py
---------------

This script generates 2D histograms from the GRECO sample, which includes

- Figure 23: GRECO Muon Fraction from FiniteReco Z vs Rho at L6
- Figure 24: GRECO Muon Fraction from PegLeg Z vs Rho at L7

To run this script:
```
$ python plot_greco2D.py --outdir $outdir
```

The script itself does two things. First, it defines the plot stylings in forms of dictionaries. Stylings specfic for each 2D plot are defined by the `settings` variable, whereas global stylings are defined by `format_xxxx` dictionaries. Note that, each 2D plot has a set of cut lines defined. If the x/y ranges need to be changed, the definitions of cut lines may need to be adjusted. The second thing the script does is to extract the information from GRECO datasets located at `/data/user/elims/nutau_pickled_data/greco/preweighted/`. When plotting additional 2D histograms, make sure the variables are in those pickled data.

plot_dragon1D.py
---------------

This script generates 1D histograms from the DRAGON sample, which includes

- Figure 26: DRAGON CoG separation distribution at L5
- Figure 27: DRAGON BDT distribution at L5
- Figure 29: DRAGON track stopping Z position distribution at L6

To run this script:
```
$ python plot_dragon1D.py --outdir $outdir (--fitnorm)
```

The script itself does two things. First, it defines all the plot stylings in forms of dictionaries. Second, it extracts the information from DRAGON datasets located at `/data/user/elims/nutau_pickled_data/dragon/preweighted/`. Since no pre-Level 4 hdf5 files available for real data and Elim does not have the time to pickled them, pre-Level 4 distribution plots are not available. One must go through all i3 files to pickle the events needed (see `resources/pickler/`). Note that, when plotting additional 1D histograms, make sure the variables are in those pickled data. 

The actual plot commands are done by `histogram1D` class defined in `plotter/Histogram1D.py`. All the plot stylings / formattings are passed to this class as dictionaries. Multiple styling dictionaries are defined in `plot_dragon1D.py`. For example, `settings` variable contains many dictionaries, each of which define specific stylings for one variable; these stylings include the cut values / lcoation of legends / the x and y ranges for the distribution and ratio plots / etc. Global stylings, such as grid lines / fontsizes / ratio reference lines / cut lines / etc, are also defined in `plot_dragon1D.py`. Those input stylings are then passed to the `histogram1D` instance when it is created. Then, within the class, the stylings are applied to the plotting commands accordingly.

The `fitnorm` is an optional argument. If not provided, each histogram from a data type (numu / nue / nutau / muon) has a normalization of 1.0. If `fitnorm` is turned on, the normalizations are fitted to best match the data points. This `fitnorm` feature is added to better match MC and data at lower selection levels. The fitting is done by the `fitter1D` class in `plotter/Fitter1D.py` via scipy.minimizer.

plot_dragon2D.py
---------------

This script generates 2D histograms from the DRAGON sample, which includes

- Figure 28: DRAGON Muon Fraction from track starting Z vs Rho at L6

To run this script:
```
$ python plot_dragon2D.py --outdir $outdir
```

The script itself does two things. First, it defines the plot stylings in forms of dictionaries. Stylings specfic for each 2D plot are defined by the `settings` variable, whereas global stylings are defined by `format_xxxx` dictionaries. Note that, each 2D plot has a set of cut lines defined. If the x/y ranges need to be changed, the definitions of cut lines may need to be adjusted. The second thing the script does is to extract the information from DRAGON datasets located at `/data/user/elims/nutau_pickled_data/dragon/preweighted/`. When plotting additional 2D histograms, make sure the variables are in those pickled data.