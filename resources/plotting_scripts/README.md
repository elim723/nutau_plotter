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

The actual plot commands are done by `histogram2D` class defined in `plotter/Histogram2D.py`. All the plot stylings / formattings are passed to this class as keyword dictionaries. Multiple styling dictionaries are defined in `plot_greco2D.py`. For example, `settings` variable contains many dictionaries, each of which define specific stylings for one variable; these stylings include the cut values / the x, y, z (colorbar) ranges / etc. Global stylings, such as grid lines / fontsizes / cut lines / etc, are also defined in `plot_greco2D.py`. Those input stylings are then passed to the `histogram2D` instance when it is created. Then, within the class, the stylings are applied to the plotting commands accordingly.

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

The actual plot commands are done by `histogram2D` class defined in `plotter/Histogram2D.py`. All the plot stylings / formattings are passed to this class as keyword dictionaries. Multiple styling dictionaries are defined in `plot_dragon2D.py`. For example, `settings` variable contains many dictionaries, each of which define specific stylings for one variable; these stylings include the cut values / the x, y, z (colorbar) ranges / etc. Global stylings, such as grid lines / fontsizes / cut lines / etc, are also defined in `plot_dragon2D.py`. Those input stylings are then passed to the `histogram2D` instance when it is created. Then, within the class, the stylings are applied to the plotting commands accordingly.

plot_resolutions.py
-------------------

This script generates pegleg resolution plots:

- Figure 5: Energy resolution from both GRECO and DRAGON
- Figure 6: Zenith resolution from both GRECO and DRAGON

In fact, this script generates six pdfs, two of which are energy resolutions and the remaining four are zenith resolutions. Each pdf has four subplots, corresponding to the four neutrino flavor and interaction types: numu CC, nue CC, nutau CC, and nu NC. Only two out of the six pdfs will eventually be put into the paper: one for energy energy and the other for zenith resolution. The reason why this script produces six resolution plots is that the x/y axes are currently (10/02/2018) under discussion. In particular, energy resolutions can be presented in linear vs linear or log10 vs log 10, where x and y axes are true and reconstructed energy respectively. On the other hand, zenith resolution plots show the difference between reconstructed and true zenith values as a function of true energy. Therefore, a total of four possible axis combinations can be drawn: 1. (reco - true) cos zenith vs linear true energy, 2. (reco - true) cos zenith vs log10 true energy, 3. (reco - true) zenith in degrees vs linear true energy, and 4. (reco - true) zenith in degrees vs log10 true energy. Since no conclusion of which plots to be shown was drawn, this script plots all four possibilities. Once the decision of which plots to be shown is made, I recommand focusing on the styles of the decided axis combination and turn off the others.

If resolution values are already calculated and stored in a pickled dictionary (such as `/data/user/elims/nutau_pickled_data/pegleg_resolutions.p`), then this script can be run via the command line:
```
python plot_resolutions.py --outdir $outdir --infile /data/user/elims/nutau_pickled_data/pegleg_resolutions.p
```
The script does two things. First, it define the plot stylings and formattings dictionaries. For each of the 6x4 subplots, one can set the text / legend locations and the axis tick locations and labels. The fact that this script formats 6x4 plots make the format dictionaries look complicated. Once it was decided what axis combinations to show, I recommand just focusing on those needed. The second thing this script does is to read the input file `infile` and pass the data and defined plot stylings into the `Resolution` class in `plotter/resolution.py`. The `Resolution` instances once it is created will then pass the plot styles and formats to the plot commands accordingly.



