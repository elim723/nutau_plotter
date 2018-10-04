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
$ python plot_resolutions.py --outdir $outdir --infile /data/user/elims/nutau_pickled_data/pegleg_resolutions.p
```
The script does two things. First, it defines the plot stylings and formattings dictionaries. For each of the 6x4 subplots, one can set the text / legend locations and the axis tick locations and labels. The fact that this script formats 6x4 plots make the format dictionaries look complicated. Once it was decided what axis combinations to show, I recommand just focusing on those needed. The second thing this script does is to read the input file `infile` and passes the data and defined plot stylings into the `Resolution` class in `plotter/resolution.py`. The `Resolution` instances once it is created will then pass the plot styles and formats to the plot commands accordingly.

If no input file `infile`, one can run the same script via
```
#### first set up the icecube software 
$ bash /data/user/elims/nutau_software/debugv3/env-shell.sh
#### then run the plotting script, which will calculate
#### the new resolutions needed
$ python plot_resolutions.py --outdir $outdir 
```
Here, the script performs the calculations of resolution values for each neutrino flavor and interaction type. First, it gets all neutrino (both CC and NC) events from both GRECO and DRAGON samples at their final-ist levels within their analysis histogram ranges (see `get_events ()` in `plot_resolutions.py`). In doing so, a weight is calculated for each event by `get_weights ()`; this step assumes baseline values for all systematic values documented in `events_collecter/nufiles/nuparams_template.txt`. This is the part where IceCube software (NuFlux and prob3) is needed. Once the weights are calculated, NC weights are concatenated as one category (see `merge_nc ()`). When all events are collected, resolutions are calculated for each of the six axes combinations by `get_all_resolutions ()`. This involves a double for-loop that calculates the resolution of a given axes combination for each data type (numu CC, nue CC, nutau CC, nuNC) and for each sample (GRECO and DRAGON). For a given configuration, its resolution as a function of true energy is calculated in `get_resolutions ()`. Here, the x values (either linear or log10 true energy), y values (either linear or log10 reconstructed energy or difference between reconstructed and true zenith in cosine or degrees), and the weights are defined from a set of events (of a specific data type from a sample). These three arrays, as well as the pre-defined bin edges for the x axis, are passed to `get_statistics ()`. This is where the median and +/- 1sigma values are defined by selecting the events in the pre-defined x bins and by building a CDF in `statistics ()`. Once the values at lower 34%, median 50%, and upper 34% are found, these values are saved into a dictionary. These values from all possible configurations are then stored in a big dictionary and saved in `$outdir/pegleg_resolutions.p`, such that user can reuse the output file as the input file next time when running the plotting script.

plot_classification.py
----------------------

This script generates PID plot:

- Figure 7: Fraction of track-like events as a function of energy from both GRECO and DRAGON samples.

The output plot has two subplots: left from analysis A (GRECO) and right from analysis B (DRAGON). For each subplot, four lines are plotted, each of which represents the track fractions as a function of true energy for numu CC, nue CC, nutau CC, or nu NC. 

If track fractions are already calculated and stored in a pickled dictionary (such as `/data/user/elims/nutau_pickled_data/pid_track_fractions.p`), then this script can be run via the command line:
```
$ python plot_classification.py --outdir $outdir --infile /data/user/elims/nutau_pickled_data/pid_track_fractions.p
```
The script does two things. First, it defines the plot stylings and formattings dictionaries. For each of the two subplots, one can set the legend location and the axis tick locations and labels. Second, this script reads the input file `infile` and passes the data and defined plot stylings into the `Resolution` class in `plotter/resolution.py`. The `Resolution` instance, once it is created, will then pass the plot styles and formats to the plot commands accordingly.

If no input file `infile` (e.g. one wants to re-calculate the track fractions), then this script is executed without the `infile` option; i.e.
```
#### first set up the icecube software 
$ bash /data/user/elims/nutau_software/debugv3/env-shell.sh
#### then run the plotting script, which will calculate
#### the new track fractions needed
$ python plot_classification.py --outdir $outdir
```
This script then does the calculations of resolution values for each neutrino flavor and interaction type. First, it gets all neutrino (both CC and NC) events from both GRECO and DRAGON samples at their final-ist levels within their analysis histogram ranges (see `get_events ()` in `plot_resolutions.py`). In doing so, a weight is calculated for each event by `get_weights ()`. This step assumes baseline values for all systematic values stated in `events_collecter/nufiles/nuparams_template.txt`; one would need to define a new `nufile` if he / she wants to use other systematic values. This is also the part where IceCube software (NuFlux and prob3) is needed. Once the weights are calculated, NC weights are concatenated as one category by `merge_nc ()`. When all events are collected, resolutions are calculated for each of the six axes combinations by `get_all_resolutions ()`. This involves a double for-loop that calculates the resolution of a given axes combination for each data type (numu CC, nue CC, nutau CC, nuNC) and for each sample (GRECO and DRAGON). For a given configuration, its resolution as a function of true energy is calculated in `get_resolutions ()`. Here, the x values (either linear or log10 true energy), y values (either linear or log10 reconstructed energy or difference between reconstructed and true zenith in cosine or degrees), and the weights are defined from a set of events (of a specific data type from a sample). These three arrays, as well as the pre-defined bin edges for the x axis `loge_edges`, are passed to `get_statistics ()`. This is where the median and +/- 1sigma values are defined by selecting the events in each pre-defined x bin and by building a CDF in `statistics ()` for that x bin. Once the values at lower 34%, median 50%, and upper 34% are found, these values are saved into a dictionary. These values from all x bins from possible configurations are then stored in a big dictionary and saved in `$outdir/pegleg_resolutions.p`. Next time of plotting, user can reuse the output file as the input file `infile` when running this plot_classification.py script.

plot_greco_numu.py
------------------

This script generates 90% contour plot:

- Figure 13: 90% contour from standard atmospheric oscillation parameters via numu disapperance.

To generate the contour plots, run
```
$ python plot_greco_numu.py --outdir $outdir
```
This scipt does three things. First, it defines the plot styles via dictionaries. These dictionaries includes `format_global` for fontsizes, and grid styles, `format_marker` for best fit point, `format_legend` for its locations / number of columns / fontsize, and `format_axes` for the ticks and ranges. Note that the color and line style for each experiment is defined in `plotter/defaults.py`. Second, it reads the data points measured by different experiments; Elim has put them into a pickled dictionary in `/data/user/elims/nutau_pickled_data/numu_contours.p`. This dictionary includes the result from GRECO sample, as well as the current latest results from T2K (2017), NOvA (2017), MINOS (2016), SK (2017), and IceCube (2017) from DRAGON sample. Last, the 90% contour points and plot styles are passed to the `Contour` class defined in `plotter/contour.py` for plotting.

plot_sig_to_bg.py
-----------------

This script generates the signal to sqrt background plot:

- Figure 8: Signal nutau (CC + NC) divided by the sqrt-root of background (not nutau) events per bin.

To generate the signal to sqrt (background), run
```
$ python plot_sig_to_bg.py --outdir $outdir
```
This script does two things. First, the plot style parameters are stored in multiple style dictionaries. Next, for a given sample (hardcoded to be GRECO), reconstructed observables and weights for all events at the final-ist selection level from all data types are loaded. The script then passes both events and style dictionaries to the `Histogram3D` class defined in `plotter/histogram3D.py` for plotting.

A special feature `stretch_cscd` is available for the GRECO signal to sqrt background plot. In GRECO analyses, the cascade and track histograms have different number of bins. For cascade, there are 8 energy and 10 zenith bins (i.e. 80 bins), whereas the track histogram have only 60 bins (6 energy and 10 zenith). The energy bins below 10 GeV in the track histogram are dropped because track-like events in GRECO analyses are defined as events with track lengths longer than 50 meters. Due to the minimum ionizing muon hypothesis in the PegLeg reconstruction, such a track-like muon with a track length longer than 50 meters must have a minimum energy of 10 GeV. Therefore, those energy bins in the track histogram are not included in the GRECO analyses. Because of the dropped bins, two styles are available for GRECO-only muon histograms by turning on/off the `stretch_cscd` style-key in the `settings` style-dictionary in the script. When `stretch_cscd` is set to be False, the GRECO muon histograms have the same x scales as the one in DRAGON but simply not display the first two energy bins in the track histogram. When `stretch_cscd` is True, the cascade-to-track width ratio is 28:21, meaning that the cascade histogram is stretched to fill up the empty space. One could easily make `stretch_cscd` an option in the command line if needed.

plot_muon_histograms.py
-----------------------

This script generates the muon histograms from both samples:

- Figure 12: Event distributions of the atmospheric muon background.

This scripts can be run in three ways depending on which sample(s) the user want to include. The current option includes both GRECO (top orange) and DRAGON (bottom red) samples using the following command.
```
$ python plot_muon_histograms.py --outdir $outdir
```
If one wants to include just DRAGON sample, add a flag like this.
```
$ python plot_muon_histograms.py --dragon --outdir $outdir
```
This will produce the DRAGON's muon histograms, given 2.5 years of effective livetime, in both cascade (left) and track (right) histograms, each of which has 8 x 8 = 64 bins. Similarly, the muon histograms for GRECO analyses, with an effective livetime of 2.27 years, are done via a `--greco` flag.
```
$ python plot_muon_histograms.py --greco --outdir $outdir
```
The `stretch_cscd` option explained in `plot_sig_to_bg.py` is also available for GRECO-only muon histogram plot.

This script does several things. First, the plot style parameters are stored in multiple style dictionaries. Next, the script collects events from both samples, including setting the same names for the different PID variables. Note that, since DRAGON muon histogram is obtained from inverted background, all events are weighted equally. To simplify the weight calculation, Elim simply rescaled event weights such that the total muon rate at Level 6 (LE) from the giant table in the appendix. The number per bin matches the one produced by FeiFei/Philip. Then, the script passes the events and style dictionaries into histogram classes based on how many sample user wants to include. If either `--dragon` or `--greco` flag is turned on, all data is passed to the `Histogram3D` class defined in `plotter/histogram3D.py` (see `plot_histo1x2 ()`). If none of the flags are given, data is passed to the `Histogram6D` class defined in `plotter/histogram6D.py` (see `plot_histo2x2 ()`), which has a different canvas layout as the `Histogram3D` class.

Note that, if one wants to change the weights in either sample, a new event dictionary must be created. Please follow the instructions in `resources/events_collecter/` to make a new event dictionary.

plot_sys_effects.py
-------------------
The last available script is to produce the systematic effects on histograms.

- Figure 10: Effect of selected systematic uncertainties on nominal event distribution.

To run this script, an input file with all weights pre-calculated is required. Such a file was prepared and located at `/data/user/elims/nutau_pickled_data/greco/syseffects_events.p`. 
```
$python plot_sys_effects.py --outdir $outdir/ --infile /data/user/elims/nutau_pickled_data/greco/syseffects_events.p
```
This script does two things. First, it defines all style-dictionaries. If GRECO sample is given, the `stretch_cscd` option explained in `plot_sig_to_bg.py` is also available for GRECO-only systematic effect plots. Then, the script loads the data from the input file and passes all the information to the `HistogramND` class defined in `plotter/histogramND.py`.

The prepared input file mentioned here are made specifically for the one shown in the current draft; nue_numu_ratio + 1 sigma, up_horizontal_ratio + 1 sigma, hole-ice forward + 1, delta m^2_31 set at 2.52 m eV^2, and axial mass resonance + 1 sigma. If one wants to change those values and/or the parameters included, a new input file can be created by `resources/events_collecter/collect_sysevents.py`. Please refer to the README in `resources/events_collecter/` and the comments in that python script for more instruction. One may also want to change the `deltas` variable in `plotter/defaults.py` where the values of the systematic parameters for the systematic effect plot are defined.