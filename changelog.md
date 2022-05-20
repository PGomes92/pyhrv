Update Version 0.4.1
--------------------
**Added**
- Request from Issue #22 & Issue #15: added ``mode`` to the Nonlinear Domain functions ``poincare()`` and ``dfa()``:
   - ``normal``: Returns computed output data and plot figures
   - ``dev``: Returns computed output data only

**Fixed**
  - Fixes Issue #28 where the pyhrv.utils._check_limits() function was incorrectly formatting strings for error and warning messages (and ironically causing an error due to this)
  - Fixed an issue where reports where not being generated correctly in Python 3.8+
  - Fixed an issue where the heart rate heatplot threw an "StopIteration" error in Python 3.8+
  - Fixes Issue #13 (missing import in ``nonlinear.py``)
  - Fixed an issue where some high level functions where using rpeak indices instead of temporal rpeaks locations when a raw ECG signal was provided

**Docs Update**
  - Fixes Issue #29 where multiple sample codes using biosppy where using the wrong r-peaks series for NNI computation
  - Fixes Issue #26 where docs where linking to outdated sample data

Update Version 0.4.0
--------------------
#### Major Changes
- updated hrv_keys.json and fixed some typos
- added ``modes`` to the Frequency Domain functions ``welch_psd()``, ``ar_psd()``, and ``lomb_psd()`:
   -	``normal``: Returns frequency domain parameters and PSD plot figure in a ReturnTuple object
	-	``dev``: Returns frequency domain parameters, frequency and power arrays, no plot figure
	-	``devplot``: Returns frequency domain parameters, frequency array, power array, and the plot figure
- Frequency Domain Comparison Features:
   - added ``pyhrv.frequency_domain.psd_comparison()`` function to plot multiple PSDs (welch, ar, or lomb) from multiple 
   NNI segments in single plot figure
    
        **Sample Plots:**
        [Welch](./SampleFigures/SamplePSDComparisonWelch.png),
        [Autoregressive](./SampleFigures/SamplePSDComparisonAR.png), and
        [Lomb-Scargle](./SampleFigures/SamplePSDComparisonLomb.png)     
        **Docs:**
        https://pyhrv.readthedocs.io/en/latest/_pages/api/frequency.html#d-psd-comparison-plot-psd-comparison
   - added ``pyhrv.frequency_domain.psd_waterfall()`` function to plot multiple PSDs (welch, ar, or lomb) from multiple 
   NNi segments in a 3D Plot (Sample plots: [Welch](./SampleFigures/SamplePSDWaterfallWelch.png), [Autoregressive](
   ./SampleFigures/SamplePSDWaterfallAR.png), [Lomb](./SampleFigures/SamplePSDWaterfallLomb.png))
- Added ``pyhrv.tools.heart_rate_heatplot()`` function for graphical visualization & classification of HR performance 
based on normal HR ranges by age and gender (Sample Plots: [plot 1](./SampleFigures/SampleHRHeatplot1.png), [plot 2](
./SampleFigures/SampleHRheatplot2), and [plot 3](./SampleFigures/SampleHRHeatplot3.png)).
- Added ``pyhrv.tools.radar_chart()`` function which plots a radar chart from a series of user-selected HRV parameters 
based on a
- removed ``overlap`` input argument from the ``pyhrv.utils.segmentation()`` functions as it had now effect anymore due to the previous updates of this function.
- ``overlap`` input argument has also been removed from the ``pyhrv.utils.sdnn_index()`` and ``pyhrv.utils.sdann()`` functions for the same reason (both use the ``pyhrv.utils.segmentation()``)
- Moved, fixed, and improved ``hrv_report()``
- Overall improved stability of the pyhrv package
- restructured pyhrv package (see especially ``tools.py`` and the ``utils.py`` modules):

```
    pyhrv                           # Toolbox
    ├── files                       # pyHRV support files
    |   ├── quickstart              # Figures for the quickstart guide (see pyhrv module README.md)
    |   ├── hr_heatplot.json        # HR reference normal values for the hr_heatplot() function
    |   ├── hrv_keys.json           # HRV keys to access the parameter results stored in
    |   ├── references.txt          # Publications references on which pyHRV functions are based
    |   ├── SampleECG.txt           # Sample ECG signal (for testing purposes)
    |   ├── SampleExport.json       # Sample HRV results export (for demonstration purposes)
    |   ├── SampleNNISeriesLong.npy # 60 minute sample NNI series (for testing purposes)
    |   ├── SampleNNISeriesShort.npy# 5 minute sample NNI series (for testing purposes)
    |   ├── SampleReport.csv        # Sample report in .csv format
    |   ├── SampleReport.pdf        # Sample report in .pdf format
    |   └── SampleReport.txt        # Sample report in .txt format
    |      
    ├── report                      # Subpackage for PDF, TXT, and CSV report generation
    |   ├── build                   # Default path for generated PDF, TXT, and CSV reports 
    |   ├── figure                  # Working directory where figures for the PDF report are temporarily stored
    |   ├── templates               # LaTeX templates for the PDF reports (feel free to customize)
    |   ├── __init__.py             # Report init file
    |   ├── main.tex                # Main LaTeX file
    |   ├── parameters.tex          # File containing all the variables in which the HRV parameters will be stored
    |   ├── pyhrv.png               # pyHRV logo for the PDF header
    |   ├── README.md               # Report README
    |   └── settings.tex            # Settings file for the LaTeX-based project
    |   
    ├── __init__.py                 # pyHRV init file
    ├── __version__.py              # pyHRV version
    ├── hrv.py                      # HRV function
    ├── frequency_domain_.py        # Frequency Domain functions
    ├── nonlinear.py                # Nonlinear Parameters functions
    ├── time_domain.py              # Time Domain functions
    ├── tools.py                    # HRV tools functions (Tachogram, ECG, radar chart, HR heatplot, and other 
    |                               # Comparison functions
    └── utils.py                    # General data verification, conversion, file handling, etc.

```

- the ``tools.py`` now only contains tools specifically designed for HRV analysis (support)
- the ``utils.py`` now only contains support functions required for the computation of the tools and computation of 
HRV parameters

#### Minor Changes
- pyhrv.nonlinear.dfa(): added ``dfa_alpha1_beats`` return key with short term fluctuation and ``dfa_alpha2_beats` return key with long term fluctuation values
- updated hrv_keys.json file
- added ``complete`` as valid value for the ``interval`` input parameter of the ``tools.plot_ecg`` and ``tools
.tachogram``
functions
- pyhrv.time_domain.tinn() & pyhrv.time_domain.geometrical_parameters() functions no issue a warning due to the current issue of the TINN function of providing wrong results (see issue #5)
- removed 'overlap' input parameter from ``sdnn`` and ``sdnn_index`` functions 
- updated main README
- Fixes #6, #7 and #9

Update Version 0.3.2
--------------------
- IMPORTANT: fixes installation issues under Python3 using pip install (full support must yet be tested)
- prepared mode-dependent returns for functions with plotting features

Update Version 0.3.1
--------------------
- fixes index out of range error in tools.segmentation() function
- fixes BioSPPy import in hrv.py

Update Version 0.3
------------------
- added ReadTheDocs documentation
- added new frequency domain method (Autoregressive)
- fixed missing histogram figure visualization under macOS
- fixed 'pyhrv.tools.check_input()' function where input parameter was changed during the computation
- updated docstrings and added documentation
- added references
- improved hrv() function
- improved nonlinear() function

Update Version 0.2
------------------
- changed toolbox name to 'pyhrv'
- added install process using `pip`
- added frequency domain parameters and method (Welch's method and Lomb-Scargle)
- added `frequency_domain()` module level function
- added frequency parameter computation in `hrv()` package level function
- rearranged file structure in repository
- added `references.txt` file
- improved `hrv_report()` function
- updated `hrv_keys().json` file
- minor bug fixes
- added 50 sample NNI series extracted from the [MIT-BIH NSRDB](https://physionet.org/physiobank/database/nsrdb/) database
- updated sample HRV report and exports
- fixed bug in `pyhrv.tools.segmentation()` function where NNI overlapping from one segment to the other got dropped

Update Version 0.1.2
--------------------
- added new time domain HRV parameters -> geometrical parameters

  * Triangular Index:  _hrv.time_domain.triangular_index()_
  * TINN:  _hrv.time_domain.tinn()_
  * Histogram (helper function powered by NumPy and Matplotlib):  _hrv.time_domain._get_histogram()
  * Geometrical Parameters (calls _tinn()_ and _triangular_index()_ as those parameters are rarely used alone): _hrv.time_domain.geometrical_parameters()_

- updated _hrv.time_domain.time_domain()_ to also return the new geometrical parameters
- added new nonlinear parameters -> sample entropy and detrended fluctuation analysis (DFA)

  * Sample Entropy:  _hrv.nonlinear_parameters.sample_entrop()_
  * DFA: _hrv.nonlinear_parameters.dfa()_
  * Both functions powered by the [nolds](https://github.com/CSchoel/nolds)

- updated _hrv.nonlinear_parameters.nonlinear()_ to also return the new nonlinear parameters
- added new tools function -> HRV Report
  * Entirely new HRV report generator:  _hrv.tools.hrv_report()_
  * generates reports in .txt or .csv format

- added_ _check_fname()_ function to detect existing files for the generated HRV exports and HRV reports and automatically generate a new file name to avoid overwriting existing files (e.g. if 'Sample.txt' exists, it will 'increment' the new file name to 'Sample_1.txt' to avoid overwriting 'Sample.txt')

- added 'SampleReport.txt' and 'SampleReport.csv' to the './hrv/files/'
- updated 'SampleExport.json' in './hrv/files/'
- fixed minor fixed in some tools functions
- improved compatibility with Python 3
- updated example sections in each module

Update Version 0.1.1
--------------------
- restructured entire repository
- tested functionality and validated results for the _time_domain_, _nonlinear_parameters_, and _tools_ modules
- fixed example scripts at the end of each module
- added new functions to the _tools_ module: _hrv_import()_,
_hrv_export_(), _check_input()_
- fixed bugs in the _tools.segmentation()_ function

Version Version 0.1
-------------------
- time domain functions
- poincaré nonlinear functions
- basic tools for HRV functions