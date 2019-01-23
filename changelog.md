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