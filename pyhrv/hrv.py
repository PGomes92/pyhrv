# -*- coding: utf-8 -*-
"""
Heart Rate Variability
----------------------

This module contains the hrv() functions which allows a full HRV parameter computation using only a single line
function.

Notes
-----
..  This module is part of the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..	This package is a contribution to the open-source biosignal processing toolbox 'BioSppy':
	https://github.com/PIA-Group/BioSPPy

Author
------
..  Pedro Gomes, Master Student, University of Applied Sciences Hamburg

Thesis Supervisors
------------------
..  Hugo Silva, PhD, Instituto de Telecomunicacoes, PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

References (all submodules (will be updated soon))
--------------------------------------------------
..  Task Force of the European Society of Cardiology and the North American Society of Pacing and Electrophysiology,
	"Heart rate variability - Standards of measurement, physiological interpretation, and clinical use."
	European Heart Journal 17, pp. 354-381, 1996
	[Electrophysiology1996]
..  T.P. Hutchinson, "Statistics and graphs for heart rate variability: pNN50 or pNN20?" Physiological Measurement
	24, N9 - N14, 2003
	[Hutchinson2003]
.. 	J. Mietus, C. Peng, I. Henry, L. Goldsmith, and A. Goldberger, “The pNNx files: re-examining a widely used
	heart rate variability measure,” Heart, vol. 88, no. 4, pp. 378–380, 2002.
	[Mietus2002]
..	M. B. Tayel and E. I. Alsaba, “Poincaré Plot for Heart Rate Variability,”
	Int. J. Medical, Heal. Biomed. Bioeng. Pharm. Eng., vol. 9, no. 9, pp. 708–711, 2015.
	[Tayel2015]
..	F. Shaffer and J. P. Ginsberg, “An Overview of Heart Rate Variability Metrics and Norms,”
	Front. Public Heal., vol. 5, no. September, pp. 1–17, 2017.
	[Shaffer2017]

Last Update
-----------
13-10-2018

:copyright: (c) 2018 by Pedro Gomes (HAW Hamburg)
:license: BSD 3-clause, see LICENSE for more details.
"""
# Compatibility
from __future__ import absolute_import

# BioSppy import
import biosppy
from biosppy import utils
import matplotlib.pyplot as plt

# Import toolbox functions
import pyhrv
import pyhrv.tools as tools
from pyhrv.time_domain import time_domain
from pyhrv.frequency_domain import frequency_domain
from pyhrv.nonlinear import nonlinear


def hrv(signal=None,
		nn=None,
		rpeaks=None,
		sampling_rate=1000.,
		interval=[0, 10],
		plot_ecg=True,
		plot_tachogram=True,
		show=False,
		fbands=None,
		kwargs_ecg_plot={},
		kwargs_tachogram={},
		kwargs_time={},
		kwargs_nonlinear={},
		kwargs_welch=None,
		kwargs_lomb=None):
	"""Computes all HRV parameters of the HRV toolkit (see list below).

	Parameters
	----------
	signal : array_like
		ECG signal.
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	sampling_rate : int, float
		Sampling rate used for the ECG acquisition in (Hz).
	plot_ecg : bool, optional
		If True, plots ECG signal with specified interval ('signal' must not be None).
	plot_tachogram : bool, optional
		If True, plots tachogram with specified interval.
	fbands : dict, optional
		Dictionary with frequency bands for frequency domain (tuples or list).
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))
	show : bool, optional
		If true, shows all plots (default: True).
	kwargs_ecg_plot : dict, optional
		**kwargs for the plot_ecg() function (see 'tools.py' module)
	kwargs_tachogram={} : dict, optional
		**kwargs for the plot_tachogram() function (see 'tools.py' module)
	kwargs_time : dict
		**kwargs for the time_domain() function (see 'time_domain.py' module)
	kwargs_welch : dict, optional
		Dictionary containing the kwargs for the 'welch_psd' function (see docstring of the 'welch_psd' for more
		information)
	kwargs_lomb : dict, optional
		Dictionary containing the kwargs for the 'lomb_psd' function (see docstring of the 'lomb_psd' for more
		information)

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All time domain results.

	Returned Parameters - Time Domain
	---------------------------------
	..	NNI parameters (# of NNI, mean, min, max) (keys: 'nn_counter', 'nn_mean', 'nn_min', 'nn_max')
	..	NNI differences (mean, min, max, standard deviation) (keys: 'nn_diff_mean', 'nn_diff_min', 'nn_diff_max',
		'nn_diff_std')
	..	HR parameters (mean, min, max, standard deviation) (keys: 'hr_mean', 'hr_min', 'hr_max', 'hr_std')
	..	SDNN (key: 'sdnn')
	..	SDNN index (key: 'sdnn_index')
	..	SDANN (key: 'sdann')
	..	RMSSD (key: 'rmssd')
	..	SDSD (key: 'sdsd')
	..	nn50 & pNN50 (keys: 'nn50', 'pnn50')
	..	nnXX (XX = custom threshold) if specified (keys: 'nnXX', 'pnnXX')

	Returned Parameters - Frequency Domain
	--------------------------------------
	(below, X = one of the methods 'fft' or 'lomb')
	..	Peak frequencies of all frequency bands (key: 'X_peak')
	..	Absolute frequencies of all frequency bands (key: 'X_abs')
	..	Relative frequencies of all frequency bands (key: 'X_rel')
	..	Logarithmic frequencies of all frequency bands (key: 'X_log')
	..	Normalized frequencies of all frequency bands (key: 'X_norms')
	..	LF/HF ratio (key: 'X_ratio')
	..	Total power over all frequency bands (key: 'X_total')
	..	Interpolation method used for NNI interpolation (FFT/Welch's method only) (key: 'fft_interpolation')
	..	Resampling frequency used for NNI interpolation (FFT/Welch's method only) (key: 'fft_resampling_frequency')
	..	Spectral window used for PSD estimation of the Welch's method
		(key: 'fft_spectral_window)'

	Returned Parameters - Nonlinear
	-------------------------------
	..	SD1	(key: 'sd1')
	..	SD2 (key: 'sd2')
	..	SD2/SD1 (key: 'sd_ratio')
	..	Area of the fitted ellipse (key: 'ellipse_area')
	..	Sample Entropy (key: 'sampen')
	..	Detrended Fluctuations Analysis (short and long term fluctuations (key: 'dfa_short', 'dfa_long')

	Returned Figures
	----------------
	..	ECG plot (key: 'ecg_plot')
	..	Tachogram (key: 'tachogram_plot')
	..	Poincaré plot (key: 'poincare_plot')

	Notes
	-----
	..	Results are stored in a biosppy.utils.ReturnTuple object and need to be accessed with the respective keys as
		done with dictionaries (see list of parameters and keys above)
	..	Provide at least one type of input data (ecg_signal, nn, or rpeaks).
	..	Input data will be prioritized in the following order: 1. ecg_signal, 2. nn, 3. rpeaks.
	..	SDNN Index and SDANN: In some cases, the NN interval may start in a segment (or time interval) N and end only
		in the successive segment N+1. In this case, use the 'overlap' parameter to select if the first element of the
		segment should be dropped or not:
		..	If True: overlap allowed, returns all NNI but the cumulative sum of the NNI in a segment can be greater
			than the specified duration.
		..	If False: no overlap allowed, first NNI will be dropped and the cumulative sum of the NNI in a segment
			will always be < specified duration.

	Raises
	------
	TypeError
		If no input data for 'signal', 'nn' and 'rpeaks' provided.

	"""
	# Check input
	if signal is not None:
		signal, rpeaks = biosppy.signals.ecg.ecg(signal=signal, sampling_rate=sampling_rate, show=False)[1:3]
	elif nn is None and rpeaks is None:
		raise TypeError('No input data provided. Please specify input data.')

	nn = tools.check_input(nn, rpeaks)

	version = utils.ReturnTuple(('v.' + pyhrv.__version__, ), ('version', ))

	# Compute time domain results
	t_results = time_domain(nn=nn, show=False, **kwargs_time)

	# Compute frequency domain results
	f_results = frequency_domain(nn=nn, fbands=fbands, kwargs_welch=kwargs_welch, kwargs_lomb=kwargs_lomb, show=False)

	# Compute nonlinear parameters
	n_results = nonlinear(nn=nn, show=False, **kwargs_nonlinear)

	# Prepare output
	results = tools.join_tuples(t_results, f_results, n_results)

	# Plot ECG signal
	if plot_ecg and signal is not None:
		ecg_plot = tools.plot_ecg(signal, show=False, interval=interval, **kwargs_ecg_plot)
		results = tools.join_tuples(results, ecg_plot)

	# Plot Tachogram
	if plot_tachogram:
		tachogram_plot = tools.tachogram(nn=nn, show=False, interval=interval, **kwargs_tachogram)
		results = tools.join_tuples(results, tachogram_plot)

	if show:
		plt.show()

	# Output
	return results

if __name__ == '__main__':
	"""
	Example Script - Computing all HRV parameters of this toolbox u
	"""
	# Import
	import numpy as np

	# Load sample NNI series
	nni = np.load('./samples/series_1.npy')

	# Compute HRV results using all the default values
	hrv_results = hrv(nn=nni)

	# Print results to the console
	for key in hrv_results.keys():
		print(key)
		print("--> Result: %s" % str(hrv_results[key]))
