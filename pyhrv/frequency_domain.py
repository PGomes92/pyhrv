# -*- coding: utf-8 -*
"""
pyHRV - Frequency Domain Module
--------------------------------

This module provides function to compute frequency domain HRV parameters using R-peak locations
or NN intervals extracted from an ECG lead I-like signal. The implemented Power Spectral Estimation (PSD)
estimation methods are:

	* Welch's Method
	* Lomb-Scargle
	* Autoregressive

Notes
-----
..  This module is part of the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..	This module is a contribution to the open-source biosignal processing toolbox 'BioSPPy':
	https://github.com/PIA-Group/BioSPPy
..	You find the API reference for this module here:
	https://pyhrv.readthedocs.io/en/latest/_pages/api/frequency.html
.. 	See 'references.txt' for a full detailed list of references

Author
------
..  Pedro Gomes, Master Student, University of Applied Sciences Hamburg

Thesis Supervisors
------------------
..  Hugo Silva, PhD, Instituto de Telecomunicacoes & PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

Last Update
-----------
19-11-2018

:copyright: (c) 2018 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.
"""
# Compatibility
from __future__ import absolute_import, division, print_function

# Imports
import warnings
import spectrum
import numpy as np
import scipy as sp
import matplotlib as mpl
from scipy.signal import welch, lombscargle
from matplotlib import pyplot as plt

# biosppy imports
import biosppy
from biosppy import utils

# Local imports/HRV toolbox imports
import pyhrv.tools as tools

# Surpress Lapack bug 0038 warning from scipy (may occur with older versions of the packages above)
warnings.filterwarnings(action="ignore", module="scipy")


def welch_psd(nni=None,
			  rpeaks=None,
			  fbands=None,
			  nfft=2**12,
			  detrend=True,
			  window='hamming',
			  show=True,
			  show_param=True,
			  legend=True):
	"""Computes a Power Spectral Density (PSD) estimation from the NNI series using the Welch’s method
	and computes all frequency domain parameters from this PSD according to the specified frequency bands.

	References: [Electrophysiology1996], [Umberto2017], [Welch2017]
	Docs:		https://pyhrv.readthedocs.io/en/latest/_pages/api/frequency.html#welch-s-method-welch-psd

	Parameters
	----------
	rpeaks : array
		R-peak locations in [ms] or [s]
	nni : array
		NN-Intervals in [ms] or [s]
	fbands : dict, optional
		Dictionary with frequency bands (2-element tuples or list)
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.000Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))
	nfft : int, optional
		Number of points computed for the FFT result (default: 2**12)
	detrend : bool optional
		If True, detrend NNI series by subtracting the mean NNI (default: True)
	window : scipy window function, optional
		Window function used for PSD estimation (default: 'hamming')
	show : bool, optional
		If true, show PSD plot (default: True)
	show_param : bool, optional
		If true, list all computed PSD parameters next to the plot (default: True)
	legend : bool, optional
		If true, add a legend with frequency bands to the plot (default: True)

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All results of the Welch's method's PSD estimation (see list and keys below)

	Returned Parameters & Keys
	--------------------------
	..	Peak frequencies of all frequency bands in [Hz] (key: 'fft_peak')
	..	Absolute powers of all frequency bands in [ms^2][(key: 'fft_abs')
	..	Relative powers of all frequency bands [%] (key: 'fft_rel')
	..	Logarithmic powers of all frequency bands [-] (key: 'fft_log')
	..	Normalized powers of all frequency bands [-] (key: 'fft_norms')
	..	LF/HF ratio [-] (key: 'fft_ratio')
	..	Total power over all frequency bands in [ms^2] (key: 'fft_total')
	..	Interpolation method used for NNI interpolation (key: 'fft_interpolation')
	..	Resampling frequency used for NNI interpolation (key: 'fft_resampling_frequency')
	..	Spectral window used for PSD estimation of the Welch's method (key: 'fft_spectral_window)'

	Notes
	-----
	..	The returned BioSPPy ReturnTuple object contains all frequency band parameters in parameter specific tuples
		of length 4 when using the ULF frequency band or of length 3 when NOT using the ULF frequency band.
		The structures of those tuples are shown in this example below (fft_results = ReturnTuple object returned by
		this function):

			Using ULF, VLF, LF and HF frequency bands:
				fft_results['fft_peak'] = (ulf_peak, vlf_peak, lf_peak, hf_peak)

			Using VLF, LF and HF frequency bands:
				fft_results['fft_peak'] = (vlf_peak, lf_peak, hf_peak)

	..	If 'show_param' is true, the parameters (incl. frequency band limits) will be listed next to the graph and no
		legend with frequency band limits will be added to the plot graph itself, i.e. the effect of 'show_param'
		will be used over the 'legend' effect.
	..	Only one type of input data is required.
	.. 	If both 'nni' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nni' and the 'nni' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in [s] format will be converted to [ms] format.

	"""
	# Check input values
	nn = tools.check_input(nni, rpeaks)

	# Verify or set default frequency bands
	fbands = _check_freq_bands(fbands)

	# Resampling (with 4Hz) and interpolate
	# Because RRi are unevenly spaced we must interpolate it for accurate PSD estimation.
	fs = 4
	t = np.cumsum(nn)
	t -= t[0]
	f_interpol = sp.interpolate.interp1d(t, nn, 'cubic')
	t_interpol = np.arange(t[0], t[-1], 1000./fs)
	nn_interpol = f_interpol(t_interpol)

	# Subtract mean value from each sample for surpression of DC-offsets
	if detrend:
		nn_interpol = nn_interpol - np.mean(nn_interpol)

	# Adapt 'nperseg' according to the total duration of the NNI series (5min threshold = 300000ms)
	if t.max() < 300000:
		nperseg = nfft
	else:
		nperseg = 300

	# Compute power spectral density estimation (where the magic happens)
	frequencies, powers = welch(
		x=nn_interpol,
		fs=fs,
		window=window,
		nperseg=nperseg,
		nfft=nfft,
		scaling='density'
	)

	# Compute frequency parameters
	params, freq_i = _compute_parameters('fft', frequencies, powers, fbands)

	# Plot PSD
	figure = _plot_psd('fft', frequencies, powers, freq_i, params, show, show_param, legend)
	figure = utils.ReturnTuple((figure, ), ('fft_plot', ))

	# Metadata
	args = (nfft, window, fs, 'cubic')
	names = ('fft_nfft', 'fft_window', 'fft_resampling_frequency', 'fft_interpolation', )
	meta = utils.ReturnTuple(args, names)

	# Output
	return tools.join_tuples(params, figure, meta)


def lomb_psd(
		nni=None,
		rpeaks=None,
		fbands=None,
		nfft=2**8,
		ma_size=None,
		show=True,
		show_param=True,
		legend=True
	):
	"""Computes a Power Spectral Density (PSD) estimation from the NNI series using the Lomb-Scargle Periodogram
	and computes all frequency domain parameters from this PSD according to the specified frequency bands.

	References: [Lomb1976], [Scargle1982], [Kuusela2014], [Laguna1995]
	Docs:		https://pyhrv.readthedocs.io/en/latest/_pages/api/frequency.html#lomb-scargle-periodogram-lomb-psd

	Parameters
	----------
	rpeaks : array
		R-peak locations in [ms] or [s]
	nni : array
		NN-Intervals in [ms] or [s]
	fbands : dict, optional
		Dictionary with frequency bands (2-element tuples or list)
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))´
	nfft : int, optional
		Number of points computed for the FFT result (default: 2**8)
	ma_size : int, optional
		Window size of the optional moving average filter (default: None)
	show : bool, optional
		If true, show PSD plot (default: True)
	show_param : bool, optional
		If true, list all computed PSD parameters next to the plot (default: True)
	legend : bool, optional
		If true, add a legend with frequency bands to the plot (default: True)

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All results of the Lomb-Scargle PSD estimation (see list and keys below)

	Returned Parameters & Keys
	--------------------------
	..	Peak frequencies of all frequency bands in [Hz] (key: 'lomb_peak')
	..	Absolute powers of all frequency bands in [ms^2][(key: 'lomb_abs')
	..	Relative powers of all frequency bands [%] (key: 'lomb_rel')
	..	Logarithmic powers of all frequency bands [-] (key: 'lomb_log')
	..	Normalized powers of all frequency bands [-] (key: 'lomb_norms')
	..	LF/HF ratio [-] (key: 'lomb_ratio')
	..	Total power over all frequency bands in [ms^2] (key: 'lomb_total')
	.. 	Number of PSD samples (key: 'lomb_nfft')
	.. 	Moving average filter order (key: 'lomb_ma')

	Notes
	-----
	..	The returned BioSPPy ReturnTuple object contains all frequency band parameters in parameter specific tuples
		of length 4 when using the ULF frequency band or of length 3 when NOT using the ULF frequency band.
		The structures of those tuples are shown in this example below (lomb_results = ReturnTuple object returned by
		this function):

			Using ULF, VLF, LF and HF frequency bands:
				lomb['fft_peak'] = (ulf_peak, vlf_peak, lf_peak, hf_peak)

			Using VLF, LF and HF frequency bands:
				lomb['fft_peak'] = (vlf_peak, lf_peak, hf_peak)

	..	If 'show_param' is true, the parameters (incl. frequency band limits) will be listed next to the graph and no
		legend with frequency band limits will be added to the plot graph itself, i.e. the effect of 'show_param'
		will be used over the 'legend' effect.
	..	Only one type of input data is required.
	.. 	If both 'nni' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nni' and the 'nni' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in [s] format will be converted to [ms] format.

	"""
	# Check input
	nn = tools.check_input(nni, rpeaks)

	# Verify or set default frequency bands
	fbands = _check_freq_bands(fbands)
	t = np.cumsum(nn)
	t -= t[0]

	# Compute PSD according to the Lomb-Scargle method
	# Specify frequency grid
	frequencies = np.linspace(0, 0.41, nfft)
	# Compute angular frequencies
	a_frequencies = np.asarray(2 * np.pi / frequencies)
	powers = np.asarray(lombscargle(t, nn, a_frequencies, normalize=True))
	# ms^2 to s^2
	powers = powers * 10**6

	# Fix power = inf at f=0
	powers[0] = 2

	# Apply moving average filter
	if ma_size is not None:
		powers = biosppy.signals.tools.smoother(powers, size=ma_size)['signal']

	# Compute frequency parameters
	params, freq_i = _compute_parameters('lomb', frequencies, powers, fbands)

	# Plot parameters
	figure = _plot_psd('lomb', frequencies, powers, freq_i, params, show, show_param, legend)
	figure = utils.ReturnTuple((figure, ), ('lomb_plot', ))

	# Define metadata
	meta = utils.ReturnTuple((nfft, ma_size, ), ('lomb_nfft', 'lomb_ma'))

	# Complete output
	return tools.join_tuples(params, figure, meta)


def ar_psd(nni=None,
		   rpeaks=None,
		   fbands=None,
		   nfft=2**12,
		   order=16,
		   show=True,
		   show_param=True,
		   legend=True):
	"""Computes a Power Spectral Density (PSD) estimation from the NNI series using the Autoregressive method
	and computes all frequency domain parameters from this PSD according to the specified frequency bands.

	References: [Electrophysiology1996], [Kuusela2014], [Kallas2012], [Boardman2002]
				(additionally recommended: [Miranda2012])
	Docs:		https://pyhrv.readthedocs.io/en/latest/_pages/api/frequency.html#autoregressive-method-ar-psd

	Parameters
	----------
	rpeaks : array
		R-peak locations in [ms] or [s]
	nni : array
		NN-Intervals in [ms] or [s]
	fbands : dict, optional
		Dictionary with frequency bands (2-element tuples or list)
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))´
	nfft : int, optional
		Number of points computed for the entire AR result (default: 2**12)
	order : int, optional
		Autoregressive model order (default: 16)
	show : bool, optional
		If true, show PSD plot (default: True)
	show_param : bool, optional
		If true, list all computed PSD parameters next to the plot (default: True)
	legend : bool, optional
		If true, add a legend with frequency bands to the plot (default: True)

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All results of the Autoregressive PSD estimation (see list and keys below)

	Returned Parameters & Keys
	--------------------------
	..	Peak frequencies of all frequency bands in [Hz] (key: 'ar_peak')
	..	Absolute powers of all frequency bands in [ms^2][(key: 'ar_abs')
	..	Relative powers of all frequency bands [%] (key: 'ar_rel')
	..	Logarithmic powers of all frequency bands [-] (key: 'ar_log')
	..	Normalized powers of all frequency bands [-] (key: 'ar_norms')
	..	LF/HF ratio [-] (key: 'ar_ratio')
	..	Total power over all frequency bands in [ms^2] (key: 'ar_total')
	..	Interpolation method (key: 'ar_interpolation')
	..	Resampling frequency (key: 'ar_resampling_frequency')
	.. 	AR model order (key: 'ar_order')
	.. 	Number of PSD samples (key: 'ar_nfft')

	Notes
	-----
	..	The returned BioSPPy ReturnTuple object contains all frequency band parameters in parameter specific tuples
		of length 4 when using the ULF frequency band or of length 3 when NOT using the ULF frequency band.
		The structures of those tuples are shown in this example below (lomb_results = ReturnTuple object returned by
		this function):

			Using ULF, VLF, LF and HF frequency bands:
				lomb['ar_peak'] = (ulf_peak, vlf_peak, lf_peak, hf_peak)

			Using VLF, LF and HF frequency bands:
				lomb['ar_peak'] = (vlf_peak, lf_peak, hf_peak)

	..	If 'show_param' is true, the parameters (incl. frequency band limits) will be listed next to the graph and no
		legend with frequency band limits will be added to the plot graph itself, i.e. the effect of 'show_param'
		will be used over the 'legend' effect.
	..	Only one type of input data is required.
	.. 	If both 'nni' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nni' and the 'nni' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in [s] format will be converted to [ms] format.

	"""
	# Check input
	nn = tools.check_input(nni, rpeaks)

	# Verify or set default frequency bands
	fbands = _check_freq_bands(fbands)

	# Resampling (with 4Hz) and interpolate
	# Because RRi are unevenly spaced we must interpolate it for accurate PSD estimation.
	fs = 4
	t = np.cumsum(nn)
	t -= t[0]
	f_interpol = sp.interpolate.interp1d(t, nn, 'cubic')
	t_interpol = np.arange(t[0], t[-1], 1000./fs)
	nn_interpol = f_interpol(t_interpol)

	# Compute autoregressive PSD
	ar = spectrum.pyule(data=nn_interpol, order=order, NFFT=nfft, sampling=fs, scale_by_freq=False)
	ar()

	# Get frequencies and powers
	frequencies = np.asarray(ar.frequencies())
	psd = np.asarray(ar.psd)
	powers = np.asarray(10 * np.log10(psd) * 10**3) 	# * 10**3 to compensate with ms^2 to s^2 conversion
														# in the upcoming steps

	# Compute frequency parameters
	params, freq_i = _compute_parameters('ar', frequencies, powers, fbands)

	# Define metadata
	meta = utils.ReturnTuple((nfft, order, fs, 'cubic'), ('ar_nfft', 'ar_order', 'ar_resampling_frequency',
														  'ar_interpolation'))
	params = tools.join_tuples(params, meta)

	# Plot PSD
	figure = _plot_psd('ar', frequencies, powers, freq_i, params, show, show_param, legend)
	figure = utils.ReturnTuple((figure, ), ('ar_plot', ))

	# Complete output
	return tools.join_tuples(params, figure)


def _compute_parameters(method, frequencies, power, freq_bands):
	"""Computes PSD HRV parameters from the PSD frequencies and powers.

	References: [Electrophysiology1996], [Basak2014]
	Docs:		https://pyhrv.readthedocs.io/en/latest/_pages/api/frequency.html#frequency-parameters

	Parameters
	----------
	method : str
		Method identifier ('fft', 'ar', 'lomb')
	frequencies
		Series of frequencies of the power spectral density computation.
	power : array
		Series of power-values of the power spectral density computation.
	freq_indices : array
		Indices of the frequency samples within each frequency band.
	freq_bands : dict, optional
		Dictionary with frequency bands (tuples or list).
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All results of the Lomb-Scargle PSD estimation (see list and keys below)

	Returned Parameters & Keys
	--------------------------
	(below, X = method identifier 'fft', 'ar' or 'lomb'
	..	Peak frequencies of all frequency bands in [Hz] (key: 'X_peak')
	..	Absolute powers of all frequency bands in [ms^2] (key: 'X_abs')
	..	Relative powers of all frequency bands in [%] (key: 'X_rel')
	..	Logarithmic powers of all frequency bands [-] (key: 'X_log')
	..	Normalized powers of all frequency bands [-](key: 'X_norms')
	..	LF/HF ratio [–] (key: 'X_ratio')
	..	Total power over all frequency bands in [ms^] (key: 'X_total')

	Raises
	------
	ValueError
		If parameter computation could not be made due to the lack of PSD samples ('nfft' too low)

	"""
	# Compute frequency resolution
	df = (frequencies[1] - frequencies[0])

	# Get indices of freq values within the specified freq bands
	ulf_i, vlf_i, lf_i, hf_i = _get_frequency_indices(frequencies, freq_bands)
	ulf_f, vlf_f, lf_f, hf_f = _get_frequency_arrays(frequencies, ulf_i, vlf_i, lf_i, hf_i)

	# Absolute powers
	if freq_bands['ulf'] is not None:
		ulf_power = np.sum(power[ulf_i]) * df
	vlf_power = np.sum(power[vlf_i]) * df
	lf_power = np.sum(power[lf_i]) * df
	hf_power = np.sum(power[hf_i]) * df
	abs_powers = (vlf_power, lf_power, hf_power, ) if freq_bands['ulf'] is None else (ulf_power, vlf_power, lf_power,
																					  hf_power,)
	total_power = np.sum(abs_powers)

	# Peak frequencies
	if freq_bands['ulf'] is not None:
		ulf_peak = ulf_f[np.argmax(power[ulf_i])]

	# Compute Peak values and catch exception caused if the number of PSD samples is too low
	try:
		vlf_peak = vlf_f[np.argmax(power[vlf_i])]
		lf_peak = lf_f[np.argmax(power[lf_i])]
		hf_peak = hf_f[np.argmax(power[hf_i])]
		peaks = (vlf_peak, lf_peak, hf_peak,) if freq_bands['ulf'] is None else (ulf_peak, vlf_peak, lf_peak, hf_peak,)
	except ValueError as e:
		if 'argmax of an empty sequence' in str(e):
			raise ValueError("'nfft' is too low: not enough PSD samples to compute the frequency parameters. Try to "
							 "increase 'nfft' to avoid this error.")

	# Relative, logarithmic powers & LF/HF ratio
	rels = tuple([float(x) / total_power * 100 for x in abs_powers])
	logs = tuple([float(np.log(x)) for x in abs_powers])
	ratio = float(lf_power) / hf_power

	# Normalized powers
	norms = tuple([100 * x / (lf_power + hf_power) for x in [lf_power, hf_power]])

	# Prepare parameters for plot
	args = (freq_bands, peaks, abs_powers, rels, logs, norms, ratio, total_power)
	names = (
		'%s_bands' % method, '%s_peak' % method, '%s_abs' % method,
		'%s_rel' % method, '%s_log' % method, '%s_norm' % method,
		'%s_ratio' % method, '%s_total' % method)

	# Output
	params = utils.ReturnTuple(args, names)
	freq_i = utils.ReturnTuple((ulf_i, vlf_i, lf_i, hf_i), ('ulf', 'vlf', 'lf', 'hf'))
	return params, freq_i


def _check_freq_bands(freq_bands):
	"""Checks provided frequency bands and re-orders the frequency
	boundaries if necessary.

	Example:
		If
			vlf = (0.003, 0.06) and lf = (0.04, 0.15)
		the max frequency of VLF > min frequency of lf (overlapping boundaries).

		Fixed frequency bands:
			vlf = (0.004, 0.04) and lf = (0.06, 0.15)

	Note, that the frequency bands should be reviewed in such case as some intervals will not be taken into the
	parameter computations.

	Parameters
	----------
	freq_bands : dict, optional
		Dictionary with frequency bands (tuples or list).
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))
	Returns
	-------
	freq_bands : dict
		Dictionary with valid frequency bands.
	"""
	if freq_bands is None:
		# Set default values
		ulf = None
		vlf = (0.000, 0.04)
		lf = (0.04, 0.15)
		hf = (0.15, 0.4)
		args = (ulf, vlf, lf, hf)
		names = ('ulf', 'vlf', 'lf', 'hf')
	else:
		# Check available data
		args_ = []
		names_ = []

		# ULF band
		ulf = freq_bands['ulf'] if 'ulf' in freq_bands.keys() else (0, 0)
		args_.append(ulf)
		names_.append('ulf')

		# VLF band
		vlf = freq_bands['vlf'] if 'vlf' in freq_bands.keys() else (0.003, 0.04)
		args_.append(vlf)
		names_.append('vlf')

		# LF band
		lf = freq_bands['lf'] if 'lf' in freq_bands.keys() else (0.04, 0.15)
		args_.append(lf)
		names_.append('lf')

		# HF band
		hf = freq_bands['hf'] if 'hf' in freq_bands.keys() else (0.15, 0.4)
		args_.append(hf)
		names_.append('hf')

		# Check if freq_band limits are valid
		# Rule: top frequency of a lower frequency band must not be higher than the lower frequency of a higher
		# frequency band
		invalid = False
		args_ = [list(x) for x in args_ if x is not None]
		for i, val in enumerate(args_[:-1]):
			if val != (0, 0):
				if args_[i][1] > args_[i+1][0]:
					subs = args_[i][1]
					args_[i][1] = args_[i+1][0]
					args_[i+1][0] = subs
					invalid = True
			else:
				args_[i] = None

		if invalid:
			raise ValueError("Invalid or overlapping frequency band limits.")

		args = args_
		names = names_

	return utils.ReturnTuple(args, names)


def _get_frequency_indices(freq, freq_bands):
	"""Returns list of lists where each list contains all indices of the PSD frequencies within a frequency band.
	Parameters
	----------
	freq : array
		Frequencies of the PSD.
	freq_bands : dict, optional
		Dictionary with frequency bands (tuples or list).
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))

	Returns
	-------
	indices : list of lists
		Lists with all indices of PSD frequencies of each frequency band.

	"""
	indices = []
	for key in freq_bands.keys():
		if freq_bands[key] is None:
			indices.append(None)
		else:
			indices.append(np.where((freq_bands[key][0] <= freq) & (freq <= freq_bands[key][1])))

	if len(indices) == 3:
		return None, indices[0], indices[1], indices[2]
	else:
		return indices


def _get_frequency_arrays(freq, ulf_i, vlf_i, lf_i, hf_i):
	"""Returns arrays with all frequencies within each frequency band.

	Parameters
	----------
	freq : array
		Frequencies of the PSD.
	ulf_i : array
		Indices of all freuquencies within the ULF band.
	vlf_i : array
		Indices of all freuquencies within the ULF band.
	lf_i : array
		Indices of all freuquencies within the ULF band.
	hf_i : array
		Indices of all freuquencies within the ULF band.

	Returns
	-------
	ulf_f : array
		Frequencies of the ULF band.
	vlf_f : array
		Frequencies of the VLF band.
	lf_f : array
		Frequencies of the LF band.
	hf_f : array
		Frequencies of the HF band.

	"""
	ulf_f = np.asarray(freq[ulf_i]) if ulf_i is not None else None
	vlf_f = np.asarray(freq[vlf_i])
	lf_f = np.asarray(freq[lf_i])
	hf_f = np.asarray(freq[hf_i])
	return ulf_f, vlf_f, lf_f, hf_f


def _plot_psd(method, freq, power, freq_indices, parameters, show, show_param, legend):
	"""Plots the PSD graph from a series of frequencies and power values.

	Parameters
	----------
	method : str
		Method identifier ('fft', 'ar', 'lomb')
	freq : array
		Series of frequencies of the power spectral density computation.
	power : array
		Series of power-values of the power spectral density computation.
	freq_indices : array
		Indices of the frequency samples within each frequency band.
	freq_bands : dict, optional
		Dictionary with frequency bands (tuples or list).
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))
	show : bool, optional
		If true, show PSD plot.
	show_param : bool, optional
		If true, list all computed PSD parameters next to the plot.
	legend : bool, optional
		If true, add a legend with frequency bands to the plot.

	Returns
	-------
	fig_psd : matplotlib figure
		Plot figure of the PSD graph.
	"""
	# Variables
	power = power / 10 ** 6
	fbands = parameters['%s_bands' % method]
	colors = {'ulf': 'b', 'vlf': 'yellowgreen', 'lf': 'salmon', 'hf': 'lightskyblue'}
	df = freq[1] - freq[0]

	if show_param:
		# Add second subplot with all computed parameters
		fig_psd = plt.figure(figsize=(12, 5))

		ax = fig_psd.add_subplot(121)
		ax2 = fig_psd.add_subplot(122)

		# Prepare parameter listing
		data = []
		index = 0

		for band in ['ulf', 'vlf', 'lf', 'hf']:
			if fbands[band] is not None:
				# Add frequency band specific data
				data.append(mpl.patches.Patch(facecolor=colors[band], label='%s: %.3fHz - %.3fHz' %
					(band.upper(), fbands[band][0], fbands[band][1])))
				data.append(
					mpl.patches.Patch(facecolor='white', label='Peak: %0.3f [$Hz$]' %
						parameters['%s_peak' % method][index]))
				data.append(
					mpl.patches.Patch(facecolor='white', label='Abs:  %0.3f [$ms^2$]' %
						parameters['%s_abs' % method][index]))
				data.append(
					mpl.patches.Patch(facecolor='white', label='Rel:  %0.3f [%%]' %
						parameters['%s_rel' % method][index]))
				data.append(
					mpl.patches.Patch(facecolor='white', label='Log:  %0.3f [$-$]' %
						parameters['%s_log' % method][index]))

				if band == 'lf':
					data.append(mpl.patches.Patch(facecolor='white', label='Norm: %0.3f [$-$]' %
						parameters['%s_norm' % method][0]))
					data.append(mpl.patches.Patch(facecolor='white', label=''))
				elif band == 'hf':
					data.append(mpl.patches.Patch(facecolor='white', label='Norm: %0.3f [$-$]' %
						parameters['%s_norm' % method][1]))
					data.append(mpl.patches.Patch(facecolor='white', label=''))

				# Spacings, total power and LF/HF ratio to format
				if band == 'ulf':
					data.append(mpl.patches.Patch(facecolor='white', label=''))
					data.append(mpl.patches.Patch(facecolor='white', label=''))

				if band == 'hf':
					spacing = 2 if fbands['ulf'] is not None else 8
					for i in range(spacing):
						data.append(mpl.patches.Patch(facecolor='white', label=''))

				if band == 'vlf':
					data.append(mpl.patches.Patch(facecolor='white', label=''))
					data.append(mpl.patches.Patch(facecolor='white', label=''))

				if (fbands['ulf'] is not None and  band == 'vlf') or (fbands['ulf'] is None and  band == 'lf'):
					data.append(mpl.patches.Patch(fc='white', label='Total Power: %.3f [$ms^2$]' % parameters[
						'%s_total' % method]))
					data.append(mpl.patches.Patch(fc='white', label='LF/HF: %.3f [-]' %
						parameters['%s_ratio' % method]))
				index += 1
		ax2.legend(handles=data, ncol=2, frameon=False)
		ax2.axis('off')
	else:
		fig_psd = plt.figure()
		ax = fig_psd.add_subplot(111)

	# Highlight the individual frequency bands
	for band in fbands.keys():
		if fbands[band] is not None:
			ax.fill_between(freq[freq_indices[band]], power[freq_indices[band]],
				facecolor=colors[band], label='%s: %.3fHz - %.3fHz' % (band.upper(), fbands[band][0], fbands[band][1]))

			# Add lines
			if band != 'hf':
				ax.vlines(fbands[band][1], 0, max(power) * (1 + 0.05),
					linestyle='--', alpha=0.5, linewidth=0.5)

	# Plot PSD function as line (not for Lomb as it tends to decrease the clarity of the plot)
	if method in ['fft', 'ar']:
		ax.plot(freq, power, color='grey', linewidth=0.5)

	# Add legend
	if legend and not show_param:
		ax.legend()

	# Finalize plot customization
	if method == 'fft':
		ax.set_title("PSD - Welch's Method")
	elif method == 'ar':
		ax.set_title("PSD - Autoregressive (Order %i)" % parameters['ar_order'])
	elif method == 'lomb':
		ax.set_title("PSD - Lomb-Scargle Periodogram")

	ax.grid(alpha=0.3)
	ax.set_xlabel('Frequency [$Hz$]')
	ax.set_ylabel('PSD [$s^2/Hz$]')
	ax.axis([0, fbands['hf'][1], 0, max(power) * (1 + 0.05)])

	if show:
		plt.show()

	return fig_psd


def frequency_domain(nni=None,
					 rpeaks=None,
					 signal=None,
					 sampling_rate=1000.,
					 fbands=None,
					 show=False,
					 show_param=True,
					 legend=True,
					 kwargs_welch=None,
					 kwargs_lomb=None,
					 kwargs_ar=None):
	"""Computes PSDs (Welch, Lomb, and Autoregressive), the parameters of the frequency domain, and returns
	the results in a ReturnTuple object.

	References: [Electrophysiology1996], [Kuusela2014], [Kallas2012], [Boardman2002], [Lomb1976], [Scargle1982],
				[Laguna1995], [Umberto2017], [Welch2017], [Basak2014]
	Docs:		https://pyhrv.readthedocs.io/en/latest/_pages/api/frequency.html#module-level-function-frequency-domain

	Parameters
	----------
	nni : array_like
		NN intervals in (ms) or (s)
	rpeaks : array_like
		R-peak times in (ms) or (s)
	signal : array_like
		ECG signal
	sampling_rate : int, float
		Sampling rate used for the ECG acquisition in (Hz)
	fbands : dict, optional
		Dictionary with frequency bands (tuples or list)
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.003Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))
	show : bool, optional
		If True, show PSD plots
	show_param : bool, optional
		If True, list all computed parameters next to the plots
	legend : bool, optional
		If True, adds legends to the PSD plots

	kwargs_welch : dict, optional
		Dictionary containing the kwargs for the 'welch_psd()' function:
			..	nfft : int, optional
					Number of points computed for the FFT result (default: 2**12)
			..	detrend : bool optional
					If True, detrend NNI series by subtracting the mean NNI (default: True)
			..	window : scipy window function, optional
					Window function used for PSD estimation (default: 'hamming')

	kwargs_lomb : dict, optional
		Dictionary containing the kwargs for the 'lomb_psd()' function:
			..	nfft : int, optional
					Number of points computed for the FFT result (default: 2**8)
			..	ma_size : int, optional
					Window size of the optional moving average filter (default: None)

	kwargs_ar : dict, optional
		Dictionary containing the kwargs for the 'ar_psd()' function:
			..	nfft : int, optional
					Number of points computed for the entire AR result (default: 2**12)
			..	order : int, optional
					Autoregressive model order (default: 16)
	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All results of the frequency domain PSD estimation methods (see list and keys below)

	Returned Parameters & Keys
	--------------------------
	(below, X = one of the methods 'fft', 'ar' or 'lomb')
	..	Peak frequencies of all frequency bands in [Hz] (key: 'X_peak')
	..	Absolute powers of all frequency bands in [ms^2] (key: 'X_abs')
	..	Relative powers of all frequency bands in [%] (key: 'X_rel')
	..	Logarithmic powers of all frequency bands in [-] (key: 'X_log')
	..	Normalized powers of all frequency bands in [-] (key: 'X_norms')
	..	LF/HF ratio in [-] (key: 'X_ratio')
	..	Total power over all frequency bands in [ms^2] (key: 'X_total')
	..	Number of PSD samples (key: 'X_nfft')
	..	FFT-specific: Interpolation method used for NNI interpolation (key: 'fft_interpolation')
	..	FFT-specific: Resampling frequency used for NNI interpolation (key: 'fft_resampling_frequency')
	..	FFT-specific: Window function used for PSD estimation of the Welch's method (key: 'fft_window)
	..	Lomb-specific: Moving average window size (key: 'lomb_ma')
	..  AR-specific: Model order (key: 'ar_order')
	..	AR-specific: Interpolation method used for NNI interpolation (key: 'ar_interpolation')
	..	AR-specific: Resampling frequency used for NNI interpolation (key: 'ar_resampling_frequency')

	Notes
	-----
	..	The returned BioSPPy ReturnTuple object contains all frequency band parameters in parameter specific tuples
		of length 4 when using the ULF frequency band or of length 3 when NOT using the ULF frequency band.
		The structures of those tuples are shown in this example below (fft_results = ReturnTuple object returned by
		this function):

			Using ULF, VLF, LF and HF frequency bands:
				fft_results['fft_peak'] = (ulf_peak, vlf_peak, lf_peak, hf_peak)

			Using VLF, LF and HF frequency bands:
				fft_results['fft_peak'] = (vlf_peak, lf_peak, hf_peak)

	..	If 'show_param' is true, the parameters (incl. frequency band limits) will be listed next to the graph and no
		legend with frequency band limits will be added to the plot graph itself, i.e. the effect of 'show_param'
		will be used over the 'legend' effect
	..	'fbands' and 'show' will be selected for all methods. Individually adding specific 'fbands' to 'kwargs_welch' or
		'kwargs_lomb' will have no effect
	.. 	Select the 'nfft' individually for each method using the kwargs dictionaries of the respective method(s)

	"""
	# Check input
	if signal is not None:
		rpeaks = biosppy.ecg.ecg(signal=signal, sampling_rate=sampling_rate, show=False)[2]
	elif nni is None and rpeaks is None:
		raise TypeError('No input data provided. Please specify input data.')

	# Get NNI series
	nn = tools.check_input(nni, rpeaks)

	# Check for kwargs for the 'welch_psd' function and compute the PSD
	if kwargs_welch is not None:
		if type(kwargs_welch) is not dict:
			raise TypeError("Expected <type 'dict'>, got %s: 'kwargs_welch' must be a dictionary containing "
							"parameters (keys) and values for the 'welch_psd' function." % type(kwargs_welch))

		# Supported kwargs
		available_kwargs = ['fbands', 'detrend', 'show', 'show_param', 'legend', 'window', 'nfft']

		# Unwrwap kwargs dictionary for Welch specific parameters
		detrend = kwargs_welch['detrend'] if 'detrend' in kwargs_welch.keys() else True
		window = kwargs_welch['window'] if 'window' in kwargs_welch.keys() else 'hamming'
		nfft = kwargs_welch['nfft'] if 'nfft' in kwargs_welch.keys() else 2**12

		unsupported_kwargs = []
		for args in kwargs_welch.keys():
			if args not in available_kwargs:
				unsupported_kwargs.append(args)

		# Throw warning if additional unsupported kwargs have been provided
		if unsupported_kwargs:
			warnings.warn("Unknown kwargs for 'welch_psd': %s. These kwargs have no effect." %
						  unsupported_kwargs, stacklevel=2)

		# Compute Welch's PSD with custom parameter settings
		welch_results = welch_psd(nn, fbands=fbands, detrend=detrend, show=False, show_param=show_param,
								  legend=legend, nfft=nfft, window=window)
	else:
		# Compute Welch's PSD with default values
		welch_results = welch_psd(nn, show=False, fbands=fbands, legend=legend, show_param=show_param)

	# Check for kwargs for the 'welch_psd' function and compute the PSD
	if kwargs_lomb is not None:
		if type(kwargs_lomb) is not dict:
			raise TypeError("Expected <type 'dict'>, got %s: 'kwargs_lomb' must be a dictionary containing "
							"parameters (keys) and values for the 'kwargs_lomb' function." % type(kwargs_lomb))

		# Supported kwargs
		available_kwargs = ['fbands', 'ma_size', 'show', 'show_param', 'legend', 'nfft', '']

		# Unwrwap kwargs dictionary
		nfft = kwargs_lomb['nfft'] if 'nfft' in kwargs_lomb.keys() else 2**8
		ma_size = kwargs_lomb['ma_size'] if 'ma_size' in kwargs_lomb.keys() else None

		unsupported_kwargs = []
		for args in kwargs_lomb.keys():
			if args not in available_kwargs:
				unsupported_kwargs.append(args)

		# Throw warning if additional unsupported kwargs have been provided
		if unsupported_kwargs:
			warnings.warn("Unknown kwargs for 'lomb_psd': %s. These kwargs have no effect."
						  % unsupported_kwargs, stacklevel=2)

		# Compute Welch's PSD with custom parameter settings
		lomb_results = lomb_psd(nn, fbands=fbands, ma_size=ma_size, show=False, show_param=show_param,
								legend=legend, nfft=nfft)
	else:
		# Compute Welch's PSD with default values
		lomb_results = lomb_psd(nn, show=False, fbands=fbands, legend=legend, show_param=show_param)

	# Check for kwargs for the 'ar_psd' function and compute the PSD
	if kwargs_ar is not None:
		if type(kwargs_ar) is not dict:
			raise TypeError("Expected <type 'dict'>, got %s: 'kwargs_ar' must be a dictionary containing "
							"parameters (keys) and values for the 'ar_psd' function." % type(kwargs_ar))

		# Supported kwargs
		available_kwargs = ['fbands', 'show', 'order', 'show_param', 'legend', 'window', 'nfft']

		# Unwrwap kwargs dictionary for Welch specific parameters
		nfft = kwargs_ar['nfft'] if 'nfft' in kwargs_ar.keys() else 2**12
		order = kwargs_ar['order'] if 'order' in kwargs_ar.keys() else 16

		unsupported_kwargs = []
		for args in kwargs_ar.keys():
			if args not in available_kwargs:
				unsupported_kwargs.append(args)

		# Throw warning if additional unsupported kwargs have been provided
		if unsupported_kwargs:
			warnings.warn("Unknown kwargs for 'welch_psd': %s. These kwargs have no effect." %
						  unsupported_kwargs, stacklevel=2)

		# Compute Autoregressive PSD with custom parameter settings
		ar_results = ar_psd(nn, fbands=fbands, order=order, show=False, show_param=show_param, legend=legend, nfft=nfft)
	else:
		# Compute Autoregressive PSD with default values
		ar_results = ar_psd(nn, show=False, fbands=fbands, legend=legend, show_param=show_param)

	# If plots should be shown (show all plots at once)
	if show:
		plt.show()

	# Output
	return tools.join_tuples(welch_results, lomb_results, ar_results)


if __name__ == "__main__":
	"""
	Example Script - HRV Frequency Domain Analysis
	
	"""
	# Load sample NNI series
	nni = np.load('./files/SampleNNISeries.npy')

	# Compute all frequency domain parameters and all methods
	results = frequency_domain(nni=nni)

	# Print results
	print("===========================")
	print("FREQUENCY DOMAIN PARAMETERS")
	print("===========================")

	print("Peak Frequencies:")
	print("VLF:		%f [Hz]" % results['fft_peak'][0])
	print("LF :		%f [Hz]" % results['fft_peak'][1])
	print("HLF:		%f [Hz]" % results['fft_peak'][2])

	print("Absolute Powers:")
	print("VLF:		%f [ms^2]" % results['fft_abs'][0])
	print("LF :		%f [ms^2]" % results['fft_abs'][1])
	print("HLF:		%f [ms^2]" % results['fft_abs'][2])

	print("Relative Powers:")
	print("VLF:		%f [%%]" % results['fft_rel'][0])
	print("LF :		%f [%%]" % results['fft_rel'][1])
	print("HLF:		%f [%%]" % results['fft_rel'][2])

	print("Logarithmic Powers:")
	print("VLF:		%f [-]" % results['fft_log'][0])
	print("LF :		%f [-]" % results['fft_log'][1])
	print("HLF:		%f [-]" % results['fft_log'][2])
	print("Total Power	:	%f [ms^2]" % results['fft_total'])
	print("LF/HF ratio	: 	%f [-]" % results['fft_ratio'])

	psd_plot = results['fft_plot']
	plt.show()

	# Alternatively compute the methods individually
	welch_psd(nni)
	lomb_psd(nni)
	ar_psd(nni)
