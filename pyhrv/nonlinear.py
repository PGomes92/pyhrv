# -*- coding: utf-8 -*-
"""
pyHRV - Nonlinear Parameters
----------------------------

This module provides functions to compute HRV nonlinear domain parameters using R-peak locations
and/or NN interval series extracted from an ECG lead I-like signal (e.g. ECG, SpO2 or BVP sensor data).

Notes
-----
..  Up to v.0.3 this work has been developed within the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)"
..	You find the API reference for this module here:
	https://pyhrv.readthedocs.io/en/latest/_pages/api/nonlinear.html
.. 	See 'references.txt' for a full detailed list of references

Author
------
..  Pedro Gomes, pgomes92@gmail.com

Contributors (and former Thesis Supervisors)
--------------------------------------------
..  Hugo Silva, PhD, Instituto de Telecomunicacoes & PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

Last Update
-----------
12-11-2019

:copyright: (c) 2019 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.

"""
# Compatibility
from __future__ import division, print_function, absolute_import

# Third party libraries
import nolds
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# BioSPPy imports
import biosppy
from biosppy.signals.ecg import ecg

# Local imports/pyHRV toolbox imports
import pyhrv


def poincare(nni=None,
			 rpeaks=None,
			 show=True,
			 figsize=None,
			 ellipse=True,
			 vectors=True,
			 legend=True,
			 marker='o',
			 mode='normal'):
	"""Creates Poincaré plot from a series of NN intervals or R-peak locations and derives the Poincaré related
	parameters SD1, SD2, SD1/SD2 ratio, and area of the Poincaré ellipse.

	References: [Tayel2015][Brennan2001]

	Parameters
	----------
	nni : array
		NN intervals in [ms] or [s]
	rpeaks : array
		R-peak times in [ms] or [s]
	show : bool, optional
		If true, shows Poincaré plot (default: True)
	show : bool, optional
		If true, shows generated plot
	figsize : array, optional
		Matplotlib figure size (width, height) (default: (6, 6))
	ellipse : bool, optional
		If true, shows fitted ellipse in plot (default: True)
	vectors : bool, optional
		If true, shows SD1 and SD2 vectors in plot (default: True)
	legend : bool, optional
		If True, adds legend to the Poincaré plot (default: True)
	marker : character, optional
		NNI marker in plot (default: 'o')
		mode : string, optional
	Return mode of the function; available modes:
		'normal'	Returns frequency domain parameters and PSD plot figure in a ReturnTuple object
		'dev'		Returns frequency domain parameters, frequency and power arrays, no plot figure

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	poincare_plot : matplotlib figure object
		Poincaré plot figure
	sd1 : float
		Standard deviation (SD1) of the major axis
	sd2 : float, key: 'sd2'
		Standard deviation (SD2) of the minor axis
	sd_ratio: float
		Ratio between SD2 and SD1 (SD2/SD1)
	ellipse_area : float
		Area of the fitted ellipse

	"""
	# Check input values
	nn = pyhrv.utils.check_input(nni, rpeaks)

	# Prepare Poincaré data
	x1 = np.asarray(nn[:-1])
	x2 = np.asarray(nn[1:])

	# SD1 & SD2 Computation
	sd1 = np.std(np.subtract(x1, x2) / np.sqrt(2))
	sd2 = np.std(np.add(x1, x2) / np.sqrt(2))

	# Area of ellipse
	area = np.pi * sd1 * sd2

	# Dev:
	# Output computed SD1 & SD2 without plot
	if mode == 'dev':
		# Output
		args = (sd1, sd2, sd2 / sd1, area)
		names = ('sd1', 'sd2', 'sd_ratio', 'ellipse_area')
		return biosppy.utils.ReturnTuple(args, names)

	# Normal:
	# Same as dev but with plot
	if mode == 'normal':
		if figsize is None:
			figsize = (6, 6)
		fig = plt.figure(figsize=figsize)
		fig.tight_layout()
		ax = fig.add_subplot(111)

		ax.set_title(r'$Poincar\acute{e}$')
		ax.set_ylabel('$NNI_{i+1}$ [ms]')
		ax.set_xlabel('$NNI_i$ [ms]')
		ax.set_xlim([np.min(nn) - 50, np.max(nn) + 50])
		ax.set_ylim([np.min(nn) - 50, np.max(nn) + 50])
		ax.grid()
		ax.plot(x1, x2, 'r%s' % marker, markersize=2, alpha=0.5, zorder=3)

		# Compute mean NNI (center of the Poincaré plot)
		nn_mean = np.mean(nn)

		# Draw poincaré ellipse
		if ellipse:
			ellipse_ = mpl.patches.Ellipse((nn_mean, nn_mean), sd1 * 2, sd2 * 2, angle=-45, fc='k', zorder=1)
			ax.add_artist(ellipse_)
			ellipse_ = mpl.patches.Ellipse((nn_mean, nn_mean), sd1 * 2 - 1, sd2 * 2 - 1, angle=-45, fc='lightyellow', zorder=1)
			ax.add_artist(ellipse_)

		# Add poincaré vectors (SD1 & SD2)
		if vectors:
			arrow_head_size = 3
			na = 4
			a1 = ax.arrow(
				nn_mean, nn_mean, (-sd1 + na) * np.cos(np.deg2rad(45)), (sd1 - na) * np.sin(np.deg2rad(45)),
				head_width=arrow_head_size, head_length=arrow_head_size, fc='g', ec='g', zorder=4, linewidth=1.5)
			a2 = ax.arrow(
				nn_mean, nn_mean, (sd2 - na) * np.cos(np.deg2rad(45)), (sd2 - na) * np.sin(np.deg2rad(45)),
				head_width=arrow_head_size, head_length=arrow_head_size, fc='b', ec='b', zorder=4, linewidth=1.5)
			a3 = mpl.patches.Patch(facecolor='white', alpha=0.0)
			a4 = mpl.patches.Patch(facecolor='white', alpha=0.0)
			ax.add_line(mpl.lines.Line2D(
				(min(nn), max(nn)),
				(min(nn), max(nn)),
				c='b', ls=':', alpha=0.6))
			ax.add_line(mpl.lines.Line2D(
				(nn_mean - sd1 * np.cos(np.deg2rad(45)) * na, nn_mean + sd1 * np.cos(np.deg2rad(45)) * na),
				(nn_mean + sd1 * np.sin(np.deg2rad(45)) * na, nn_mean - sd1 * np.sin(np.deg2rad(45)) * na),
				c='g', ls=':', alpha=0.6))

			# Add legend
			if legend:
				ax.legend(
					[a1, a2, a3, a4],
					['SD1: %.3f$ms$' % sd1, 'SD2: %.3f$ms$' % sd2, 'S: %.3f$ms^2$' % area, 'SD1/SD2: %.3f' % (sd1/sd2)],
					framealpha=1)

		# Show plot
		if show:
			plt.show()

		# Output
		args = (fig, sd1, sd2, sd2/sd1, area)
		names = ('poincare_plot', 'sd1', 'sd2', 'sd_ratio', 'ellipse_area')
		return biosppy.utils.ReturnTuple(args, names)




def sample_entropy(nni=None, rpeaks=None, dim=2, tolerance=None):
	"""Computes the sample entropy (sampen) of the NNI series.

	Parameters
	----------
	nni : array
		NN intervals in [ms] or [s].
	rpeaks : array
		R-peak times in [ms] or [s].
	dim : int, optional
		Entropy embedding dimension (default: 2).
	tolerance : int, float, optional
		Tolerance distance for which the vectors to be considered equal (default: std(NNI) * 0.2).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	sample_entropy : float
		Sample entropy of the NNI series.

	Raises
	------
	TypeError
		If 'tolerance' is no numeric value.

	"""
	# Check input values
	nn = pyhrv.utils.check_input(nni, rpeaks)

	if tolerance is None:
		tolerance = np.std(nn, ddof=-1) * 0.2
	else:
		try:
			tolerance = float(tolerance)
		except:
			raise TypeError('Tolerance level cannot be converted to float.'
							'Please verify that tolerance is a numeric (int or float).')

	# Compute Sample Entropy
	sampen = float(nolds.sampen(nn, dim, tolerance))

	# Output
	args = (sampen, )
	names = ('sampen', )
	return biosppy.utils.ReturnTuple(args, names)


def dfa(nn=None, rpeaks=None, short=None, long=None, show=True, figsize=None, legend=True, mode='normal'):
	"""Conducts Detrended Fluctuation Analysis for short and long-term fluctuation of an NNI series.

	References: [Joshua2008][Kuusela2014][Fred2017]
	Docs:		https://pyhrv.readthedocs.io/en/latest/_pages/api/nonlinear.html#sample-entropy-sample-entropy

	Parameters
	----------
	nn : array
		NN intervals in [ms] or [s].
	rpeaks : array
		R-peak times in [ms] or [s].
	short : array, 2 elements
		Interval limits of the short term fluctuations (default: None: [4, 16]).
	long : array, 2 elements
		Interval limits of the long term fluctuations (default: None: [17, 64]).
	show : bool
		If True, shows DFA plot (default: True)
	legend : bool
		If True, adds legend with alpha1 and alpha2 values to the DFA plot (default: True)
	mode : string, optional
		Return mode of the function; available modes:
		'normal'	Returns frequency domain parameters and PSD plot figure in a ReturnTuple object
		'dev'		Returns frequency domain parameters, frequency and power arrays, no plot figure

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	dfa_short : float
		Alpha value of the short term fluctuations
	dfa_long : float
		Alpha value of the long term fluctuations
	dfa_plot : matplotlib plot figure
		Matplotlib plot figure of the DFA

	"""
	# Check input values
	nn = pyhrv.utils.check_input(nn, rpeaks)

	# Check intervals
	short = pyhrv.utils.check_interval(short, default=(4, 16))
	long = pyhrv.utils.check_interval(long, default=(17, 64))

	# Create arrays
	short = range(short[0], short[1] + 1)
	long = range(long[0], long[1] + 1)

	# Prepare plot
	if mode == 'normal':
		if figsize is None:
			figsize = (6, 6)
		fig = plt.figure(figsize=figsize)
		ax = fig.add_subplot(111)
		ax.set_title('Detrended Fluctuation Analysis (DFA)')
		ax.set_xlabel('log n [beats]')
		ax.set_ylabel('log F(n)')

	# Compute alpha values
	try:
		alpha1, dfa_short = nolds.dfa(nn, short, debug_data=True, overlap=False)
		alpha2, dfa_long = nolds.dfa(nn, long, debug_data=True, overlap=False)
	except ValueError:
		# If DFA could not be conducted due to insufficient number of NNIs, return an empty graph and 'nan' for alpha1/2
		warnings.warn("Not enough NNI samples for Detrended Fluctuations Analysis.")

		# Update plot
		if mode == 'normal':
			ax.axis([0, 1, 0, 1])
			ax.text(0.5, 0.5, '[Insufficient number of NNI samples for DFA]', horizontalalignment='center',
					verticalalignment='center')
		alpha1, alpha2 = 'nan', 'nan'
	else:
		# Plot DFA results if number of NNI were sufficent to conduct DFA
		# Plot short term DFA
		vals, flucts, poly = dfa_short[0], dfa_short[1], np.polyval(dfa_short[2], dfa_short[0])
		label = r'$ \alpha_{1}: %0.2f$' % alpha1

		# Update plot
		if mode == 'normal':
			ax.plot(vals, flucts, 'bo', markersize=1)
			ax.plot(vals, poly, 'b', label=label, alpha=0.7)

			# Plot long term DFA
			vals, flucts, poly = dfa_long[0], dfa_long[1], np.polyval(dfa_long[2], dfa_long[0])
			label = r'$ \alpha_{2}: %0.2f$' % alpha2
			ax.plot(vals, flucts, 'go', markersize=1)
			ax.plot(vals, poly, 'g', label=label, alpha=0.7)

			# Add legend
			if legend:
				ax.legend()
			ax.grid()

	# Normal Mode:
	# Returns results & plot figure
	if mode == 'normal':
		# Plot axis
		if show:
			plt.show()
		args = (fig, alpha1, alpha2, short, long)
		return biosppy.utils.ReturnTuple(args, ('dfa_plot', 'dfa_alpha1', 'dfa_alpha2', 'dfa_alpha1_beats', 'dfa_alpha2_beats'))

	# Dev Mode:
	# Returns results only
	if mode == 'dev':
		args = (alpha1, alpha2, short, long)
		return biosppy.utils.ReturnTuple(args, ('dfa_alpha1', 'dfa_alpha2', 'dfa_alpha1_beats', 'dfa_alpha2_beats'))


def nonlinear(nni=None,
			  rpeaks=None,
			  signal=None,
			  sampling_rate=1000.,
			  show=False,
			  kwargs_poincare=None,
			  kwargs_sampen=None,
			  kwargs_dfa=None,):
	"""Computes all time domain parameters of the HRV time domain module
		and returns them in a ReturnTuple object.

	References: [Peng1995][Willson2002]

	Parameters
	----------
	nni : array
		NN intervals in [ms] or [s].
	rpeaks : array
		R-peak times in [ms] or [s].
	signal : array
		ECG signal.
	sampling_rate : int, float
		Sampling rate used for the ECG acquisition in (Hz).
	show : bool, optional
		If True, shows DFA plot.

	kwargs_poincare : dict, optional
		Dictionary containing the kwargs for the nonlinear 'poincare()' function:
			..	ellipse : bool, optional
					If true, shows fitted ellipse in plot (default: True).
			..	vectors : bool, optional
					If true, shows SD1 and SD2 vectors in plot (default: True).
			..	legend : bool, optional
					If True, adds legend to the Poincaré plot (default: True).
			..	marker : character, optional
					NNI marker in plot (default: 'o').

	kwargs_dfa : dict, optional
		Dictionary containing the kwargs for the nonlinear 'dfa()' function:
			..	short : array, 2 elements
					Interval limits of the short term fluctuations (default: None: [4, 16]).
			..	long : array, 2 elements
					Interval limits of the long term fluctuations (default: None: [17, 64]).
			..	legend : bool
					If True, adds legend with alpha1 and alpha2 values to the DFA plot (default: True)

	kwargs_sampen : dict, optional
		Dictionary containing the kwargs for the nonlinear 'sample_entropy()' function:
			..	dim : int, optional
					Entropy embedding dimension (default: 2).
			..	tolerance : int, float, optional
					Tolerance distance for which the vectors to be considered equal (default: std(NNI) * 0.2).

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All time domain results.

	Returned Parameters
	-------------------
	..	SD1	in [ms] (key: 'sd1')
	..	SD2 in [ms] (key: 'sd2')
	..	SD2/SD1 [-] (key: 'sd_ratio')
	..	Area of the fitted ellipse in [ms^2] (key: 'ellipse_area')
	..	Sample Entropy [-] (key: 'sampen')
	..	Detrended Fluctuations Analysis [-] (short and long term fluctuations (key: 'dfa_short', 'dfa_long')

	Returned Figures
	----------------
	..	Poincaré plot (key: 'poincare_plot')

	Notes
	-----
	..	Results are stored in a biosppy.utils.ReturnTuple object and need to be accessed with the respective keys as
		done with dictionaries (see list of parameters and keys above)
	..	Provide at least one type of input data (signal, nn, or rpeaks).
	..	Input data will be prioritized in the following order: 1. signal, 2. nn, 3. rpeaks.
	..	NN and R-peak series provided in [s] format will be converted to [ms] format.
	..	Currently only calls the poincare() function.

	Raises
	------
	TypeError
		If no input data for 'nn', 'rpeaks', and 'signal' provided.

	"""
	# Check input
	if signal is not None:
		t, signal, rpeaks = biosppy.ecg.ecg(signal=signal, sampling_rate=sampling_rate, show=False)[:3]
		rpeaks = t[rpeaks]
	elif nni is None and rpeaks is None:
		raise TypeError('No input data provided. Please specify input data.')

	# Get NNI series
	nn = pyhrv.utils.check_input(nni, rpeaks)

	# Unwrap kwargs_poincare dictionary & compute Poincaré
	if kwargs_poincare is not None:
		if type(kwargs_poincare) is not dict:
			raise TypeError("Expected <type 'dict'>, got %s: 'kwargs_poincare' must be a dictionary containing "
							"parameters (keys) and values for the 'poincare()' function." % type(kwargs_poincare))

		# Supported kwargs
		available_kwargs = ['ellipse', 'vectors', 'legend', 'marker']

		# Unwrwap kwargs dictionaries
		ellipse = kwargs_poincare['ellipse'] if 'ellipse' in kwargs_poincare.keys() else True
		vectors = kwargs_poincare['vectors'] if 'vectors' in kwargs_poincare.keys() else True
		legend = kwargs_poincare['legend'] if 'legend' in kwargs_poincare.keys() else True
		marker = kwargs_poincare['marker'] if 'marker' in kwargs_poincare.keys() else 'o'

		unsupported_kwargs = []
		for args in kwargs_poincare.keys():
			if args not in available_kwargs:
				unsupported_kwargs.append(args)

		# Throw warning if additional unsupported kwargs have been provided
		if unsupported_kwargs:
			warnings.warn("Unknown kwargs for 'poincare()': %s. These kwargs have no effect."
						  % unsupported_kwargs, stacklevel=2)

		# Compute Poincaré plot with custom configuration
		p_results = poincare(nn, show=False, ellipse=ellipse, vectors=vectors, legend=legend, marker=marker)
	else:
		# Compute Poincaré plot with default values
		p_results = poincare(nn, show=False)

	# Unwrap kwargs_sampen dictionary & compute Sample Entropy
	if kwargs_sampen is not None:
		if type(kwargs_sampen) is not dict:
			raise TypeError("Expected <type 'dict'>, got %s: 'kwargs_sampen' must be a dictionary containing "
							"parameters (keys) and values for the 'sample_entropy()' function." % type(kwargs_sampen))

		# Supported kwargs
		available_kwargs = ['dim', 'tolerance']

		# Unwrwap kwargs dictionaries
		dim = kwargs_sampen['dim'] if 'dim' in kwargs_sampen.keys() else 2
		tolerance = kwargs_sampen['tolerance'] if 'tolerance' in kwargs_sampen.keys() else None

		unsupported_kwargs = []
		for args in kwargs_sampen.keys():
			if args not in available_kwargs:
				unsupported_kwargs.append(args)

		# Throw warning if additional unsupported kwargs have been provided
		if unsupported_kwargs:
			warnings.warn("Unknown kwargs for 'sample_entropy()': %s. These kwargs have no effect."
						  % unsupported_kwargs, stacklevel=2)

		# Compute Poincaré plot with custom configuration
		s_results = sample_entropy(nn, dim=dim, tolerance=tolerance)
	else:
		# Compute Poincaré plot with default values
		s_results = sample_entropy(nn)

	# Unwrap kwargs_dfa dictionary & compute Poincaré
	if kwargs_dfa is not None:
		if type(kwargs_dfa) is not dict:
			raise TypeError("Expected <type 'dict'>, got %s: 'kwargs_dfa' must be a dictionary containing "
							"parameters (keys) and values for the 'dfa()' function." % type(kwargs_dfa))

		# Supported kwargs
		available_kwargs = ['short', 'legend', 'long']

		# Unwrwap kwargs dictionaries
		short = kwargs_dfa['short'] if 'short' in kwargs_dfa.keys() else None
		long = kwargs_dfa['long'] if 'long' in kwargs_dfa.keys() else None
		legend = kwargs_dfa['legend'] if 'legend' in kwargs_dfa.keys() else True

		unsupported_kwargs = []
		for args in kwargs_dfa.keys():
			if args not in available_kwargs:
				unsupported_kwargs.append(args)

		# Throw warning if additional unsupported kwargs have been provided
		if unsupported_kwargs:
			warnings.warn("Unknown kwargs for 'dfa()': %s. These kwargs have no effect."
						  % unsupported_kwargs, stacklevel=2)

		# Compute DFA plot with custom configuration
		d_results = dfa(nn, show=False, short=short, long=long, legend=legend)
	else:
		# Compute DFA plot with default values
		d_results = dfa(nn, show=False)

	# Join Results
	results = pyhrv.utils.join_tuples(p_results, s_results, d_results)

	# Plot
	if show:
		plt.show()

	# Output
	return results


if __name__ == "__main__":
	"""
	Example Script - Nonlinear Parameters
	"""
	import numpy as np

	# Load sample NNI series
	nni = pyhrv.utils.load_sample_nni()

	# Compute Poincaré
	res1 = poincare(nni, show=False)

	# Compute Triangular Index
	res2 = sample_entropy(nni)

	# Compute Detrended Fluctuation Analysis
	res3 = dfa(nni, show=False)

	# Join results
	results = pyhrv.utils.join_tuples(res1, res2, res3)

	# Results
	print("=========================")
	print("NONLINEAR ANALYSIS")
	print("=========================")
	print("Poincaré Plot")
	print("SD1:				%f [ms]" % results['sd1'])
	print("SD2:				%f [ms]" % results['sd2'])
	print("SD2/SD1: 		%f [ms]" % results['sd_ratio'])
	print("Area S:			%f [ms]" % results['ellipse_area'])
	print("Sample Entropy:	%f" % results['sampen'])
	print("DFA (alpha1):	%f	[ms]" % results['dfa_alpha1'])
	print("DFA (alpha2):	%f	[ms]" % results['dfa_alpha2'])

	# Alternatively use the nonlinear() function to compute all the nonlinear parameters at once
	nonlinear(nni=nni)
