# -*- coding: utf-8 -*-
"""
Heart Rate Variability Toolbox - Nonlinear Parameters
-----------------------------------------------------

This module provides functions to compute HRV nonlinear domain parameters using R-peak locations
and/or NN interval series extracted from an ECG lead I-like signl (e.g. ECG, SpO2 or BVP sensor data).


Notes
-----
..  This module is part of the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..	This module is a contribution to the open-source biosignal processing toolbox 'BioSppy':
	https://github.com/PIA-Group/BioSPPy

Author
------
..  Pedro Gomes, Master Student, University of Applied Sciences Hamburg

Thesis Supervisors
------------------
..  Hugo Silva, PhD, Instituto de Telecomunicacoes & PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

Last Update
-----------
12-09-2018

:copyright: (c) 2018 by Pedro Gomes
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

# BioSppy imports
from biosppy import utils

# Local imports/HRV toolbox imports
import pyhrv.tools as tools


def poincare(nn=None, rpeaks=None, show=True, figsize=None, ellipse=True, vectors=True, legend=True, marker='o'):
	"""Creates Poincaré plot from a series of NN intervals or R-peak locations and derives the Poincaré related
	parameters SD1, SD2, SD1/SD2 ratio, and area of the Poincaré ellipse.

	References: [Tayel2015]

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	show : bool, optional
		If true, shows Poincaré plot (default: True).
	figsize : array_like, optional
		Matplotlib figure size (width, height) (default: (6, 6)).
	ellipse : bool, optional
		If true, shows fitted ellipse in plot (default: True).
	vectors : bool, optional
		If true, shows SD1 and SD2 vectors in plot (default: True).
	legend : bool, optional
		If True, adds legend to the Poincaré plot (default: True).
	marker : character, optional
		NNI marker in plot (default: 'o').

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	poincare_plot : matplotlib figure object
		Poincaré plot figure.
	sd1 : float
		SD1 value.
	sd2 : float, key: 'sd2'
		SD2 value.
	sd_ratio: float
		Ratio between SD2 and SD1 (SD2/SD1).
	ellipse_area : float, key
		Area of the fitted ellipse.

	Raises
	------
	TypeError
		If no input data for 'rpeaks' or 'nn' provided.
	"""
	# Check input values
	nn = tools.check_input(nn, rpeaks)

	# Prepare figure
	if figsize is None:
		figsize = (6, 6)
	fig = plt.figure(figsize=figsize)
	fig.tight_layout()
	ax = fig.add_subplot(111)
	
	# Prepare Poincaré data
	x1 = np.asarray(nn[:-1])
	x2 = np.asarray(nn[1:])

	# SD1 & SD2 Computation
	sd1 = np.std(np.subtract(x1, x2) / np.sqrt(2))
	sd2 = np.std(np.add(x1, x2) / np.sqrt(2))

	# Area of ellipse
	area = np.pi * sd1 * sd2

	# Plot
	ax.set_title(r'$Poincar\acute{e}$')
	ax.set_ylabel('$RR_{i+1}$ [ms]')
	ax.set_xlabel('$RR_i$ [ms]')
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
	return utils.ReturnTuple(args, names)


def sample_entropy(nn=None, rpeaks=None, dim=2, tolerance=None):
	"""Computes the sample entropy (sampen) of the NNI series.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	dim : int, optional
		Entropy embedding dimension (defaults: 2).
	tolerance : int, float, optional
		Tolerance distance for two vectors to be considered equal (default: std(NNI) * 0.2).

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
	nn = tools.check_input(nn, rpeaks)

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
	return utils.ReturnTuple(args, names)


def dfa(nn=None, rpeaks=None, short=None, long=None, show=True, figsize=None):
	"""Conducts detrended fluctuation analysis for short and long-term fluctuations.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	short : array_like, 2 elements
		Interval limits of the short term fluctuations (default: [4, 16]).
	long : array_like, 2 elements
		Interval limits of the long term fluctuations (default: [17, 64]).
	show : true

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	dfa_short : float
		Alpha value of the short term fluctuations.
	dfa_long : float
		Alpha value of the long term fluctuations.

	"""
	# Check input values
	nn = tools.check_input(nn, rpeaks)

	# Check intervals
	short = tools.check_interval(short, default=(4, 16))
	long = tools.check_interval(long, default=(17, 64))

	# Create arrays
	short = range(short[0], short[1] + 1)
	long = range(long[0], long[1] + 1)

	# Prepare plot
	if figsize is None:
		figsize = (6, 6)
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(111)
	ax.set_title('Detrended Fluctuation Analysis (DFA)')
	ax.set_xlabel('log n [beats]')
	ax.set_ylabel('log F(n)')

	# try:
	# Compute alpha values
	try:
		alpha1, dfa_short = nolds.dfa(nn, short, debug_data=True, overlap=False)
		alpha2, dfa_long = nolds.dfa(nn, long, debug_data=True, overlap=False)
	except ValueError:
		# If DFA could not be conducted due to insufficient number of NNIs, return an empty graph and 'nan' for alpha1/2
		warnings.warn("Not enough NNI samples for Detrended Fluctuations Analysis.")
		ax.axis([0, 1, 0, 1])
		ax.text(0.5, 0.5, '[Insufficient number of NNI samples for DFA]', horizontalalignment='center',
				verticalalignment='center')
		alpha1, alpha2 = 'nan', 'nan'
	else:
		# Plot DFA results if number of NNI were sufficent to conduct DFA
		# Plot short term DFA
		vals, flucts, poly = dfa_short[0], dfa_short[1], np.polyval(dfa_short[2], dfa_short[0])
		label = r'$ \alpha_{1}: %0.2f$' % alpha1
		ax.plot(vals, flucts, 'bo', markersize=1)
		ax.plot(vals, poly, 'b', label=label, alpha=0.7)

		# Plot long term DFA
		vals, flucts, poly = dfa_long[0], dfa_long[1], np.polyval(dfa_long[2], dfa_long[0])
		label = r'$ \alpha_{2}: %0.2f$' % alpha2
		ax.plot(vals, flucts, 'go', markersize=1)
		ax.plot(vals, poly, 'g', label=label, alpha=0.7)

		# Add legend
		ax.legend()
		ax.grid()

	# Plot axis
	if show:
		plt.show()

	# Output
	args = (fig, alpha1, alpha2,)
	return utils.ReturnTuple(args, ('dfa_plot', 'dfa_alpha1', 'dfa_alpha2', ))


def nonlinear(signal=None,
			  nn=None,
			  rpeaks=None,
			  sampling_rate=1000.,
			  show=False,
			  kwargs_poincare={},
			  kwargs_sampen={},
			  kwargs_dfa={}):
	"""Computes all time domain parameters of the HRV time domain module
		and returns them in a ReturnTuple object.

	Parameters
	----------
	ecg_signal : array_like
		ECG signal.
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	sampling_rate : int, float
		Sampling rate used for the ECG acquisition in (Hz).
	show : bool, optional
		If True, shows DFA plot.
	kwargs_poincare : dict, optional
		**kwargs for the 'poincare()' function.
	kwargs_sampen : dict, optional
		**kwargs for the 'sampen()' Sample Entropy function.
	kwargs_dfa : dict, optional
		**kwargs for the 'dfa()' function.

	**kwargs_poincare
	-----------------
	show : bool, optional
		If true, shows Poincaré plot (default: True).
	figsize : array_like, optional
		Matplotlib figure size (width, height) (default: (5, 5)).
	ellipse : bool, optional
		If true, shows fitted ellipse in plot (default: True).
	vectors : bool, optional
		If true, shows SD1 and SD2 vectors in plot (default: True).
	legend : bool, optional
		If True, adds legend to the Poincaré plot (default: True).
	marker : character, optional
		NNI marker in plot (default: 'o').

	**kwargs_sampen
	---------------
	dim : int, optional
		Entropy embedding dimension (defaults: 2).
	tolerance : int, float, optional
		Tolerance distance for two vectors to be considered equal (default: std(NNI) * 0.2).

	**kwargs_dfa
	------------
	short : array_like, 2 elements
		Interval limits of the short term fluctuations (default: [4, 12]).
	long : array_like, 2 elements
		Interval limits of the long term fluctuations (default: [13, 64]).
	figsize : array_like, 2 elements

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All time domain results.

	Returned Parameters
	-------------------
	..	SD1	(key: 'sd1')
	..	SD2 (key: 'sd2')
	..	SD2/SD1 (key: 'sd_ratio')
	..	Area of the fitted ellipse (key: 'ellipse_area')
	..	Sample Entropy (key: 'sampen')
	..	Detrended Fluctuations Analysis (short and long term fluctuations (key: 'dfa_short', 'dfa_long')

	Returned Figures
	----------------
	..	Poincaré plot (key: 'poincare_plot')

	Notes
	-----
	..	Results are stored in a biosppy.utils.ReturnTuple object and need to be accessed with the respective keys as
		done with dictionaries (see list of parameters and keys above)
	..	Provide at least one type of input data (ecg_signal, nn, or rpeaks).
	..	Input data will be prioritized in the following order: 1. ecg_signal, 2. nn, 3. rpeaks.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.
	..	Currently only calls the poincare() function.

	Raises
	------
	TypeError
		If no input data for 'nn', 'rpeaks', and 'signal' provided.

	"""
	# Check input
	if signal is not None:
		rpeaks = ecg(signal=signal, sampling_rate=sampling_rate, show=False)[2]
	elif nn is None and rpeaks is None:
		raise TypeError('No input data provided. Please specify input data.')

	# Get NNI series
	nn = tools.check_input(nn, rpeaks)

	# Call all nonlinear parameter functions (currently only Poincaré)
	results = poincare(nn, show=show, **kwargs_poincare)
	results = tools.join_tuples(results, sample_entropy(nn, **kwargs_sampen))
	results = tools.join_tuples(results, dfa(nn, show=show, **kwargs_dfa))

	# Output
	return results


if __name__ == "__main__":
	"""
	Example Script - Nonlinear Parameters
	"""
	from opensignalsreader import OpenSignalsReader
	from biosppy.signals.ecg import ecg

	# Load OpenSignals (r)evolution ECG sample file
	acq = OpenSignalsReader('./files/SampleECG.txt')
	signal = acq.signal('ECG')

	# # Filter data & get r-peak locations (ms)
	r_peaks = ecg(signal, sampling_rate=acq.sampling_rate, show=False)[2]
	nni = tools.nn_intervals(r_peaks)
	nni = np.load('./samples/series_1.npy')

	# Compute HRV parameters
	# Compute Poincaré
	res1 = poincare(nni, show=False)

	# Compute Triangular Index
	res2 = sample_entropy(nni)

	# Compute Detrended Fluctuation Analysis
	res3 = dfa(nni, show=False)

	# Join results
	results = tools.join_tuples(res1, res2, res3)

	# Results
	print("=========================")
	print("NONLINEAR ANALYSIS")
	print("=========================")
	print("Poincaré Plot")
	print("SD1:		%f (ms)" % results['sd1'])
	print("SD2:		%f (ms)" % results['sd2'])
	print("SD2/SD1: 	%f (ms)" % results['sd_ratio'])
	print("Area S:		%f (ms)" % results['ellipse_area'])
	print("Sample Entropy:	%" % results['sampen'])
	print("DFA:			%f	(ms)")

	# Alternatively use the nonlinear() function to compute all the nonlinear parameters using a single function
	nonlinear(nn=nni, show=True)
