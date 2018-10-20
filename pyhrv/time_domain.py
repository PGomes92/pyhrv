# -*- coding: utf-8 -*-
"""
Heart Rate Variability Toolbox - Time Domain Module
---------------------------------------------------

This module provides functions to compute HRV time domain  parameters using R-peak locations
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
..  Hugo Silva, PhD, Instituto de Telecomunicoes & PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

Last Update
-----------
12-09-2018

:copyright: (c) 2018 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.
"""
# Compatibility
from __future__ import division, print_function

# Third party libraries
import warnings
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# BioSppy imports
from biosppy.signals.ecg import ecg
from biosppy import utils

# Local imports/HRV toolbox imports
import pyhrv.tools as tools

# Suppress SciPy warnings
# warnings.filterwarnings('error', category=FutureWarning)


def nn_parameters(nn=None, rpeaks=None):
	"""Computes basic statistical parameters from a series of NN intervals (# of intervals, mean, min, max).

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	nn_counter : int
		Number of NN intervals.
	nn_mean : float
		Mean NN interval (ms).
	nn_min : float
		Minimum NN interval (ms).
	nn_max : float
		Maximum NN interval (ms).

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# output
	args = (int(nn.size), nn.mean(), nn.min(), nn.max())
	names = ('nn_counter', 'nn_mean', 'nn_min', 'nn_max')
	return utils.ReturnTuple(args, names)


def nn_differences_parameters(nn=None, rpeaks=None):
	"""Computes basic statistical parameters from a series of successive NN interval differences (mean, min, max,
	standard deviation).

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	nn_diff_mean: float
		Mean NN interval difference (ms).
	nn_diff_min : float
		Minimum NN interval difference (ms).
	nn_diff_max : float
		Maximum NN interval difference (ms).
	nn_diff_std : float
		Standard deviation of the NN interval differences series (ms).

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Get NN interval differences
	nnd = tools.nn_diff(nn)

	# output
	args = (float(nnd.mean()), int(nnd.min()), int(nnd.max()), float(nnd.std(ddof=1)))
	names = ('nn_diff_mean', 'nn_diff_min', 'nn_diff_max', 'nn_diff_std')
	return utils.ReturnTuple(args, names)


def hr_parameters(nn=None, rpeaks=None):
	"""Computes basic statistical parameters from a series of Heart Rate (HR) data (mean, min, max, standard deviation).

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	bpm_mean : float
		Mean heart rate (bpm).
	bpm_min : float
		Minimum heart rate value (bpm).
	bpm_max : float
		Maximum heart rate value (bpm).
	bpm_std : float
		Standard deviation of the HR series (bpm).

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Get heart rate series
	hr = tools.heart_rate(nn)

	# Output
	args = (hr.mean(), hr.min(), hr.max(), hr.std(ddof=1))
	names = ('hr_mean', 'hr_min', 'hr_max', 'hr_std')
	return utils.ReturnTuple(args, names)


def sdnn(nn=None, rpeaks=None):
	"""Computation of the standard deviation of a NN interval series.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	sdnn : float
		Standard deviation of NN intervals (ms).

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Computation of SDNN & Output
	args = [tools.std(nn)]
	names = ['sdnn']
	return utils.ReturnTuple(args, names)


def sdnn_index(nn=None, rpeaks=None, full=False, overlap=False, duration=300):
	"""Computes the mean of the SDNN values of each segment (default: 300s segments).

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	full : bool, optional
		If True, returns last segment, even if the cumulative sum of NNI does not reach the 300s (default: False).
	overlap : bool, optional
		If True, allow to return NNI that go from the interval of one segment to the successive segment (default: False).
	duration : int, optional
		Maximum duration duration per segment in (s) (default: 300s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	sdnn_index : float
		Mean of the standard deviations of all NN intervals within 5 minutes intervals (ms)

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.
	..	In some cases, the NN interval may start in a segment (or time interval) N and end only in the successive
		segment N+1. In this case, use the 'overlap' parameter to select if the first element of the segment should be
		dropped or not:
		..	If True: overlap allowed, returns all NNI but the cumulative sum of the NNI in a segment can be greater
			than the specified duration.
		..	If False: no overlap allowed, first NNI will be dropped and the cumulative sum of the NNI in a segment
			will always be < specified duration.
	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Signal segmentation into 5 min segments
	segments, seg = tools.segmentation(nn,  full=full, overlap=overlap, duration=duration)

	if seg:
		sdnn_values = [sdnn(x) for x in segments]
		sdnn_index = np.mean(sdnn_values)
	else:
		sdnn_index = float('nan')
		if tools.WARN:
			warnings.warn("Signal duration too short for SDNN index computation.")

	# Output
	args = [sdnn_index]
	names = ['sdnn_index']
	return utils.ReturnTuple(args, names)


def sdann(nn=None, rpeaks=None, full=False, overlap=False, duration=300):
	"""Computes the standard deviation of the mean NNI value of each segment (default: 300s segments).

	Parameters
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	full : bool, optional
		If True, returns last segment, even if the cumulative sum of NNI does not reach the 300s (default: False).
	overlap : bool, optional
		If True, allow to return NNI that go from the interval of one segment to the successive segment (default: False).
	duration : int, optional
		Maximum duration duration per segment in (s) (default: 300s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	sdnn_index : float
		Standard deviations of the means of all NN intervals within 5 minutes intervals in (ms).

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.
	..	In some cases, the NN interval may start in a segment (or time interval) N and end only in the successive
		segment N+1. In this case, use the 'overlap' parameter to select if the first element of the segment should be
		dropped or not:
		..	If True: overlap allowed, returns all NNI but the cumulative sum of the NNI in a segment can be greater
			than the specified duration.
		..	If False: no overlap allowed, first NNI will be dropped and the cumulative sum of the NNI in a segment
			will always be < specified duration.
	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Signal segmentation into 5 min segments
	segments, seg = tools.segmentation(nn, full=full, overlap=overlap, duration=duration)

	if seg:
		mean_values = [np.mean(x) for x in segments]
		sdann_ = tools.std(mean_values)
	else:
		sdann_ = float('nan')
		if tools.WARN:
			warnings.warn("Signal duration too short for SDANN computation.")

	# Output
	args = [sdann_]
	names = ['sdann']
	return utils.ReturnTuple(args, names)


def rmssd(nn=None, rpeaks=None):
	"""Computes root mean of squared differences of successive NN Intervals.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	rmssd : float
		RMSSD value in (ms).

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Compute RMSSD
	nnd = tools.nn_diff(nn)
	rmssd = np.sum(x**2 for x in nnd)
	rmssd = np.sqrt(1. / nnd.size * rmssd)

	# Output
	args = [rmssd]
	names = ['rmssd']
	return utils.ReturnTuple(args, names)


def sdsd(nn=None, rpeaks=None):
	"""Computation of the standard deviation of differences of successive NN intervals.
	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	sdsd : float
		Standard deviation of successive differences of NN intervals (ms)

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Compute NN differences
	nnd = tools.nn_diff(nn)

	# Computation of SDNN
	sdsd_ = tools.std(nnd)

	# Output
	args = [sdsd_]
	names = ['sdsd']
	return utils.ReturnTuple(args, names)


def nnXX(nn=None, rpeaks=None, threshold=None):
	"""Find number of NN interval differences greater than a specified threshold and ratio between number of intervals
	> threshold and total number of NN interval differences.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	threshold : int
		Threshold for nnXX values in (ms).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	nnXX: int
		Number of NN interval differences greater than the specified threshold (-).
	pnnXX : float
		Ratio between nnXX and total number of NN interval differences (-).

	Notes
	-----
	..	Results are stored in a biosppy.utils.ReturnTuple object and need to be accessed with the respective keys as
		done with dictionaries (see list of parameters and keys above)
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input values
	nn = tools.check_input(nn, rpeaks)

	# Count NN20
	nnd = tools.nn_diff(nn)
	nnxx = sum(i > threshold for i in nnd)
	pnnxx = nnxx / len(nnd) * 100

	# Output
	args = (nnxx, pnnxx)
	names = ('nn%i' % threshold, 'pnn%i' % threshold)
	return utils.ReturnTuple(args, names)


def nn50(nn=None, rpeaks=None):
	"""Find number of NN interval differences which are greater 50ms (NN50) and ratio between NN50 and total amount of
	NN intervals.

	References: [Hutchinson2003] [Mietus2002]

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	nn50 : int
		Number of NN interval differences greater 50ms.
	pnn50 : float
		Ratio between NN50 and total number of NN intervals.

	Raises
	------
	TypeError
		If no input data for 'rpeaks' or 'nn' provided.

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input values
	if nn is not None:
		return nnXX(nn, threshold=50)
	elif r_peaks is not None:
		return nnXX(rpeaks=rpeaks, threshold=50)
	else:
		raise TypeError("No data for r_peaks or nn_intervals provided. Please specify input data.")


def nn20(nn=None, rpeaks=None):
	"""Find number of NN interval differences which are greater 20ms (NN20) and ratio between NN20 and total amount of
	NN intervals.

	References: [Hutchinson2003] [Mietus2002]

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).

	Returns
	-------
	nn20 : int
		Number of NN interval differences greater 20ms.
	pnn20 : float
		Ratio between NN20 and total number of NN intervals.

	Raises
	------
	TypeError
		If no input data for 'rpeaks' or 'nn' provided.

	Notes
	-----
	..	Provide at least one type of input data.
	.. 	If both 'nn' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nn' data will be computed
		from the 'rpeaks'.
	..	NN and R-peak series provided in (s) format will be converted to (ms) format.

	"""
	# Check input values
	if nn is not None:
		return nnXX(nn, threshold=20)
	elif r_peaks is not None:
		return nnXX(rpeaks=rpeaks, threshold=20)
	else:
		raise TypeError("No data for r_peaks or nn_intervals provided. Please specify input data.")


def tinn(nn=None, rpeaks=None, binsize=7.815, plot=True, show=True, figsize=None, legend=True):
	"""Computes TINN based on the NN intervals histogram.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	binsize : int, float
		Bin size of the histogram bins (default: 7.8125ms).
	plot : bool
		If True, creates histogram plot using matplotlib, else uses numpy (data only, no plot).
	show : bool, optional
		If true, shows histogram (default: True).
	figsize : array_like, optional
		Matplotlib figure size (width, height) (default: (6, 6)).
	legend : bool, optional
		If True, adds legend to the histogram (default: True).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	tinn_histogram : matplotlib figure object
		Histogram figure (only if input parameter 'plot' is True).
	tinn_n : float
		N value of the TINN computation.
	tinn_m : float
		M value of the TINN computation.
	tinn : float
		TINN value.

	Raises
	------
	TypeError (via 'check_input()')
		If no input data for 'rpeaks' or 'nn' provided.

	Notes
	-----
	..	Default bin size set to recommended bin size of 1/128 (with 128Hz being the minimum recommended sampling
		frequency) as recommended by the HRV guidelines REFERENCE.
	..	'show' has only effect if 'plot' is also True.
	.. 	'legend' has only effect if 'plot' is also True.
	..	'figsize' has only effect if 'plot' is also True.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Get Histogram data (with or without histogram plot figure)
	if plot:
		fig, ax, D, bins = _get_histogram(nn, figsize=figsize, binsize=binsize, legend=legend, plot=plot)
	else:
		D, bins = _get_histogram(nn, figsize=figsize, binsize=binsize, legend=legend, plot=plot)

	# Use only all bins except the last one to avoid indexing error with 'D' (otherwise bins.size = D.size + 1)
	# bins = np.asarray(bins[:-1])

	# Get bins of the triangular's N side (left side of the bin with the highest value of the distribution)
	n_bins = [bin for bin in bins if bin < bins[np.argmax(D)]]

	# Get bins of the triangular's M side (right side of the bin with the highest value of the distribution)
	m_bins = [bin for bin in bins if bin > bins[np.argmax(D)]]

	# Set a maximum error
	min_error = 2 ** 14
	N = 0
	M = 0

	# Compute triangle and error for each possible N value within the bins
	for n in n_bins:
		# Compute triangle and error for each possible M value within the bins
		for m in m_bins:
			# Get bin indices and time array that are valid for q(t) (i.e. q(t)!=0 for N < t < M)
			qi = np.zeros(bins.size)
			for i, bin in enumerate(bins):
				qi[i] = (True if n <= bin <= m else False)
			t = bins[[i for i, q in enumerate(qi) if q]]

			# Compute linear function that describes the N side of the triangle (q(N) = 0 to q(X) = max(D(X)))
			qn = interp1d([t[0], bins[np.argmax(D)]], [0, np.max(D)], 'linear', bounds_error=False)
			qn = qn(bins)

			# Compute linear function that describes the M side of the triangle (q(X) = max(D(X)) to q(M) = 0)
			qm = interp1d([bins[np.argmax(D)], t[-1]], [np.max(D), 0], 'linear', bounds_error=False)
			qm = qm(bins)

			# Join the linear functions of both sides to single array
			q = np.zeros(len(bins))
			for i, val in enumerate(bins):
				if str(qn[i]) != 'nan':
					q[i] = qn[i]
				elif str(qm[i]) != 'nan':
					q[i] = qm[i]
				else:
					q[i] = 0

			# Compute squared error
			error = np.sum([(D[i] - q[i]) ** 2 for i, _ in enumerate(bins)])

			# Save N and M value if error is < smaller than before
			if error < min_error:
				N, M, min_error = n, m, error
				qf = q

	# Compute TINN
	tinn = M - N

	# If plot figure required, add interpolated triangle and other specified plot characteristics
	if plot:
		# Add triangle to the plot
		ax.plot([N, bins[np.argmax(D)]], [0, D.max()], 'r--', linewidth=0.8)
		ax.plot([bins[np.argmax(D)], M], [D.max(), 0], 'r--', linewidth=0.8)

		# Add legend
		if legend:
			h = mpl.patches.Patch(facecolor='skyblue')
			tri = mpl.lines.Line2D([0, 0], [0, 0], linestyle='--', linewidth=0.8, color='r')
			x = mpl.patches.Patch(facecolor='g', alpha=0.0)
			dx = mpl.patches.Patch(facecolor='g', alpha=0.0)
			n = mpl.patches.Patch(facecolor='white', alpha=0.0)
			m = mpl.patches.Patch(facecolor='white', alpha=0.0)
			tinn_ = mpl.patches.Patch(facecolor='white', alpha=0.0)
			ax.legend(
				[h, tri, x, dx, n, m, tinn_],
				['Histogram D(NNI)', 'Triangular Interpol.', 'D(X): %i' % D.max(), 'X: %.3f$ms$' % bins[np.argmax(D)],
				 'N: %.3f$ms$' % N, 'M: %.3fms' % M, 'TINN: %.3fms' % tinn],
				loc=0
			)

		# Show plot
		if show:
			plt.show()

		# Output
		args = (fig, N, M, tinn,)
		names = ('tinn_histogram', 'tinn_n', 'tinn_m', 'tinn',)
	else:
		# Output
		args = (N, M, tinn,)
		names = ('tinn_n', 'tinn_m', 'tinn',)

	return utils.ReturnTuple(args, names)


def triangular_index(nn=None, rpeaks=None, binsize=7.815, plot=True, show=True, figsize=None, legend=True):
	"""Computes triangular index based on the NN intervals histogram.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	binsize : int, float
		Bin size of the histogram bins (default: 7.8125ms).
	plot : bool
		If True, creates histogram plot using matplotlib, else uses numpy (data only, no plot).
	show : bool, optional
		If true, shows histogram (default: True).
	figsize : array_like, optional
		Matplotlib figure size (width, height) (default: (6, 6)).
	legend : bool, optional
		If True, adds legend to the histogram (default: True).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	tri_histogram : matplotlib figure object
		Histogram figure (only if input parameter 'plot' is True).
	tri_index : float
		Triangular index.

	Raises
	------
	TypeError
		If no input data for 'rpeaks' or 'nn' provided.

	Notes
	-----
	..	Default bin size set to recommended bin size of 1/128 (with 128Hz being the minimum recommended sampling
		frequency) as recommended by the HRV guidelines REFERENCE.
	..	'show' has only effect if 'plot' is also True.
	.. 	'legend' has only effect if 'plot' is also True.
	..	'figsize' has only effect if 'plot' is also True.

	"""
	# Check input
	nn = tools.check_input(nn, rpeaks)

	# If histogram should be plotted
	if plot:
		# Get histogram values
		fig, ax, D, bins = _get_histogram(nn, figsize=figsize, binsize=binsize, legend=legend, plot=plot)

		# Compute Triangular index: number of nn intervals / maximum value of the distribution
		tri_index = nn.size / D.max()

		# Add legend
		if legend:
			h = mpl.patches.Patch(facecolor='skyblue')
			x = mpl.patches.Patch(facecolor='g', alpha=0.0)
			dx = mpl.patches.Patch(facecolor='g', alpha=0.0)
			tri = mpl.patches.Patch(facecolor='white', alpha=0.0)
			ax.legend(
				[h, x, dx, tri],
				['Histogram D(NNI)', 'D(X): %i' % D.max(), 'X: %.3f' % bins[np.argmax(D)],
				 'TriIndex: %.3f' % tri_index],
				loc=0
			)

		# Show plot
		if show:
			fig.show()

		# Output
		args = (fig, tri_index,)
		names = ('tri_histogram', 'tri_index',)

	# If histogram should not be plotted
	else:
		D, bins = _get_histogram(nn, figsize=figsize, binsize=binsize, legend=legend, plot=plot)

		# Compute Triangular index: number of nn intervals / maximum value of the distribution
		tri_index = nn.size / D.max()

		# Output
		args = (tri_index, )
		names = ('tri_index', )

	return utils.ReturnTuple(args, names)


def _get_histogram(nn=None, plot=True, figsize=None, binsize=None, legend=True):
	"""Prepares NNI histogram data for all geometrical functions.

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	plot : bool
		If True, creates histogram plot using matplotlib, else uses numpy (data only, no plot).
	figsize : array_like, optional
		Matplotlib figure size (width, height) (default: (6, 6)).
	binsize : int, float
		Bin size of the histogram bins.
	legend : bool
		If True, highlights D(X) marker to the plot to be added to the legends (default=True).

	Returns
	-------
	fig : matplotlib figure object
		Figure of the histogram plot (only if input parameter 'plot' is True).
	vals : array_like
		Histogram distribution values.
	bins : array_like
		Histogram bins.

	Raises
	------
	TypeError
		If no input data provided for 'nn'.
	TypeError
		If no input data provided for 'binsize'.

	Notes
	-----
	..	'figsize' has only effect if 'plot' is also True.
	.. 	'legend' has only effect if 'plot' is also True.

	"""
	# Check input data & confirm numpy
	if nn is None:
		raise TypeError("No input data provided for 'nn'.")
	else:
		nn = np.asarray(nn)

	if binsize is None:
		raise TypeError("No input data provided for 'binsize'")

	# Create bins array
	bins = np.arange(0, np.max(nn) + binsize, binsize)

	# Get histogram plot and data
	if plot:
		# Check figsize
		if figsize is None:
			figsize = (6, 6)

		# Prepare plot figure
		fig = plt.figure(figsize=figsize)
		ax = fig.add_subplot(111)
		vals, bins, patches = ax.hist(nn, bins, density=False, align='left', facecolor='skyblue', edgecolor='black')
		bins = bins[:-1]

		# Highlight bin of the histograms maximum value with a different color and prepare legend
		if legend:
			ax.vlines(bins[np.argmax(vals)], 0, (vals.max() * 1.1),
					  linestyles='--', color='g', linewidth=0.6)
			pos = (bins[np.argmax(vals)], vals.max() * 1.11)
			ax.annotate('D(X)', xy=pos, xytext=pos, ha='center', color='g')

		# Configure figure and plot
		ax.axis([nn.min() - (3 * binsize), nn.max() + (3 * binsize), 0, vals.max() * 1.15])
		ax.set_xlabel('NNI Bins [ms]')
		ax.set_ylabel('D(NNI) [-]')
		ax.set_title('NNI Histogram')
		return fig, ax, vals, bins

	else:
		vals, bins = np.histogram(nn, bins, density=False)
		return vals, bins[:-1]


def geometrical_parameters(nn=None, rpeaks=None, binsize=7.815, plot=True, show=True, figsize=None, legend=True):
	"""Creates NNI histogram with specified binsize (default: 7.815ms) and computes geometrical parameters (triangular
	index, TINN, N, and M).

	Parameters
	----------
	nn : array_like
		NN intervals in (ms) or (s).
	rpeaks : array_like
		R-peak times in (ms) or (s).
	binsize : int, float
		Bin size of the histogram bins (default: 7.8125ms).
	plot : bool
		If True, creates histogram plot using matplotlib, else uses numpy (data only, no plot).
	show : bool, optional
		If true, shows histogram (default: True).
	figsize : array_like, optional
		Matplotlib figure size (width, height) (default: (6, 6)).
	legend : bool, optional
		If True, adds legend to the histogram (default: True).

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	nn_histogram : matplotlib figure object
		Histogram figure (only if input parameter 'plot' is True).
	tri_index : float
		Triangular index.
	tinn_n : float
		N value of the TINN computation.
	tinn_m : float
		M value of the TINN computation.
	tinn : float
		TINN value.

	Raises
	------
	TypeError (via 'check_input()')
		If no input data for 'rpeaks' or 'nn' provided.

	Notes
	-----
	..	Default bin size set to recommended bin size of 1/128 (with 128Hz being the minimum recommended sampling
		frequency) as recommended by the HRV guidelines REFERENCE.
	..	'show' has only effect if 'plot' is also True.
	.. 	'legend' has only effect if 'plot' is also True.
	..	'figsize' has only effect if 'plot' is also True.

	"""

	# Check input
	nn = tools.check_input(nn, rpeaks)

	# Get Histogram data & plot (optional)
	if plot:
		fig, ax, D, bins = _get_histogram(nn, figsize=figsize, binsize=binsize, legend=legend, plot=plot)
	else:
		fig = None

	# Get TINN values without plot figure
	tinn_vals = tinn(nn=nn, rpeaks=rpeaks, binsize=binsize, show=False, legend=False, figsize=figsize, plot=False)

	# Get triangular index without plot figure
	trindex = triangular_index(nn=nn, rpeaks=rpeaks, binsize=binsize, show=False, legend=False, plot=False)['tri_index']

	# Histogram plot & settings
	if plot:
		# Plot triangular interpolation
		N, M = tinn_vals['tinn_n'], tinn_vals['tinn_m']
		ax.plot([N, bins[np.argmax(D)]], [0, D.max()], 'r--', linewidth=0.8)
		ax.plot([bins[np.argmax(D)], M], [D.max(), 0], 'r--', linewidth=0.8)

		# Add Legend
		if legend:
			l1 = mpl.patches.Patch(facecolor='skyblue', label='Histogram D(NNI)')
			l2 = mpl.lines.Line2D([0, 0], [0, 0], linestyle='--', linewidth=0.8, color='r', label='Tri. Interpol.')
			l3 = mpl.patches.Patch(facecolor='g', alpha=0.0, label='D(X): %i' % D.max())
			l4 = mpl.patches.Patch(facecolor='g', alpha=0.0, label='X: %.3f$ms$' % bins[np.argmax(D)])
			l5 = mpl.patches.Patch(facecolor='white', alpha=0.0, label='N: %.3f$ms$' % tinn_vals['tinn_n'])
			l6 = mpl.patches.Patch(facecolor='white', alpha=0.0, label='M: %.3fms' % tinn_vals['tinn_m'])
			l7 = mpl.patches.Patch(facecolor='white', alpha=0.0, label='TINN: %.3fms' % tinn_vals['tinn'])
			l8 = mpl.patches.Patch(facecolor='white', alpha=0.0, label='Tri. Index: %.3f' % trindex)
			ax.legend(handles=[l1, l2, l3, l4, l5, l6, l7, l8], loc=0, ncol=1)

		# Show plot
		if show:
			plt.show()

	# Output
	args = (fig, tinn_vals['tinn_n'], tinn_vals['tinn_m'], tinn_vals['tinn'], trindex)
	names = ('nn_histogram', 'tinn_n', 'tinn_m', 'tinn', 'tri_index')
	return utils.ReturnTuple(args, names)


def time_domain(signal=None,
				nn=None,
				rpeaks=None,
				sampling_rate=1000.,
				threshold=None,
				plot=True,
				show=False,
				binsize=7.8125,
				**kwargs):
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
	threshold : int
		Custom threshold in (ms) for the NNXX and pNNXX parameters.
	plot : bool
		If True, creates histogram plot using matplotlib, else uses numpy (data only, no plot) - (geometrical params).
	figsize : array_like, optional
		Matplotlib figure size for the histogram (width, height) (default: (6, 6)) - (geometrical params).
	binsize : int, float
		Bin size of the histogram bins - (geometrical params).
	legend : bool
		If True, highlights D(X) marker to the plot to be added to the legends (default=True) - (geometrical params).

	**kwargs
	--------
	full : bool, optional
		If True, returns last segment, even if the cumulative sum of NNI does not reach the 300s (default: False).
	overlap : bool, optional
		If True, allow to return NNI that go from the interval of one segment to the successive segment (default: False).
	duration : int, optional
		Maximum duration duration per segment in (s) (default: 300s).

	Returns
	-------
	results : biosppy.utils.ReturnTuple object
		All time domain results (see list and keys below)

	Returned Parameters
	-------------------
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
	..	nn20 & pNN20 (keys: 'nn20', 'pnn20')
	..	nnXX (XX = custom threshold) if specified (keys: 'nnXX', 'pnnXX')
	..	Triangular Index (key: 'tri_index')
	.. 	TINN (key: 'tinn', 'tinn_n', 'tinn_m')
	..	NNI histogram (key: 'nn_histogram')

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
	..	Default bin size set to recommended bin size of 1/128 (with 128Hz being the minimum recommended sampling
		frequency) as recommended by the HRV guidelines REFERENCE.
	..	'show' has only effect if 'plot' is also True.
	.. 	'legend' has only effect if 'plot' is also True.
	..	'figsize' has only effect if 'plot' is also True.

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

	# Call time domain functions & wrap results in a single biosspy.utils.ReturnTuple object
	results = nn_parameters(nn)
	results = tools.join_tuples(results, hr_parameters(nn))
	results = tools.join_tuples(results, nn_differences_parameters(nn))
	results = tools.join_tuples(results, sdnn(nn))
	results = tools.join_tuples(results, sdnn_index(nn, **kwargs))
	results = tools.join_tuples(results, sdann(nn, **kwargs))
	results = tools.join_tuples(results, rmssd(nn))
	results = tools.join_tuples(results, sdsd(nn))
	results = tools.join_tuples(results, nn50(nn))
	results = tools.join_tuples(results, nn20(nn))

	# Compute custom threshold if required
	if threshold is not None:
		results = tools.join_tuples(results, nnXX(nn, threshold=int(threshold)))

	# Compute geometrical parameters
	results = tools.join_tuples(results, geometrical_parameters(nn, plot=plot, show=show, binsize=binsize))

	# Output
	return results

if __name__ == "__main__":
	"""
	Example Script - HRV Time Domain Analysis
	"""
	from opensignalsreader import OpenSignalsReader
	from tools import nn_intervals

	# Load OpenSignals (r)evolution ECG sample file
	acq = OpenSignalsReader('./samples/SampleECG.txt')
	signal = acq.signal('ECG')

	# Filter data & get r-peak locations (ms)
	signal, rpeaks = ecg(signal, sampling_rate=acq.sampling_rate, show=False)[1:3]
	nni = nn_intervals(rpeaks)

	# Time Domain results
	print("=========================")
	print("TIME DOMAIN Results")
	print("=========================")

	hr_ = hr_parameters(nni)
	print("HR Results")
	print("> Mean HR:			%f (bpm)" % hr_['hr_mean'])
	print("> Min HR:			%f (bpm)" % hr_['hr_min'])
	print("> Max HR:			%f (bpm)" % hr_['hr_max'])
	print("> Std. Dev. HR:		%f (bpm)" % hr_['hr_std'])

	nn_para_ = nn_parameters(nni)
	print("NN Results")
	print("> Mean NN:			%f (ms)" % nn_para_['nn_mean'])
	print("> Min NN:			%f (ms)" % nn_para_['nn_min'])
	print("> Max NN:			%f (ms)" % nn_para_['nn_max'])

	nn_diff_ = nn_differences_parameters(nni)
	print("∆NN Results")
	print("> Mean ∆NN:			%f (ms)" % nn_diff_['nn_diff_mean'])
	print("> Min ∆NN:			%f (ms)" % nn_diff_['nn_diff_min'])
	print("> Max ∆NN:			%f (ms)" % nn_diff_['nn_diff_max'])
	print("> Std. Dev. ∆NN:	%f (ms)" % nn_diff_['nn_diff_std'])

	print("SDNN:				%f (ms)" % sdnn(nni)['sdnn'])
	print("SDNN Index:			%f (ms)" % sdnn_index(nni)['sdnn_index'])
	print("SDANN:				%f (ms)" % sdann(nni)['sdann'])
	print("RMMSD:				%f (ms)" % rmssd(nni)['rmssd'])
	print("SDSD:				%f (ms)" % sdsd(nni)['sdsd'])
	print("NN50:				%i (-)" % nn50(nni)['nn50'])
	print("pNN50: 				%f (%%)" % nn50(nni)['pnn50'])
	print("NN20:				%i (-)" % nn20(nni)['nn20'])
	print("pNN20: 				%f (%%)" % nn20(nni)['pnn20'])

	# Compute geometrical parameters (without plot)
	print("=== Geometrical Parameters")
	geo = geometrical_parameters(nni, plot=True, show=True)
	print("Triangular Index: 	%f (-)" % geo['tri_index'])
	print("TINN:				%f (ms)" % geo['tinn'])
	print("> N:				%f (ms)" % geo['tinn_n'])
	print("> M:				%f (ms)" % geo['tinn_m'])

	# Alternatively use the individual geometrical parameter functions
	geo = triangular_index(nni, plot=False)
	geo = tinn(nni, plot=False)

	# Alternatively use the time_domain() function to compute all time domain parameters using a single function
	time_domain(signal=signal)
