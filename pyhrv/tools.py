#!/usr/bin/env python -W ignore::FutureWarning
# -*- coding: utf-8 -*-
"""
pyHRV - Heart Rate Variability Toolbox - Tools
----------------------------------------------

This module provides basic tools for HRV analysis such as the computation of NN intervals, NN interval differences,
heart rate series and other functions (ecg plotting, tachogram, signal segmentation, hrv report, hrv import & export).

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
..  Hugo Silva, PhD, Instituto de Telecomunicacoes, PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

Docs
----
..	You can find the documentation for this module here:
	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html

Last Update
-----------
12-09-2018

:copyright: (c) 2018 by Pedro Gomes (HAW Hamburg)
:license: BSD 3-clause, see LICENSE for more details.
"""
from __future__ import absolute_import, division

# Version
from pyhrv.__version__ import __version__

import os
import warnings
import json
import collections
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime as dt
import biosppy
from biosppy import utils


# Turn off toolbox triggered warnings
WARN = False
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)


def nn_intervals(rpeaks=None):
	"""Computes the NN intervals [ms] between successive R-peaks.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#nn-intervals-nn-intervals

	Parameter
	---------
	rpeaks : array
		R-peak times in [ms] or [s]

	Returns
	-------
	nni : array
		NN intervals in [ms]

	Raises
	------
	TypeError
		If no data provided for 'rpeaks'
	TypeError
		If data format is not list or numpy array
	TypeError
		If 'rpeaks' array contains non-integer or non-float value

	Notes
	-----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#nn-intervals-nn-intervals

	"""
	# Check input signal
	if rpeaks is None:
		raise TypeError("No data for R-peak locations provided. Please specify input data.")
	elif type(rpeaks) is not list and not np.ndarray:
		raise TypeError("List, tuple or numpy array expected, received  %s" % type(rpeaks))

	if all(isinstance(n, int) for n in rpeaks) is False and all(isinstance(n, float) for n in rpeaks) is False:
		raise TypeError("Incompatible data type in list or numpy array detected (only int or float allowed).")

	# Confirm numpy arrays & compute NN intervals
	rpeaks = np.asarray(rpeaks)
	nn_int = np.zeros(rpeaks.size - 1)

	for i in range(nn_int.size):
		nn_int[i] = rpeaks[i + 1] - rpeaks[i]

	return nn_format(nn_int)


def nn_format(nni=None):
	"""Checks format of the NN intervals (seconds or milliseconds) and converts s data to ms data, if necessary.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#nn-format-nn-format

	Parameters
	----------
	nni : array
		Series of NN intervals in [ms] or [s]

	Returns
	-------
	nni : array
		Series of NN intervals in [ms]

	Raises
	------
	TypeError
		If no data provided for 'nni'

	Notes
	-----
	..	You can find the documentation for this module here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#nn-format-nn-format

	"""
	# Check input
	if nni is None:
		raise TypeError("No input data provided for 'nn'. Please specify input data")
	nn_ = np.asarray(nni, dtype='float64')

	# Convert if data has been identified in [s], else proceed with ensuring the NumPy array format
	if np.max(nn_) < 10:
		nn_ = [int(x * 1000) for x in nn_]

	return np.asarray(nn_)


def nni_diff(nni=None):
	"""Computes the series of differences between successive NN intervals [ms].

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#nn-interval-differences-nn-diff

	Parameters
	----------
	nni : array
		NN intervals in [ms] or [s].

	Returns
	-------
	nni_diff_ : numpy array
		Difference between successive NN intervals in [ms].

	Raises
	------
	TypeError
		If no data provided for 'rpeaks'.
	TypeError
		If no list or numpy array is provided.
	TypeError
		If NN interval array contains non-integer or non-float value.

	Notes
	..	You can find the documentation for this module here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#nn-interval-differences-nn-diff

	"""
	# Check input signal
	if nni is None:
		raise TypeError("No data for R-peak locations provided. Please specify input data.")
	elif type(nni) is not list and type(nni) is not np.ndarray:
		raise TypeError("List or numpy array expected, received  %s" % type(nni))
	elif all(isinstance(x, int) for x in nni) and all(isinstance(x, float) for x in nni):
		raise TypeError("'nni' data contains non-int or non-float data.")
	else:
		nn = nn_format(nni)

	# Confirm numpy arrays & compute NN interval differences
	nn_diff_ = np.zeros(nn.size - 1)

	for i in range(nn.size - 1):
		nn_diff_[i] = abs(nn[i + 1] - nn[i])

	return np.asarray(nn_diff_)


def heart_rate(nni=None, rpeaks=None):
	"""Computes a series of Heart Rate values in [bpm] from a series of NN intervals or R-peaks in [ms] or [s] or the HR from a single NNI.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#heart-rate-heart-rate

	Parameters
	----------
	nni : int, float, array
		NN intervals in [ms] or [s].
	rpeaks : int, float, array
		R-peak times in [ms] or [s].

	Returns
	-------
	bpm : list, numpy array, float
		Heart rate computation [bpm].
		Float value if 1 NN interval has been provided
		Float array if series of NN intervals or R-peaks are provided.

	Raises
	------
	TypeError
		If no input data for 'rpeaks' or 'nn_intervals provided.
	TypeError
		If provided NN data is not provided in float, int, list or numpy array format.

	Notes
	-----
	..	You can find the documentation for this module here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#heart-rate-heart-rate

	"""
	# Check input
	if nni is None and rpeaks is not None:
		# Compute NN intervals if rpeaks array is given; only 1 interval if 2 r-peaks provided
		nni = nn_intervals(rpeaks) if len(rpeaks) > 2 else int(np.abs(rpeaks[1] - rpeaks[0]))
	elif nni is not None:
		# Use given NN intervals & confirm numpy if series of NN intervals is provided
		if type(nni) is list or type(nni) is np.ndarray:
			nni = nn_format(nni) if len(nni) > 1 else nni[0]
		elif type(nni) is int or float:
			nni = int(nni) if nni > 10 else int(nni) / 1000
	else:
		raise TypeError("No data for R-peak locations or NN intervals provided. Please specify input data.")

	# Compute heart rate data
	if type(nni) is int:
		return 60000. / float(nni)
	elif type(nni) is np.ndarray:
		return np.asarray([60000. / float(x) for x in nni])
	else:
		raise TypeError("Invalid data type. Please provide data in int, float, list or numpy array format.")


def plot_ecg(signal=None,
			 t=None,
			 sampling_rate=1000.,
			 interval=None,
			 rpeaks=True,
			 figsize=None,
			 title=None,
			 show=True):
	"""Plots ECG signal on a medical grade ECG paper-like figure layout.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#plot-ecg-plot-ecg

	Parameters
	----------
	signal : array
		ECG signal (filtered or unfiltered)
	t : array, optional
		Time vector for the ECG signal (default: None)
	sampling_rate : int, float, optional
		Sampling rate of the acquired signal in [Hz] (default: 1000Hz)
	interval : array, 2-element, optional
		Visualization interval of the ECG signal plot (default: None: [0s, 10s]
	rpeaks : bool, optional
		If True, marks R-peaks in ECG signal (default: True)
	figsize : array, optional
		Matplotlib figure size (width, height) (default: None: (12, 4))
	title : str, optional
		Plot figure title (default: None).
	show : bool, optional
		If True, shows the ECG plot figure(default: True)

	Returns
	-------
	fig_ecg : matplotlib figure object
		Matplotlib figure of ECG plot

	Raises
	------
	TypeError
		If no ECG data provided.

	Notes
	----
	..	The 'rpeaks' parameter will have no effect if there are more then 50 r-epaks within the visualization interval.
		In this case, no markers will be set to avoid overloading the plot
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#plot-ecg-plot-ecg

	"""
	# Check input data
	if signal is None:
		raise TypeError("No ECG data provided. Please specify input data.")
	else:
		# Confirm numpy
		signal = np.asarray(signal)

	# Compute time vector
	if t is None:
		t = time_vector(signal, sampling_rate=sampling_rate)

	# Configure interval of visualized signal
	interval = check_interval(interval, limits=[0, t[-1]], default=[0, 10])

	# Prepare figure
	if figsize is None:
		figsize = (12, 4)

	fig_ecg = plt.figure(figsize=figsize)
	ax = fig_ecg.add_subplot(111)

	# Configure axis according to according to BITalino ECG sensor ranges
	if signal.max() > 1.5:
		y_min = int(signal.min() - (signal.max() - signal.min()) * 0.2)
		y_max = int(signal.max() + (signal.max() - signal.min()) * 0.2)
		unit = '-'
		y_minor = np.linspace(y_min, y_max, 12)
		y_major = np.linspace(y_min, y_max, 4)
	elif signal.max() < 1.0:
		y_min, y_max = -1., 1.,
		unit = 'mV'
		y_minor = np.arange(-0.9, y_min, 0.1)
		y_major = np.arange(-1.0, y_max + 0.5, 0.5)
	else:
		y_min, y_max = -1.5, 1.5,
		unit = 'mV'
		y_minor = np.arange(-1.4, y_min, 0.1)
		y_major = np.arange(y_min, y_max + 0.5, 0.5)

	ax.axis([interval[0], interval[1], y_min, y_max])
	ax.set_xlabel('Time [$s$]')
	ax.set_ylabel('ECG [$%s$]' % unit)

	# Set ticks as ECG paper (box height ~= 0.1mV; width ~= 0.1s when using default values)
	n = int(interval[1] / 10)
	ax.set_xticks(np.arange(0.0, interval[1] + 0.1, float(n)/5), minor=True)
	ax.xaxis.grid(which='minor', color='salmon', lw=0.3)
	ax.set_xticks(np.arange(0, interval[1] + 0.1, n))
	ax.xaxis.grid(which='major', color='r', lw=0.7)
	ax.set_yticks(y_minor, minor=True)
	ax.yaxis.grid(which='minor', color='salmon', lw=0.3)
	ax.set_yticks(y_major)
	ax.yaxis.grid(which='major', color='r', lw=0.7)

	# Add legend
	unit = '' if unit == '-' else unit
	text_ = 'Division (x): %is\nDivision (y): %.1f%s' % (n, (np.abs(y_major[1] - y_major[0])), unit)
	ax.text(0.88, 0.85, text_, transform=ax.transAxes, fontsize=9,
		bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

	# Plot ECG signal
	ax.plot(t, signal, 'r')
	fig_ecg.tight_layout()

	# Plot r-peaks
	rps = biosppy.signals.ecg.ecg(signal=signal, sampling_rate=sampling_rate, show=False)[2]
	p = [float(signal[x]) for x in rps]
	r = t[rps]
	if rpeaks:
		ax.plot(r, p, 'g*', alpha=0.7)

	# Add title
	if title is not None:
		ax.set_title('ECG Signal - %s' % str(title))
	else:
		ax.set_title('ECG Signal')

	# Show plot
	if show:
		plt.show()

	# Output
	args = (fig_ecg, )
	names = ('ecg_plot', )
	return utils.ReturnTuple(args, names)


def tachogram(nni=None,
			  signal=None,
			  rpeaks=None,
			  sampling_rate=1000.,
			  hr=True,
			  interval=None,
			  title=None,
			  figsize=None,
			  show=True):
	"""Plots Tachogram (NNI & HR) of an ECG signal, NNI or R-peak series.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#tachogram-tachogram

	Parameters
	----------
	nni : array
		NN intervals in [ms] or [s].
	rpeaks : array
		R-peak times in [ms] or [s].
	signal : array, optional
		ECG signal.
	sampling_rate : int, float
		Sampling rate of the acquired signal in [Hz].
	hr : bool, optional
		If True, plots series of heart rate data in [bpm] (default: True).
	interval : list, optional
		Sets visualization interval of the signal (default: [0, 10]).
	title : str, optional
		Plot figure title (default: None).
	figsize : array, optional
		Matplotlib figure size (width, height) (default: (12, 4)).
	show : bool, optional
		If True, shows plot figure (default: True).

	Returns
	-------
	fig : matplotlib.pyplot figure
		Tachogram figure & graph

	Raises
	------
	TypeError
		If no input data for 'nni', 'rpeaks' or 'signal' is provided

	Notes
	-----
	..	NN intervals are derived from the ECG signal if 'signal' is provided.
	.. 	If both 'nni' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nni' data will be computed
		from the 'rpeaks'.
	..	If both 'nni' and 'signal' are provided, 'nni' will be chosen over 'signal'.
	..	If both 'rpeaks' and 'signal' are provided, 'rpeaks' will be chosen over 'signal'.

	"""
	# Check input
	if signal is not None:
		rpeaks = biosppy.signals.ecg.ecg(signal=signal, sampling_rate=sampling_rate, show=False)[2]
	elif nni is None and rpeaks is None:
		raise TypeError('No input data provided. Please specify input data.')

	# Get NNI series
	nni = check_input(nni, rpeaks)

	# Time vector back to ms
	t = np.cumsum(nni) / 1000.

	# Configure interval of visualized signal
	interval = check_interval(interval, limits=[0, t[-1]], default=[0, 10])

	# Prepare figure
	if figsize is None:
		figsize = (12, 4)
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(111)

	# X-Axis configuration
	# Set x-axis format to seconds if the duration of the signal <= 60s
	if interval[1] <= 60:
		ax.set_xlabel('Time [s]')
	# Set x-axis format to MM:SS if the duration of the signal > 60s and <= 1h
	elif 60 < interval[1] <= 3600:
		ax.set_xlabel('Time [MM:SS]')
		formatter = mpl.ticker.FuncFormatter(lambda ms, x: str(dt.timedelta(seconds=ms))[2:])
		ax.xaxis.set_major_formatter(formatter)
	# Set x-axis format to HH:MM:SS if the duration of the signal > 1h
	else:
		ax.set_xlabel('Time [HH:MM:SS]')
		formatter = mpl.ticker.FuncFormatter(lambda ms, x: str(dt.timedelta(seconds=ms)))
		ax.xaxis.set_major_formatter(formatter)

	n = int(interval[1] / 10)
	ax.set_xticks(np.arange(0, interval[1] + n, n))

	# Y-Axis configuration (min, max set to maximum of the visualization interval)
	ax.set_ylabel('NN Interval [$ms$]')
	nn_min = np.min(nni[np.argwhere(np.logical_and(interval[0] <= t, t <= interval[1]))])
	nn_max = np.max(nni[np.argwhere(np.logical_and(interval[0] <= t, t <= interval[1]))])
	ax.axis([interval[0], interval[1], nn_min * 0.9, nn_max * 1.1])

	# Plot 'x' markers only if less than 50 rpeaks are within the given data, otherwise don't add them
	if np.argwhere(t < interval[1]).size < 50:
		l1 = ax.plot(t, nni, color='g', label='NN Intervals', marker='x', linestyle='--', linewidth=0.8)
		ax.vlines(t, 200, 3000, linestyles='--', linewidth=0.5, alpha=0.7, colors='lightskyblue')
	else:
		l1 = ax.plot(t, nni, color='g', label='NN Intervals', linestyle='--', linewidth=0.8)
	lns = []

	# Plot heart rate signal
	if hr:
		ax2 = ax.twinx()
		bpm_values = heart_rate(nni)
		hr_min = heart_rate(nn_max)
		hr_max = heart_rate(nn_min)

		ax2.set_ylabel('Heart Rate [$1/min$]', rotation=270, labelpad=15)
		ax2.axis([interval[0], interval[1], hr_min * 0.9, hr_max * 1.1])

		# Plot 'x' markers only if less than 50 rpeaks are within the given data, otherwise don't add them
		if np.argwhere(t < interval[1]).size < 50:
			l2 = ax2.plot(t, bpm_values, color='red', label='Heart Rate', marker='x', linestyle='--', linewidth=0.8)
		else:
			l2 = ax2.plot(t, bpm_values, color='red', label='Heart Rate', linestyle='--', linewidth=0.8)
		lns = l1 + l2
		labs = [l.get_label() for l in lns]
		ax.legend(lns, labs, loc=1)
	else:
		ax.legend(loc=1)

	# Add title
	if title is not None:
		ax.set_title('Tachogram - %s' % str(title))
	else:
		ax.set_title('Tachogram')

	# Show plot
	if show:
		plt.show()

	# Output
	args = (fig,)
	names = ('tachogram_plot',)
	return utils.ReturnTuple(args, names)


def check_interval(interval=None, limits=None, default=None):
	"""General purpose function that checks and verifies correctness of interval limits within optionally defined
	valid interval specifications and and/or default values if no interval is specified.

	This function can be used to set visualization intervals, check overlapping frequency bands, or for other similar
	purposes and is intended to automatically catch possible error sources due to invalid intervals.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#check-interval-check-interval

	Parameters
	----------
	interval : array, 2-elements
		Input interval [min, max] (default: None)
	limits : array, 2-elements
		Minimum and maximum allowed interval limits (default: None)
	default : array, 2-elements
		Specified default interval (e.g. if 'interval' is None)

	Returns
	-------
	interval : array
		Valid interval limits.

	Raises
	------
	TypeError
		If no input data is specified
	ValueError
		If any of the input data has equal lower and upper limits

	Warnings
	--------
	..	If overlapping limits had to be fixed.
	..	If limits are overlapping (e.g. lower limit > upper limit)

	Notes
	-----
	(Warnings are only triggered if the module variable 'WARN' is set to 'True')
	..	If 'interval[0]' is out of boundaries ('interval[0]' < 'limit[0]' or 'interval[0]' >= 'limit[1]')
		-> sets 'interval[0] = limit[0]'
	..	If 'interval[1]' is out of boundaries ('interval[1]' <= 'limit[0]' or 'interval[1]' > 'limit[1]')
		-> sets 'interval[1] = limit[1]'
	..	If limits are overlapping (e.g. lower limit > upper limit) and had to be fixed.
	..	This thing is got unnecessarily complicated, but in the end I just went with it...
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#check-interval-check-interval

	"""
	if interval is None and limits is None and default is None:
		raise TypeError("No input data specified. Please verify your input data.")

	# Create local copy to prevent changing input variable
	interval = list(interval) if interval is not None else None

	# Check default limits
	if default is not None:
		default = _check_limits(default, 'default')

	# Check maximum range limits
	if limits is None and default is not None:
		limits = default
	elif limits is not None:
		limits = _check_limits(limits, 'limits')

	# Check interval limits
	if interval is None:
		if default is not None:
			return default
		elif default is None and limits is not None:
			return limits

	# If only interval is specified, but not 'min', 'max' or 'default' check if lower limit >= upper limit
	elif interval is not None and limits is None:
		return _check_limits(interval, 'interval')

	# If none of the input is 'None'
	else:
		# Check interval
		interval = _check_limits(interval, 'interval')
		if not limits[0] <= interval[0]:
			interval[0] = limits[0]
			if WARN:
				warnings.warn("Interval limits out of boundaries. Interval set to: %s" % interval, stacklevel=2)
		if not limits[1] >= interval[1]:
			interval[1] = limits[1]
			if WARN:
				warnings.warn("Interval limits out of boundaries. Interval set to: %s" % interval, stacklevel=2)
		return interval


def _check_limits(interval, name):
	"""Checks if interval limits are not overlapping or equal.

	Parameters
	----------
	interval : array, 2-element
		Interval boundaries [min, max].
	name : str
		Variable name to be used on exceptions and warnings.

	Returns
	-------
	interval : array, 2-element
		Valid and corrected interval limits (if correction is necessary).

	Raises
	------
	ValueError
		If interval limits are equal.

	Warnings
	--------
	..	If limits are overlapping (e.g. lower limit > upper limit) and had to be fixed.

	"""
	# upper limit < 0 or upper limit > max interval -> set upper limit to max
	if interval[0] > interval[1]:
		interval[0], interval[1] = interval[1], interval[0]
		vals = (name, name, interval[0], interval[1])
		if WARN:
			warnings.warn("Corrected invalid '%s' limits (lower limit > upper limit).'%s' set to: %s" % vals)
	if interval[0] == interval[1]:
		raise ValueError("'%f': Invalid interval limits as they are equal." % name)
	return interval


def segmentation(nni=None, full=True, overlap=False, duration=300):
	"""Segmentation of signal peaks into individual segments of set duration.
	(e.g. splitting R-peak locations into 5min segments for computation of SDNN index)

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#segmentation-segmentation

	Parameters
	----------
	nni: array
		Series of NN intervals [ms] or [s].
	full : bool, optional
		If True, returns last segment, even if the sum of NNI does not reach the segment duration (default: False).
	overlap : bool, optional
		If True, allow to return NNI that go from the interval of one segment to the successive segment (default: False).
	duration : int, optional
		Maximum duration duration per segment in [s] (default: 300s).

	Returns
	-------
	segments : array, array of arrays
		NN intervals in each segment/time interval. If cumulative sum of NN input data < duration, the NN input data
		will be returned.

	Raises
	------
	TypeError
		If no 'nni' data is not provided.

	Warnings
	--------
	If signal is shorter than the specified duration.

	Notes
	----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#segmentation-segmentation

	"""
	# Check input
	if nni is None:
		raise TypeError("Please specify input data.")

	# Preparations
	nni= nn_format(nni)
	tn = np.cumsum(nni)
	max_time = tn[-1]
	duration *= 1000			# convert from s to ms

	# Check if signal is longer than maximum segment duration
	if np.sum(nni) > duration:

		# Compute limits for each segment
		segments = []
		limits = np.arange(0, max_time + duration, duration)

		# Current index
		cindex = 0

		# Segment signals
		for i, _ in enumerate(range(0, limits.size - 1)):
			csegment = []
			while np.sum(csegment) < duration:
				csegment.append(nni[cindex])
				if cindex < nni.size - 1:
					cindex += 1
				else:
					break

			# Check if overlap exists (just to be sure)
			if np.sum(csegment) > duration:
				csegment = csegment[:-1]
				cindex -= 1

			segments.append(list(csegment))

		# Remove the last incomplete segment if required
		if not full:
			segments = segments[:-1]

		return segments, True
	else:
		if WARN:
			warnings.warn("Signal duration is to short for segmentation into %is. "
						  "Input data will be returned." % duration)
		return nni, False


def hrv_report(results=None, path=None, rfile=None, nni=None, info={}, file_format='txt', delimiter=';', hide=False, plots=False):
	"""Generates HRV report (in .txt or .csv format) of the provided HRV results.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#hrv-reports-hrv-report

	Parameters
	----------
	results : dict, biosppy.utils.ReturnTuple object
		Results of the HRV parameter computations
	rfile : str, file handler
		Absolute path of the output directory or report file handler
	nni : array, optional
		NN interval series in [ms] or [s]
	info : dict, optional
		Dictionary with HRV metadata:
		---------------------------------------------------------------
		Keys		:	Description
		---------------------------------------------------------------
		file		:	Name of the signal acquisition file.
		device		:	Acquisition device.
		identifier	:	ID of the acquisition device (e.g. MAC-address)
		fs			:	Sampling rate used during acquisition.
		resolution	:	Resolution used during acquisition.
		---------------------------------------------------------------
	file_format : str, optional
		Output file format, select 'txt' or 'csv' (default: 'txt')
	delimiter : str, optional
		Delimiter to separate the columns in the exported report (default: ';')
	hide : bool
		Hide parameters in report that have not been computed
	plots : bool, optional
		If True, save figures in results as .png file

	Returns
	-------
	rfile : str
		Absolute path of the output report file (may vary from the input file name)

	Raises
	------
	TypeError
		If no HRV results provided
	TypeError
		If no file or directory path provided
	TypeError
		If selected output file format is not supported
	IOError
		If selected file or directory does not exist

	Warnings
	--------
	..	New generated file name, if provided file does already exist

	Notes
	-----
	..	Existing files will not be overwritten, instead the new file will consist of the given file name with an
		(incremented) identifier (e.g. '_1') that will be added at the end of the provided file name.
	..	Plot figures will be saved in .png format.
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#hrv-reports-hrv-report

	"""
	SUPPORTED_FILE_FORMATS = ['txt', 'csv']
	plot_files = []

	# Check input
	if results is None:
		raise TypeError("No HRV results provided. Please specify input data.")
	if file_format not in SUPPORTED_FILE_FORMATS:
		raise TypeError("Unsupported file format. Only txt and csv formats are supported.")

	# Check path input
	if path is None:
		raise TypeError("No file name or directory provided. Please specify at least an output directory.")
	elif type(path) is str:
		if rfile is None:
			# Generate automatic file name
			rfile = 'hrv_report' + dt.datetime.now().strftime('_%Y-%m-%d_%H-%M-%S') + '.' + file_format
			path += rfile
		else:
			# Check if file name has an compatible extension
			_, fformat = os.path.splitext(rfile)
			if fformat != file_format or fformat == '':
				path = path + rfile + '.' + file_format
			else:
				path = path + rfile
	elif type(path) is file:
		path_ = path.name
		path.close()
		path = path_

	if 'hrv_export' in str(path):
		path = path.replace('hrv_export', 'hrv_report')

	rfile, _ = _check_fname(path, file_format, rfile)

	# Load HRV parameters metadata
	params = json.load(open(os.path.join(os.path.split(__file__)[0], './files/hrv_keys.json'), 'r'), encoding='utf-8')

	# Prepare output dictionary
	output = {'Name': rfile}
	if 'comment' in results.keys():
		output['comment'] = results['comment']
	else:
		output['comment'] = 'n/a'

	for key in results.keys():
		if not isinstance(results[key], plt.Figure):
			# Hide 'nan' and 'n/a' values in report if preferred
			if hide and str(results[key]) in ['nan', 'n/a']:
				continue
			else:
				output[key] = results[key]
		else:
			# If matplotlib figure found, plots=True, and the figure is known in 'hrv_keys.json' -> save the figure
			if plots:
				if key in params:
					pfname = os.path.splitext(rfile)[0] + '_' + str(key)
					results[key].savefig(pfname, dpi=300)
					plot_files.append(os.path.split(pfname)[1] + '.png')

	# Metadata
	tstamp = dt.datetime.now().strftime('_%Y-%m-%d_%H-%M-%S')
	mdesc = {'file': 'File', 'device': 'Device',
			'identifier': 'Identifier/MAC', 'fs': 'Sampling Rate', 'resolution': 'Resolution', 'tstamp': 'Date Time'}
	mdata = {}

	for key in mdesc.keys():
		mdata[key] = info[key] if key in info.keys() else 'n/a'
	mdata['tstamp'] = tstamp[1:]

	# Prepare text file format
	hformat = '# %s %*s %s\n'
	cformat = '%s %*s %*s\n'
	sformat = '\n	-	%s'
	titles = {
		'time': 'Time Domain',
		'frequency_fft': 'Frequency Domain - FFT Welch\'s Method',
		'frequency_ar': 'Frequency Domain - Autoregression Method',
		'frequency_lomb': 'Frequency Domain - Lomb-Scargle Method',
		'nonlinear': 'Nonlinear Methods',
		'metadata': 'METADATA',
		'plots': 'Plots'
	}

	with open(rfile, 'w+') as f:
		# Write header
		f.write('# BIOSPPY HRV REPORT - v.%s\n' % __version__)
		for key in mdata.keys():
			f.write(hformat % (mdesc[key], 20 - len(mdesc[key]), delimiter, mdata[key]))

		# Add comment
		line = 70
		f.write('\n\n%s\n%s\n%s\n%s\n' % (line * '-', 'COMMENTS', output['comment'], line * '-'))

		# Prepare content
		content = dict()
		for key in params.keys():
			if 'plot' not in str(key):
				key = str(key)

				# Prepare parameter label & units
				if key not in ['comment', 'Name']:
					label = str(params[key][1]) + (' (%s)' % params[key][2])

					# Select output parameters (to report file)
					if key in output.keys() and 'plot' not in str(key) and 'histogram' not in str(key):
						if str(output[key]) in ['nan', 'n/a'] and hide:
							continue
						else:
							para = output[key]
							out = ''
							if isinstance(para, collections.Iterable) and type(para) is not str:
								for i, val in enumerate(list(para)):
									if val is str or np.nan:
										val_ = str(val) if val not in ['n/a', 'nan'] else 'n/a'
									elif val == float:
										val_ = '%.3f' % val if val not in ['n/a', 'nan'] else 'n/a'
									out += sformat % val_
							elif type(para) is str:
								out = '%s' % output[key] if output[key] not in ['n/a', 'nan', None] else 'n/a'
							else:
								out = '%.3f' % output[key] if output[key] not in ['n/a', 'nan', None] else 'n/a'
						content[key] = cformat % (label, 50 - len(label), delimiter, 1, out)
					elif not hide:
						content[key] = cformat % (label, 50-len(label), delimiter, 1, 'n/a')
				else:
					content['comment'] = output[key]

		# Write output to report file
		current_domain = []

		# Go through all parameters in 'hrv_keys.json'
		for n in range(1, len(params.keys()) - 1):
			# Go through content
			for key in content.keys():
				# Check keys by specified order set in the last element of each entry in 'hrv_keys.json'
				if params[key][-1] == n:
					# Set parameter title in output file (Time Domain, Frequency Domain, etc.)
					if current_domain != params[key][0] and str(key) not in ['nn_intervals']:
						current_domain = params[key][0]
						if current_domain != 'plot':
							f.write('\n\n%s\n%s\n%s\n' % (line * '-', titles[current_domain], line * '-'))

					# Finally, write parameter content to output file
					f.write(content[key])

		# Add generated plot figures (file names to facilitate the link between report files and plot figures)
		if plots:
			f.write('\n\n%s\n%s\n%s\n' % (line * '-', 'Plot Figure Files', line * '-'))
			for plot in plot_files:
				f.write('%s\n' % plot)

		# Add NNI if desired
		if nni is not None:
			f.write('\n\n%s\n' % params['nn_intervals'][1])
			for i in nni:
				f.write('%.3f\n ' % i)

	return rfile


def hrv_export(results=None, path=None, efile=None, comment=None, plots=False):
	"""
	Exports HRV results into a JSON file.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#hrv-export-hrv-export

	Parameters
	----------
	results : dict, biosppy.utils.ReturnTuple object
		Results of the HRV analysis
	path : str
		Absolute path of the output directory
	efile : str, optional
		Output file name
	comment : str, optional
		Optional comment
	plots : bool, optional
		If True, save figures of the results in .png format

	Returns
	-------
	efile : str
		Absolute path of the output export file (may vary from the input data)

	Raises
	------
	TypeError
		No input data provided
	TypeError
		Unsupported data format provided (other than dict, or utils.ReturnTuple object.)
	TypeError
		If no file or directory path provided

	Notes
	-----
	..	If 'path' is a file handler, 'efile' will be ignored.
	..	Creates file with automatic name generation if only an output path is provided.
	..	Output file name may vary from input file name due changes made to avoid overwrting existing files (your
		results are important after all!).
	..	Existing files will not be overwritten, instead the new file will consist of the given file name with an
		(incremented) identifier (e.g. '_1') that will be added at the end of the provided file name.
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#hrv-export-hrv-export

	"""
	# Check input (if available & utils.ReturnTuple object)
	if results is None:
		raise TypeError("No results data provided. Please specify input data.")
	elif results is not type(dict()) and isinstance(results, utils.ReturnTuple) is False:
		raise TypeError("Unsupported data format: %s. "
						"Please provide input data as Python dictionary or biosppy.utils.ReturnTuple object." % type(results))

	if path is None:
		raise TypeError("No file name or directory provided. Please specify at least an output directory.")
	elif type(path) is str:
		if efile is None:
			# Generate automatic file name
			efile = 'hrv_export' + dt.datetime.now().strftime('_%Y-%m-%d_%H-%M-%S') + '.json'
			path += efile
		else:
			# Check if file name has an '.json' extension
			_, fformat = os.path.splitext(efile)
			if fformat != 'json':
				path = path + efile + '.json'
			else:
				path = path + efile
	elif type(path) is file:
		path_ = path.name
		path.close()
		path = path_

	efile, _ = _check_fname(path, 'json', efile)

	# Get HRV parameters
	params = json.load(open(os.path.join(os.path.split(__file__)[0], './files/hrv_keys.json'), 'r'))

	# Save plot figures
	if plots:
		for key in results.keys():
			if isinstance(results[key], plt.Figure) and key in params.keys():
				results[key].savefig(os.path.splitext(efile)[0] + '_' + str(key), dpi=300)

	# Prepare output dictionary
	output = {'Name': efile, 'Comment': str(comment)}
	for key in results.keys():
		if 'plot' not in str(key) and 'histogram' not in str(key):
			output[key] = results[key] if str(results[key]) != 'nan' else 'n/a'

	# Write output dictionary to JSON file
	json.encoder.FLOAT_REPR = lambda o: format(o, 'f')
	with open(efile, 'w+') as f:
		json.dump(output, f, sort_keys=True, indent=4, separators=(',', ': '))

	return str(efile)


def hrv_import(hrv_file=None):
	"""Imports HRV results stored in JSON files generated with the 'hrv_export()' function.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#hrv-import-hrv-import

	Parameters
	----------
	hrv_file : file object, str
		File handler or absolute string path of the HRV JSON file

	Returns
	-------
	output : biosppy.utils.ReturnTuple object
		All imported results.

	Raises
	------
	TypeError
		No input data provided.

	Notes
	-----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#hrv-import-hrv-import

	"""
	# Check input data and load JSON file content
	if hrv_file is None:
		raise TypeError("No input data provided. Please specify input data.")
	elif type(hrv_file) is str:
		data = json.load(open(hrv_file, 'r'))
	elif isinstance(hrv_file, file):
		data = json.load(hrv_file)

	results = dict()
	for key in data.keys():
		results[str(key)] = data[key] if type(data[key]) is not unicode else str(data[key])

	# Create utils.ReturnTuple object from imported data
	return utils.ReturnTuple(results.values(), results.keys())


def join_tuples(*args):
	"""Joins multiple biosppy.utils.ReturnTuple objects into one biosppy.utils.ReturnTuple object.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#join-tuples-join-tuples

	Parameters
	----------
	tuples : list, array, utils.ReturnTuple objects
		List or array containing utils.ReturnTuple objects.

	Returns
	-------
	return_tuple : biosppy.utils.ReturnTuple object
		biosppy.utils.ReturnTuple object with the content of all input tuples joined together.

	Raises
	------
	TypeError:
		If no input data is provided
	TypeError:
		If input data contains non-biosppy.utils.ReturnTuple objects

	Notes
	----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#join-tuples-join-tuples

	"""
	# Check input
	if args is None:
		raise TypeError("Please specify input data.")

	for i in args:
		if not isinstance(i, utils.ReturnTuple):
			raise TypeError("The list of tuples contains non-utils.ReturnTuple objects.")

	# Join tuples
	names = ()
	vals = ()

	for i in args:
		for key in i.keys():
			names = names + (key, )
			vals = vals + (i[key], )

	return utils.ReturnTuple(vals, names)


def _check_fname(rfile, fformat='txt', name='new_file'):
	"""Checks if file or path exists and creates new file with a given suffix or incrementing identifier at the end of
	the file name to avoid overwriting existing files.

	Parameters
	----------
	rfile : str
		Absolute file path or directory.
	fformat : str
		File format (e.g. 'txt').
	name : str
		File name for newly created file if only directory is provided.

	Returns
	-------
	rfile : str
		Absolute file path of a new file.

	Raises
	------
	IOError
		If file or directory does not exist.

	Warnings
	--------
	..	New generated file name, if provided file does already exist.

	Notes
	-----
	..	Existing files will not be overwritten, instead the new file will consist of the given file name with an
		(incremented) identifier (e.g. '_1') that will be added at the end of the provided file name.
	..	If only directory provided, the new file name will consist of the provided 'name' string.

	"""
	# Prepare path data
	fformat = '.' + fformat

	# Check if 'rfile' is a path name or path + file name
	if os.path.isdir(rfile) and not os.path.isfile(rfile + name + fformat):
		rfile = rfile + name + fformat
	elif os.path.isfile(rfile):
		old_rfile = rfile

		# Increment file identifier until an available file name has been found
		while(os.path.isfile(rfile)):
			rfile, format = os.path.splitext(rfile)

			# check for duplicate file name and create new file with (incremented) number at the end
			if rfile[-3:-1].isdigit():
				rfile = rfile[:-3] + str(int(rfile[-3:]) + 1)
			elif rfile[-2:-1].isdigit():
				rfile = rfile[:-2] + str(int(rfile[-2:]) + 1)
			elif rfile[-1].isdigit():
				rfile = rfile[:-1] + str(int(rfile[-1:]) + 1)
			elif rfile[-1].isdigit and rfile[-1] != '_':
				rfile += '_1'
			rfile += ('%s' % fformat)

		# Show warning if file does already exist
		msg = "\nFile '%s' does already exist." \
			  "New file name '%s' selected to avoid overwriting existing files." % (old_rfile, rfile)
		warnings.warn(msg, stacklevel=2)
	elif not os.path.isfile(rfile):
		rfile = os.path.splitext(rfile)[0] + fformat
		with open(rfile, 'w+'):
			pass
	else:
		raise IOError("File or directory does not exist. Please verify input data.")

	return rfile, os.path.split(rfile)


def std(array, dof=1):
	"""Computes the standard deviation of a data series.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#standard-deviation-std

	Parameters
	----------
	array : list, numpy array
		Data series.
	dof : int
		Degree of freedom (default to 1).

	Returns
	-------
	result : float
		Standard deviation of the input data series

	Raises
	------
	TypeError
		If input array is not specified

	Notes
	-----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#standard-deviation-std

	"""
	if array is None:
		raise TypeError("Please specify input data.")

	array = np.asarray(array)
	result = np.sum([(x - np.mean(array))**2 for x in array])
	result = np.sqrt(1. / (array.size - dof) * result)
	return result


def time_vector(signal=None, sampling_rate=1000.):
	"""Computes time vector for the provided input signal.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#time-vector-time-vector

	Parameters
	----------
	signal : array
		ECG signal (or any other sensor signal)
	sampling_rate : int, float, optional
		Sampling rate of the input signal in [Hz]

	Returns
	-------
	time_vector : array
		Time vector for the input signal sampled at the input 'sampling_rate'

	Raises
	------
	TypeError
		If input signal is not specified.

	Notes
	-----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#time-vector-time-vector

	"""
	if signal is None:
		raise TypeError("Please specify input signal.")

	signal = np.asarray(signal)
	t = np.arange(0, signal.size / sampling_rate, 1./sampling_rate)
	return t


def check_input(nni=None, rpeaks=None):
	"""Checks if input series of NN intervals or R-peaks are provided and, if yes, returns a NN interval series in [ms]
	format.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#check-input-check-input

	Parameters
	----------
	nni : array, int
		NN intervals in [ms] or [s].
	rpeaks : array, int
		R-peak times in [ms] or [s].

	Returns
	-------
	nni : array
		NN interval series [ms].

	Raises
	------
	TypeError
		If no input data for 'nni' and 'rpeaks' provided.

	Notes
	-----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#check-input-check-input

	"""
	# Check input
	if nni is None and rpeaks is not None:
		# Compute NN intervals if r_peaks array is given
		nni = nn_intervals(rpeaks)
	elif nni is not None:
		# Use given NN intervals & confirm numpy
		nni = nn_format(nni)
	else:
		raise TypeError("No R-peak data or NN intervals provided. Please specify input data.")
	return nni


if __name__ == "__main__":
	"""
	Example Script - HRV Tools
	"""
	import pyhrv.time_domain as td
	from opensignalsreader import OpenSignalsReader
	from biosppy.signals.ecg import ecg

	# Change path to your preferred directory
	path = './files/'

	# Load OpenSignals (r)evolution ECG sample file
	acq = OpenSignalsReader(path + 'SampleECG.txt')
	signal = acq.signal('ECG')
	sampling_rate = acq.sampling_rate

	# Filter data & get r-peak locations [ms]
	signal, rpeaks = ecg(signal, sampling_rate=sampling_rate, show=False)[1:3]

	# Plot ECG & save plot figure
	results = plot_ecg(signal, sampling_rate=sampling_rate, show=True, interval=[0, 22])

	# Plot Tachogram & save plot figure
	results = join_tuples(results, tachogram(rpeaks=rpeaks, show=False))

	# Compute Time Domain Parameters
	results = join_tuples(results, td.time_domain(rpeaks=rpeaks, plot=True, show=False))

	# Export HRV results to JSON file & save plot figures
	hrv_export(results, path=path, efile='SampleExport3', plots=True)

	# Create HRV report file & save plot figures
	hrv_report(results, path=path, rfile='SampleReport3', plots=True)

	# Import HRV results from exported JSON export file
	hrv_import(path + 'SampleExport.json')

	check_interval(None, default=[2, 8])
