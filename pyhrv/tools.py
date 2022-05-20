#!/usr/bin/env python -W ignore::FutureWarning
# -*- coding: utf-8 -*-
"""
pyHRV - Heart Rate Variability Toolbox - Tools
----------------------------------------------

This module provides support tools for HRV analysis such as the computation of HRV relevant data series (NNI, NNI
differences Heart Rate) and

Notes
-----
..  This module is part of the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..	This module is a contribution to the open-source biosignal processing toolbox 'BioSppy':
	https://github.com/PIA-Group/BioSPPy

Author
------
..  Pedro Gomes, pgomes92@gmail.com

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
12-11-2019

:copyright: (c) 2018 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.

"""
# Compatibility
from __future__ import absolute_import, division

# Imports
import os
import sys
import warnings
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime as dt
from matplotlib.projections import register_projection

# BioSPPy imports
import biosppy

# Local imports
import pyhrv
import pyhrv.time_domain
import pyhrv.frequency_domain
import pyhrv.nonlinear

# Turn off toolbox triggered warnings
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

	# if all(isinstance(n, int) for n in rpeaks) is False or all(isinstance(n, float) for n in rpeaks) is False:
	# 	raise TypeError("Incompatible data type in list or numpy array detected (only int or float allowed).")

	# Confirm numpy arrays & compute NN intervals
	rpeaks = np.asarray(rpeaks)
	nn_int = np.zeros(rpeaks.size - 1)

	for i in range(nn_int.size):
		nn_int[i] = rpeaks[i + 1] - rpeaks[i]

	return pyhrv.utils.nn_format(nn_int)


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
		nn = pyhrv.utils.nn_format(nni)

	# Confirm numpy arrays & compute NN interval differences
	nn_diff_ = np.zeros(nn.size - 1)

	for i in range(nn.size - 1):
		nn_diff_[i] = abs(nn[i + 1] - nn[i])

	return np.asarray(nn_diff_)


def plot_ecg(signal=None,
			 t=None,
			 sampling_rate=1000.,
			 interval=None,
			 rpeaks=True,
			 figsize=None,
			 title=None,
			 show=True):
	"""Plots ECG lead-I like signal on a medical grade ECG paper-like figure layout.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#plot-ecg-plot-ecg

	Parameters
	----------
	signal : array
		ECG lead-I like signal (filtered or unfiltered)
	t : array, optional
		Time vector for the ECG lead-I like signal (default: None)
	sampling_rate : int, float, optional
		Sampling rate of the acquired signal in [Hz] (default: 1000Hz)
	interval : array, 2-element, optional
		Visualization interval of the ECG lead-I like signal plot (default: None: [0s, 10s]
	rpeaks : bool, optional
		If True, marks R-peaks in ECG lead-I like signal (default: True)
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
		t = pyhrv.utils.time_vector(signal, sampling_rate=sampling_rate)

	# Configure interval of visualized signal
	if interval == 'complete':
		interval = [0, t[-1]]
	else:
		interval = pyhrv.utils.check_interval(interval, limits=[0, t[-1]], default=[0, 10])

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
	try:
		ax.set_xticks(np.arange(0.0, interval[1] + 0.1, float(n)/5), minor=True)
		ax.xaxis.grid(which='minor', color='salmon', lw=0.3)
		ax.set_xticks(np.arange(0, interval[1] + 0.1, n))
		ax.xaxis.grid(which='major', color='r', lw=0.7)
		ax.set_yticks(y_minor, minor=True)
		ax.yaxis.grid(which='minor', color='salmon', lw=0.3)
		ax.set_yticks(y_major)
		ax.yaxis.grid(which='major', color='r', lw=0.7)
	except:
		ax.grid(False)

	# Add legend
	unit = '' if unit == '-' else unit
	text_ = 'Division (x): %is\nDivision (y): %.1f%s' % (n, (np.abs(y_major[1] - y_major[0])), unit)
	ax.text(0.88, 0.85, text_, transform=ax.transAxes, fontsize=9,
		bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

	# Plot ECG lead-I like signal
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
		ax.set_title('ECG lead-I like signal - %s' % str(title))
	else:
		ax.set_title('ECG lead-I like signal')

	# Show plot
	if show:
		plt.show()

	# Output
	args = (fig_ecg, )
	names = ('ecg_plot', )
	return biosppy.utils.ReturnTuple(args, names)


def tachogram(nni=None,
			  signal=None,
			  rpeaks=None,
			  sampling_rate=1000.,
			  hr=True,
			  interval=None,
			  title=None,
			  figsize=None,
			  show=True):
	"""Plots Tachogram (NNI & HR) of an ECG lead-I like signal, NNI or R-peak series.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#tachogram-tachogram

	Parameters
	----------
	nni : array
		NN intervals in [ms] or [s].
	rpeaks : array
		R-peak times in [ms] or [s].
	signal : array, optional
		ECG lead-I like signal.
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
	..	NN intervals are derived from the ECG lead-I like signal if 'signal' is provided.
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
	nni = pyhrv.utils.check_input(nni, rpeaks)

	# Time vector back to ms
	t = np.cumsum(nni) / 1000.

	# Configure interval of visualized signal
	if interval == 'complete':
		interval = [0, t[-1]]
	else:
		interval = pyhrv.utils.check_interval(interval, limits=[0, t[-1]], default=[0, 10])

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

	try:
		n = int(interval[1] / 10)
		ax.set_xticks(np.arange(0, interval[1] + n, n))
	except Exception as e:
		ax.grid(False)

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
	args = (fig, )
	names = ('tachogram_plot', )
	return biosppy.utils.ReturnTuple(args, names)


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
			nni = pyhrv.utils.nn_format(nni) if len(nni) > 1 else nni[0]
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


def heart_rate_heatplot(nni=None,
						rpeaks=None,
						signal=None,
						sampling_rate=1000.,
						age=18,
						gender='male',
						interval=None,
						figsize=None,
						show=True):
	"""Graphical visualization & classification of HR performance based on normal HR ranges by age and gender.

	Docs: https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#heart-rate-heatplot-hr-heatplot

	Parameters
	----------
	nni : array
		NN intervals in [ms] or [s].
	rpeaks : array
		R-peak times in [ms] or [s].
	signal : array, optional
		ECG lead-I like signal.
	sampling_rate : int, float, optional
		Sampling rate of the acquired signal in [Hz].
	age : int, float
		Age of the subject (default: 18).
	gender : str
		Gender of the subject ('m', 'male', 'f', 'female'; default: 'male').
	interval : list, optional
		Sets visualization interval of the signal (default: [0, 10]).
	figsize : array, optional
		Matplotlib figure size (width, height) (default: (12, 4)).
	show : bool, optional
		If True, shows plot figure (default: True).

	Returns
	-------
	hr_heatplot : biosppy.utils.ReturnTuple object

	Raises
	------
	TypeError
		If no input data for 'nni', 'rpeaks' or 'signal' is provided

	Notes
	-----
	.. 	If both 'nni' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nni' data will be computed
		from the 'rpeaks'
	.. 	Modify the 'hr_heatplot.json' file to write own database values

	"""
	# Helper function
	def _get_classification(val, data):
		for key in data.keys():
			if data[key][0] <= int(val) <= data[key][1]:
				return key

	# Check input
	if signal is not None:
		rpeaks = biosppy.signals.ecg.ecg(signal=signal, sampling_rate=sampling_rate, show=False)[2]
	elif nni is None and rpeaks is None:
		raise TypeError('No input data provided. Please specify input data.')

	# Get NNI series
	nn = pyhrv.utils.check_input(nni, rpeaks)

	# Compute HR data and
	hr_data = heart_rate(nn)
	t = np.cumsum(nn) / 1000
	interval = pyhrv.utils.check_interval(interval, limits=[0, t[-1]], default=[0, t[-1]])

	# Prepare figure
	if figsize is None:
		figsize = (12, 5)
	fig, (ax, ax1, ax2) = plt.subplots(3, 1, figsize=figsize, gridspec_kw={'height_ratios': [12, 1, 1]})
	ax1.axis("off")
	fig.suptitle("Heart Rate Heat Plot (%s, %s)" % (gender, age))

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

	# Set gender
	if gender not in ["male", "m", "female", "f"]:
		raise ValueError("Unknown gender '%s' for this database." % gender)
	else:
		if gender == 'm':
			gender = 'male'
		elif gender == 'f':
			gender = 'female'

	# Load comparison data from database
	database = json.load(open(os.path.join(os.path.split(__file__)[0], './files/hr_heatplot.json')))

	# Get database values
	if age > 17:
		for key in database["ages"].keys():
			if database["ages"][key][0] - 1 < age < database["ages"][key][1] + 1:
				_age = database["ages"][key][0]

		color_map = database["colors"]
		data = database[gender][str(_age)]
		order = database["order"]

		# Plot with information based on reference database:
		# Create classifier counter (preparation for steps after the plot)
		classifier_counter = {}
		for key in data.keys():
			classifier_counter[key] = 0

		# Add threshold lines based on the comparison data
		for threshold in data.keys():
			ax.hlines(data[threshold][0], 0, t[-1], linewidth=0.4, alpha=1, color=color_map[threshold])
		ax.plot(t, hr_data, 'k--', linewidth=0.5)

		# Add colorized HR markers
		old_classifier = _get_classification(hr_data[0], data)
		start_index = 0
		end_index = 0
		for hr_val in hr_data:
			classifier_counter[old_classifier] += 1
			current_classifier = _get_classification(hr_val, data)
			if current_classifier != old_classifier:
				ax.plot(t[start_index:end_index], hr_data[start_index:end_index], 'o',
						markerfacecolor=color_map[old_classifier], markeredgecolor=color_map[old_classifier])
				start_index = end_index
				old_classifier = current_classifier
			end_index += 1

		# Compute distribution of HR values in %
		percentages = {}
		_left = 0
		legend = []
		ax2.tick_params(left=False)
		ax2.set_yticklabels([])
		for i in list(range(7)):
			classifier = str(order[str(i)][0])
			percentages[classifier] = float(classifier_counter[classifier]) / hr_data.size * 100
			ax2.barh(y=0, width=percentages[classifier], left=_left, color=color_map[classifier])
			_left += percentages[classifier]
			legend.append(mpl.patches.Patch(label="%s\n(%.2f%s)" % (order[str(i)][1], percentages[classifier], "$\%$"),
											fc=color_map[classifier]))
		ax.legend(handles=legend, loc=8, ncol=7)
	elif age <= 0:
		raise ValueError("Age cannot be <= 0.")
	else:
		warnings.warn("No reference data for age %i available." % age)
		ax.plot(t, hr_data, 'k--', linewidth=0.5)
		ax2.plot("", 0)

	# Set axis limits
	ax.axis([interval[0], interval[1], hr_data.min() * 0.7, hr_data.max() * 1.1])
	ax.set_ylabel('Heart Rate [$1/min$]')
	ax2.set_xlim([0, 100])
	ax2.set_xlabel("Distribution of HR over the HR classifiers [$\%$]")

	# Show plot
	if show:
		plt.show()

	# Output
	return biosppy.utils.ReturnTuple((fig, ), ('hr_heatplot', ))


def time_varying(nni=None, rpeaks=None, parameter='sdnn', window='n20', interpolation=None, show=True, mode='normal'):
	"""Computes time varying plot of a pyHRV parameter at every NNI of the input NNI (or rpeak) series using a moving
	time window or a moving NNI window.

	Parameters
	----------
	nni : array
		NN-Intervals in [ms] or [s]
	rpeaks : array
		R-peak locations in [ms] or [s]
	parameter : string
		pyHRV parameter key for which the time varying computation is to be plotted (check the hrv_keys.json file for a
		full list of available keys)
	window : string
		Time varying window configuration using the following syntax:
			'tX'	for using a moving time window, with X being the window interval before and after the current NNI
					Example:	t20 generates a time window of 20s before and 20s after each NNI for the computation
								of th pyHRV parameter
			OR
			'nX'	for using a moving NNI window, with X being the number of NNI included before and after the current
					NNI
					Example:	n20 generates a window which includes 20 NNI before and 20 NNI after the current NNI
	interpolation : int (optional)
		Frequency at which the computed parameter signal is be resampled and interpolated (for example to create a
		parameter signal with the same sampling frequency of the original ECG signal)
	show : bool, optional
		If true, show time varying plot (default: True)
	mode :

	Returns
	-------

	"""
	# Check input series
	nn = pyhrv.utils.check_input(nni, rpeaks)

	# Check if parameter is on the list of invalid parameters (computational time of these parameters are too long or
	# the parameters are input parameters for PSD functions
	invalid_parameters = ['plot', 'tinn_m', 'tinn_n', 'fft_nfft', 'fft_window', 'fft_resampling_frequency',
						  'fft_interpolation', 'ar_nfft', 'ar_order', 'lomb_nfft', 'lomb_ma']

	# Check selected parameter
	if parameter is None:
		raise TypeError("No parameter set for 'parameter'")
	elif parameter in invalid_parameters:
		raise ValueError("Parameter '%s' is not supported by this function. Please select another one." % parameter)
	elif parameter not in pyhrv.utils.load_hrv_keys_json().keys():
		raise ValueError("Unknown parameter '%s' (not a pyHRV parameter)." % parameter)

	# Check window and decode window configuration
	if window[0] != 't' and window[0] != 'n':
		raise ValueError("Invalid mode '%s'. Please select 't' for a time window or 'n' for a NNI window." % window[0])
	elif int(window[1:]) <= 0:
		raise ValueError("'window' cannot be <= 0.")
	else:
		window_mode = window[0]
		window_size = int(window[1:])

	# Internal helper function
	def _compute_parameter(array, func):
		try:
			# Try to pass the show and mode argument to to suppress PSD plots
			val = eval(func + '(nni=array, mode=\'dev\')[0][\'%s\']' % parameter)
		except TypeError as e:
			if 'mode' in str(e):
				try:
					# If functions has now mode feature but 'mode' argument, but a plotting feature
					val = eval(func + '(nni=array, plot=False)[\'%s\']' % parameter)
				except TypeError as a:
					try:
						val = eval(func + '(nni=array, show=False)[\'%s\']' % parameter)
					except TypeError as ae:
						if 'plot' in str(ae):
							# If functions has now plotting feature try regular function
							val = eval(func + '(nni=array)[\'%s\']' % parameter)
						else:
							val = eval(func + '(nni=array)[\'%s\']' % parameter)
		return val

	# Vars
	parameter_values = np.asarray([])

	# Get hrv_keys & the respective function
	hrv_keys = pyhrv.utils.load_hrv_keys_json()
	parameter_func = hrv_keys[parameter][-1]
	parameter_label = hrv_keys[parameter][1]
	parameter_unit = hrv_keys[parameter][2]

	# Beat window computation
	if window_mode == 'n':
		for i, _ in enumerate(nni):
			if i == 0:
				continue
			# Incomplete initial window
			elif i <= (window_size - 1):
				vals = nn[:(i + window_size + 1)]
				parameter_values = np.append(parameter_values, _compute_parameter(vals, parameter_func))
			# Complete Window
			elif i < (nni.size - window_size):
				vals = nn[i - window_size: i + window_size + 1]
				parameter_values = np.append(parameter_values, _compute_parameter(vals, parameter_func))
			# Incomplete ending window
			else:
				vals = nn[i - window_size:i]
				parameter_values = np.append(parameter_values, _compute_parameter(vals, parameter_func))

	# Time window computation
	elif window_mode == 't':
		t = np.cumsum(nn) / 1000
		for i, _t in enumerate(t):
			if i == 0:
				continue
			# Incomplete initial window
			elif _t <= window_size:
				# t_vals = np.where((t <= _t) & (t <== (_t + window_size)))
				indices = np.where(t <= (_t + window_size))[0]
				parameter_values = np.append(parameter_values, _compute_parameter(nn[indices], parameter_func))
			# Complete Window
			elif _t < t[-1] - window_size:
				indices = np.where(((_t - window_size) <= t) & (t <= (_t + window_size)))[0]
				parameter_values = np.append(parameter_values, _compute_parameter(nn[indices], parameter_func))
			# Incomplete end window
			else:
				indices = np.where(((_t - window_size) <= t) & (t <= t[-1]))[0]
				parameter_values = np.append(parameter_values, _compute_parameter(nn[indices], parameter_func))

	# Interpolation (optional) and time vector
	if interpolation is not None:
		t = np.cumsum(nn)
		f_interpol = sp.interpolate.interp1d(t, parameter_values, 'cubic')
		t = np.arange(t[0], t[-1], 1000. / interpolation)
		parameter_values = f_interpol(t)
		t /= 1000.
	else:
		t = np.cumsum(nn) / 1000

	# Define start and end intervals
	if window_mode == 'n':
		indices = np.arange(0, len(nn))
		start_interval = np.where(indices < window_size + 1)[0]
		valid_interval = np.where((indices >= (window_size + 1)) & (indices <= (indices[-1] - window_size)))[0]
		end_interval = np.where(indices > (indices[-1] - window_size))[0][:-1]
	elif window_mode == 't':
		start_interval = np.where(t < window_size)[0]
		valid_interval = np.where((t >= window_size) & (t <= t[-1] - window_size))[0]
		end_interval = np.where(t > t[-1] - window_size)[0][:-1]

	y_min, y_max = 0, parameter_values.max() * 1.2

	# Figure
	fig = plt.figure(figsize=(12, 4))
	ax = fig.add_subplot(111)
	_win_mode = "NNI Window: %i Intervals" % window_size if window_mode == 'n' else "Time Window: %is" % window_size
	fig.suptitle('Time Varying - %s Evolution' % parameter_label)
	ax.set_title('(%s)' % _win_mode, size=10)
	ax.set_ylabel('%s [$%s$]' % (parameter.upper(), parameter_unit))
	ax.set_xlim([0, t[-1]])
	ax.set_ylim([y_min, y_max])

	# Plot start values (except the very first NNI)
	ax.plot(t[1:window_size + 1], parameter_values[1:window_size + 1], 'r--')

	# Plot valid values
	ax.plot(t[valid_interval], parameter_values[valid_interval], 'g')

	# Plot final values
	ax.plot(t[end_interval], parameter_values[end_interval], 'r--')

	# X-Axis configuration
	# Set x-axis format to seconds if the duration of the signal <= 60s
	if t[-1] <= 60:
		ax.set_xlabel('Time [s]')
	# Set x-axis format to MM:SS if the duration of the signal > 60s and <= 1h
	elif 60 < t[-1] <= 3600:
		ax.set_xlabel('Time [MM:SS]')
		formatter = mpl.ticker.FuncFormatter(lambda ms, x: str(dt.timedelta(seconds=ms))[2:])
		ax.xaxis.set_major_formatter(formatter)
	# Set x-axis format to HH:MM:SS if the duration of the signal > 1h
	else:
		ax.set_xlabel('Time [HH:MM:SS]')
		formatter = mpl.ticker.FuncFormatter(lambda ms, x: str(dt.timedelta(seconds=ms)))
		ax.xaxis.set_major_formatter(formatter)

	# Window areas
	legends = []
	ax.vlines(t[window_size], y_min, y_max, color='r')
	ax.fill_between([0, t[window_size]], [y_max, y_max], facecolor='r', alpha=0.3)
	ax.vlines(t[parameter_values.size - window_size - 1], y_min, y_max, color='r')
	ax.fill_between([t[parameter_values.size - window_size - 1], t[-1]], [y_max, y_max], facecolor='r', alpha=0.3)
	legends.append(mpl.patches.Patch(fc='g', label='Complete Window'))
	legends.append(mpl.patches.Patch(fc='r', label='Incomplete Window', alpha=0.3))

	# Recommended minimum window size
	# TODO in future versions: add available recommended minimum durations to the HRV keys json file
	parameter_minimum = 50
	if t[window_size] < parameter_minimum:
		ax.vlines(parameter_minimum, y_min, y_max, color='orange')
		ax.fill_between([t[window_size], parameter_minimum], [y_max, y_max], color='orange', alpha=0.3)
		legends.append(mpl.patches.Patch(fc='orange', label='Recommended Minimum Window Size (%is)' %
															parameter_minimum, alpha=0.3))
	ax.legend(handles=legends, loc=8, framealpha=1., ncol=3)

	# Add overall value
	val = _compute_parameter(nn, parameter_func)
	ax.hlines(val, 0, t[-1], linestyles='--', linewidth=0.7)
	ax.text(1, val + 1, 'Overall')

	# Check mode
	if mode not in ['normal', 'dev', 'devplot']:
		warnings.warn("Unknown mode '%s'. Will proceed with 'normal' mode." % mode, stacklevel=2)
		mode = 'normal'

	if mode == 'normal':
		if show:
			plt.show()

		# Output
		args = (fig,)
		names = ("time_varying_%s" % parameter,)
		return biosppy.utils.ReturnTuple(args, names)

	elif mode == 'dev':
		return t, parameter_values, parameter

	elif mode == 'devplot':
		if mode == 'normal':
			if show:
				plt.show()

			# Output
			args = (fig, )
			names = ("time_varying_%s" % parameter, )
			return biosppy.utils.ReturnTuple(args, names), t, parameter_values, parameter


def radar_chart(nni=None,
				rpeaks=None,
				comparison_nni=None,
				comparison_rpeaks=None,
				parameters=None,
				reference_label='Reference',
				comparison_label='Comparison',
				show=True,
				legend=True):
	"""Plots a radar chart of HRV parameters to visualize the evolution the parameters computed from a NNI series
	(e.g. extracted from an ECG recording while doing sports) compared to a reference/baseline NNI series (
	e.g. extracted from an ECG recording while at rest).

	The radarchart normalizes the values of the reference NNI series with the values extracted from the baseline NNI
	series being used as the 100% reference values.

	Example: 	Reference NNI series: 	SDNN = 100ms → 100%
				Comparison NNI series: 	SDNN = 150ms → 150%

	The radar chart is not limited by the number of HRV parameters to be included in the chart; it dynamically
	adjusts itself to the number of compared parameters.

	Docs: https://pyhrv.readthedocs.io/en/latest/_pages/api/tools.html#radar-chart-radar-chart

	Parameters
	----------
	nni : array
		Baseline or reference NNI series in [ms] or [s] (default: None)
	rpeaks : array
		Baseline or referene R-peak series in [ms] or [s] (default: None)
	comparison_nni : array
		Comparison NNI series in [ms] or [s] (default: None)
	comparison_rpeaks : array
		Comparison R-peak series in [ms] or [s] (default: None)
	parameters : list
		List of pyHRV parameters (see keys of the hrv_keys.json file for a full list of available parameters).
		The list must contain more than 1 pyHRV parameters (default: None)
	reference_label : str, optional
		Plot label of the reference input data (e.g. 'ECG while at rest'; default: 'Reference')
	comparison_label : str, optional
		Plot label of the comparison input data (e.g. 'ECG while running'; default: 'Comparison')
	show : bool, optional
		If True, shows plot figure (default: True).
	legend : bool, optional
		If true, add a legend with the computed results to the plot (default: True)

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	reference_results : dict
		Results of the computed HRV parameters of the reference NNI series
		Keys: 	parameters listed in the input parameter 'parameters'
	comparison results : dict
		Results of the computed HRV parameters of the comparison NNI series
		Keys: 	parameters listed in the input parameter 'parameters'
	radar_plot :  matplotlib figure
		Figure of the generated radar plot

	Raises
	------
	TypeError
		If an error occurred during the computation of a parameter
	TypeError
		If no input data is provided for the baseline/reference NNI or R-peak series
	TypeError
		If no input data is provided for the comparison NNI or R-peak series
	TypeError
		If no selection of pyHRV parameters is provided
	ValueError
		If less than 2 pyHRV parameters were provided

	Notes
	-----
	.. 	If both 'nni' and 'rpeaks' are provided, 'rpeaks' will be chosen over the 'nn' and the 'nni' data will be computed
		from the 'rpeaks'
	.. 	If both 'comparison_nni' and 'comparison_rpeaks' are provided, 'comparison_rpeaks' will be chosen over the
		the 'comparison_nni' and the nni data will be computed from the 'comparison_rpeaks'

	"""
	# Helper function & variables
	para_func = pyhrv.utils.load_hrv_keys_json()
	unknown_parameters, ref_params, comp_params = [], {}, {}

	def _compute_parameter(nni_series, parameter):

		# Get function name for the requested parameter
		func = para_func[parameter][-1]

		try:
			# Try to pass the show and mode argument to to suppress PSD plots
			index = 0
			if parameter.endswith('_vlf'):
				parameter = parameter.replace('_vlf', '')
			elif parameter.endswith('_lf'):
				index = 1
				parameter = parameter.replace('_lf', '')
			elif parameter.endswith('_hf'):
				index = 2
				parameter = parameter.replace('_hf', '')
			val = eval(func + '(nni=nni_series, mode=\'dev\')[0][\'%s\']' % (parameter))
			val = val[index]
		except TypeError as e:
			if 'mode' in str(e):
				try:
					# If functions has now mode feature but 'mode' argument, but a plotting feature
					val = eval(func + '(nni=nni_series, plot=False)[\'%s\']' % parameter)
				except TypeError as a:
					if 'plot' in str(a):
						# If functions has now plotting feature try regular function
						val = eval(func + '(nni=nni_series)[\'%s\']' % parameter)
					else:
						raise TypeError(e)
		return val

	# Check input data
	if nni is None and rpeaks is None:
		raise TypeError("No input data provided for baseline or reference NNI. Please specify the reference NNI series.")
	else:
		nn = pyhrv.utils.check_input(nni, rpeaks)

	if comparison_nni is not None and comparison_rpeaks is not None:
		raise TypeError("No input data provided for comparison NNI. Please specify the comarison NNI series.")
	else:
		comp_nn = pyhrv.utils.check_input(comparison_nni, comparison_rpeaks)

	if parameters is None:
		raise TypeError("No input list of parameters provided for 'parameters'. Please specify a list of the parameters"
						"to be computed and compared.")
	elif len(parameters) < 2:
		raise ValueError("Not enough parameters selected for a radar chart. Please specify at least 2 HRV parameters "
						 "listed in the 'hrv_keys.json' file.")

	# Check for parameter that require a minimum duration to be computed & remove them if the criteria is not met
	if nn.sum() / 1000. <= 600 or comp_nn.sum() / 1000. <= 600:
		for p in ['sdann', 'sdnn_index']:
			if p in parameters:
				parameters.remove(p)
				warnings.warn("Input NNI series are too short for the computation of the '%s' parameter. This "
							  "parameter has been removed from the parameter list." % p, stacklevel=2)

	# Register projection of custom RadarAxes class
	register_projection(pyhrv.utils.pyHRVRadarAxes)

	# Check if the provided input parameter exists in pyHRV (hrv_keys.json) & compute available parameters
	for p in parameters:
		p = p.lower()
		if p not in para_func.keys():
			# Save unknown parameters
			unknown_parameters.append(p)
		else:
			# Compute available parameters
			ref_params[p] = _compute_parameter(nn, p)
			comp_params[p] = _compute_parameter(comp_nn, p)

			# Check if any parameters could not be computed (returned as None or Nan) and remove them
			# (avoids visualization artifacts)
			if np.isnan(ref_params[p]) or np.isnan(comp_params[p]):
				ref_params.pop(p)
				comp_params.pop(p)
				warnings.warn("The parameter '%s' could not be computed and has been removed from the parameter list."
							  % p)

	# Raise warning pointing out unknown parameters
	if unknown_parameters != []:
		warnings.warn("Unknown parameters '%s' will not be computed." % unknown_parameters, stacklevel=2)

	# Prepare plot
	colors = ['lightskyblue', 'salmon']
	if legend:
		fig, (ax_l, ax) = plt.subplots(1, 2, figsize=(12, 6), subplot_kw=dict(projection='radar'))
	else:
		fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={'projection': 'radar'})
	theta = np.linspace(0, 2 * np.pi, len(ref_params.keys()), endpoint=False)
	ax.theta = theta

	# Prepare plot data
	ax.set_varlabels([para_func[s][1].replace(' ', '\n') for s in ref_params.keys()])
	ref_vals = [100 for x in ref_params.keys()]
	com_vals = [comp_params[p] / ref_params[p] * 100 for p in ref_params.keys()]

	# Plot data
	for i, vals in enumerate([ref_vals, com_vals]):
		ax.plot(theta, vals, color=colors[i])
		ax.fill(theta, vals, color=colors[i], alpha=0.3)

	title = "HRV Parameter Radar Chart\nReference NNI Series (%s) vs. Comparison NNI Series (%s)\n" % (colors[0], colors[1]) \
			+ r"(Chart values in $\%$, Reference NNI parameters $\hat=$100$\%$)"

	# Add legend to second empty plot
	if legend:
		ax_l.set_title(title, horizontalalignment='center')
		legend = []

		# Helper function
		def _add_legend(label, fc="white"):
			return legend.append(mpl.patches.Patch(fc=fc, label="\n" + label))

		# Add list of computed parameters
		_add_legend(reference_label, colors[0])
		for p in ref_params.keys():
			_add_legend("%s:" % para_func[p][1])

		# Add list of comparison parameters
		_add_legend(comparison_label, colors[1])
		for p in ref_params.keys():
			u = para_func[p][2] if para_func[p][2] != "-" else ""
			_add_legend("%.2f%s vs. %.2f%s" % (ref_params[p], u, comp_params[p], u))

		# Add relative differences
		_add_legend("")
		for i, _ in enumerate(ref_params.keys()):
			val = com_vals[i] - 100
			_add_legend("+%.2f%s" % (val, r"$\%$") if val > 0 else "%.2f%s" % (val, r"$\%$"))

		ax_l.legend(handles=legend, ncol=3, frameon=False, loc=7)
		ax_l.axis('off')
	else:
		ax.set_title(title, horizontalalignment='center')

	# Show plot
	if show:
		plt.show()

	# Output
	args = (ref_params, comp_params, fig, )
	names = ('reference_results', 'comparison_results', 'radar_plot', )
	return biosppy.utils.ReturnTuple(args, names)


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
		Unsupported data format provided (other than dict, or biosppy.utils.ReturnTuple object.)
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
	# Check input (if available & biosppy.utils.ReturnTuple object)
	if results is None:
		raise TypeError("No results data provided. Please specify input data.")
	elif results is not type(dict()) and isinstance(results, biosppy.utils.ReturnTuple) is False:
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

	efile, _ = pyhrv.utils.check_fname(path, 'json', efile)

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
		if isinstance(results[key], biosppy.utils.ReturnTuple):
			output[key] = dict(results[key])
		elif isinstance(results[key], tuple):
			output[key] = list(results[key])
		elif isinstance(results[key], str):
			output[key] = results[key]
		elif isinstance(results[key], range):
			output[key] = list(results[key])
		elif results[key] is None:
			output[key] = 'n/a'
		elif 'plot' not in str(key) and 'histogram' not in str(key):
			output[key] = float(results[key]) if str(results[key]) != 'nan' else 'n/a'

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
		results[str(key)] = data[key] if type(data[key]) is not str else str(data[key])

	# Create biosppy.utils.ReturnTuple object from imported data
	return biosppy.utils.ReturnTuple(results.values(), results.keys())


if __name__ == "__main__":
	"""
	Example Script - HRV Tools
	"""
	import pyhrv
	from biosppy.signals.ecg import ecg

	# Load a Sample Signal
	nni = pyhrv.utils.load_sample_nni()
	heart_rate_heatplot(nni)

	# # Load OpenSignals (r)evolution ECG sample file
	# signal = np.loadtxt('./files/SampleECG.txt')[:, -1]
	#
	# # Filter data & get r-peak locations [ms]
	# signal, rpeaks = ecg(signal, show=False)[1:3]
	#
	# # Plot ECG for the interval of 0s and 22s
	# plot_ecg(signal, interval=[0, 22])
	#
	# # Plot Tachogram for the interval of 0s and 22s
	# tachogram(nni, interval=[0, 22])
	#
	# # Heart Rate Heatplot to highlight HR performance compared to a sports database
	# heart_rate_heatplot(nni, gender='male', age=28)
	#
	# # Time Varying is designed to show the evolution of HRV parameters over time using a moving window
	# # Define a moving window of 3 NNIs before and after the current NNI using the NNI window indicator 'n'
	# time_varying(nni, parameter='sdnn', window='n3')
	#
	# # Define a moving window of 3 seconds before and after the current NNI using the time window indicator 't'
	# time_varying(nni, parameter='sdnn', window='t3')
	#
	# # Radar charts are created dynamically, depending on the number of parameters used as input
	# # For this example, let's split he test NNI series into two segments & select a list of 6 parameters
	# ref_nni = nni[:100]
	# comp_nni = nni[100:200]
	# params = ['nni_mean', 'nni_max', 'sdnn', 'rmssd', 'sdsd', 'nn50', 'nn20']
	# radar_chart(ref_nni, comparison_nni=comp_nni, parameters=params)
	#
	# # Now with only 3 parameters
	# params = ['nni_mean', 'sdnn', 'rmssd']
	# radar_chart(ref_nni, comparison_nni=comp_nni, parameters=params)
	#
	# # Export and import HRV results into and from JSON files:
	# # First, compute hrv parameters
	# results = pyhrv.hrv(nni, show=False)
	#
	# hrv_export(results, path='./files/', efile='SampleExport')
	# hrv_import('./files/SampleExport.json')
