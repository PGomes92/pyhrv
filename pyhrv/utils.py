# -*- coding: utf-8 -*-
"""
pyHRV - Utilities (utils)
-------------------------

This module contains general purpose functions (utilities) to support the features of the pyHRV toolbox incl. loading
NNI sample data, check input data of HRV functions, segment arrays, check for duplicate file names to avoid accidental
overwriting of files, and other utilities.

Notes
-----
..  Up to v.0.3 this work has been developed within the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..	This module has been added in v.0.4
..	You can find the API reference for this module here:
	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html

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
from __future__ import absolute_import, division

# Imports
import os
import warnings
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.projections.polar import PolarAxes

# BioSPPy import
import biosppy
import pyhrv 


def load_sample_nni(series='short'):
	"""Returns a short-term (5min) or long-term (60min) series of sample NNI found in the pyhrv/files/ directory.

	Docs:	
	
	Parameters
	----------
	series : string, optional
		If 'long', returns a 60min NNI series, if 'short', returns a 5min NNI series

	Returns
	-------
	nni_series : array
		Sample NNI series

	Raises
	------
	ValueError
		If an unknown value for the 'series' input parameter is provided (other than 'short' or 'long')

	Note
	----
	.. 	These sample series were extracted from the MIT-BIH NSRDB Database from physionet.org:
		https://physionet.org/physiobank/database/nsrdb/

	"""
	if series == 'long':
		return np.load(os.path.join(os.path.split(__file__)[0], './files/SampleNNISeriesLong.npy'))
	elif series == 'short':
		return np.load(os.path.join(os.path.split(__file__)[0], './files/SampleNNISeriesShort.npy'))
	else:
		raise ValueError("Unknown input value '%s'. Please select 'short' or 'long'." % series)


def load_hrv_keys_json():
	"""Loads the content of the 'hrv_keys.json' file found in the 'pyhrv/files/' directory.

	Returns
	-------
	hrv_keys : dict
		Content of the pyhrv/files/hrv_keys.json file in a dictionary
	"""
	return json.load(open(os.path.join(os.path.split(__file__)[0], './files/hrv_keys.json'), 'r'))


def check_input(nni=None, rpeaks=None):
	"""Checks if input series of NN intervals or R-peaks are provided and, if yes, returns a NN interval series in [ms]
	format.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#check-input-check-input

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

	"""
	# Check input
	if nni is None and rpeaks is not None:
		# Compute NN intervals if r_peaks array is given
		nni = pyhrv.tools.nn_intervals(rpeaks)
	elif nni is not None:
		# Use given NN intervals & confirm numpy
		nni = nn_format(nni)
	else:
		raise TypeError("No R-peak data or NN intervals provided. Please specify input data.")
	return nni


def nn_format(nni=None):
	"""Checks format of the NN intervals (seconds or milliseconds) and converts s data to ms data, if necessary.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#nn-format-nn-format

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
		https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#nn-format-nn-format

	"""
	# Check input
	if nni is None:
		raise TypeError("No input data provided for 'nn'. Please specify input data")
	nn_ = np.asarray(nni, dtype='float64')

	# Convert if data has been identified in [s], else proceed with ensuring the NumPy array format
	if np.max(nn_) < 10:
		nn_ = [int(x * 1000) for x in nn_]

	return np.asarray(nn_)


def check_interval(interval=None, limits=None, default=None):
	"""General purpose function that checks and verifies correctness of interval limits within optionally defined
	valid interval specifications and and/or default values if no interval is specified.

	This function can be used to set visualization intervals, check overlapping frequency bands, or for other similar
	purposes and is intended to automatically catch possible error sources due to invalid intervals.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#check-interval-check-interval

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
			warnings.warn("Interval limits out of boundaries. Interval s€et to: %s" % interval, stacklevel=2)
		if not limits[1] >= interval[1]:
			interval[1] = limits[1]
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
		warnings.warn("The provided lower limit of the parameter '%s' is greater than the upper limit. Limits have been switched." % name)
	if interval[0] == interval[1]:
		raise ValueError("The provided lower and upper limits of the parameter '%s' are invalid as they are identical." % name)
	return interval


def segmentation(nni=None, full=True, duration=300, warn=True):
	"""Segmentation of signal peaks into individual segments of set duration.
	(e.g. splitting R-peak locations into 5min segments for computation of SDNN index)

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#segmentation-segmentation

	Parameters
	----------
	nni : array
		Series of NN intervals [ms] or [s].
	full : bool, optional
		If True, returns last segment, even if the sum of NNI does not reach the segment duration (default: False).
	duration : int, optional
		Segment duration in [s] (default: 300s)
	warn : bool, optional
		If True, raise a warning message if a segmentation could not be conducted (duration > NNI series duration)

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

	"""
	# Check input
	if nni is None:
		raise TypeError("Please specify input data.")

	# Preparations
	nni = nn_format(nni)
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
		if warn:
			warnings.warn("Signal duration is to short for segmentation into %is. Input data will be returned." % duration)
		return nni, False


def join_tuples(*args):
	"""Joins multiple biosppy.biosppy.utils.ReturnTuple objects into one biosppy.biosppy.utils.ReturnTuple object.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#join-tuples-join-tuples

	Parameters
	----------
	tuples : list, array, biosppy.utils.ReturnTuple objects
		List or array containing biosppy.utils.ReturnTuple objects.

	Returns
	-------
	return_tuple : biosppy.biosppy.utils.ReturnTuple object
		biosppy.biosppy.utils.ReturnTuple object with the content of all input tuples joined together.

	Raises
	------
	TypeError:
		If no input data is provided
	TypeError:
		If input data contains non-biosppy.biosppy.utils.ReturnTuple objects

	Notes
	----
	..	You can find the documentation for this function here:
		https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#join-tuples-join-tuples

	"""
	# Check input
	if args is None:
		raise TypeError("Please specify input data.")

	for i in args:
		if not isinstance(i, biosppy.utils.ReturnTuple):
			raise TypeError("The list of tuples contains non-biosppy.utils.ReturnTuple objects.")

	# Join tuples
	names = ()
	vals = ()

	for i in args:
		for key in i.keys():
			names = names + (key, )
			vals = vals + (i[key], )

	return biosppy.utils.ReturnTuple(vals, names)


def check_fname(rfile, fformat='txt', name='new_file'):
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

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#standard-deviation-std

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

	"""
	if array is None:
		raise TypeError("Please specify input data.")

	array = np.asarray(array)
	result = np.sum([(x - np.mean(array))**2 for x in array])
	result = np.sqrt(1. / (array.size - dof) * result)
	return result


def time_vector(signal=None, sampling_rate=1000.):
	"""Computes time vector for the provided input signal.

	Docs:	https://pyhrv.readthedocs.io/en/latest/_pages/api/utils.html#time-vector-time-vector

	Parameters
	----------
	signal : array
		ECG lead-I like signal (or any other sensor signal)
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

	"""
	if signal is None:
		raise TypeError("Please specify input signal.")

	signal = np.asarray(signal)
	t = np.arange(0, signal.size / sampling_rate, 1./sampling_rate)
	return t


class pyHRVRadarAxes(PolarAxes):
	"""Custom child class inheriting from the PolarAxes class for the radar_chart() function"""
	name = 'radar'

	def __init__(self, *args, **kwargs):
		super(pyHRVRadarAxes, self).__init__(*args, **kwargs)
		self.set_theta_zero_location('N')
		self.theta = []

	def plot(self, *args, **kwargs):
		"""Override plot so that line is closed by default"""
		lines = super(pyHRVRadarAxes, self).plot(*args, **kwargs)
		for line in lines:
			self._close_line(line)

	def _close_line(self, line):
		x, y = line.get_data()
		if x[0] != x[-1]:
			x = np.concatenate((x, [x[0]]))
			y = np.concatenate((y, [y[0]]))
			line.set_data(x, y)

	def set_varlabels(self, labels):
		self.set_thetagrids(np.degrees(self.theta), labels)

	def _gen_axes_patch(self):
		"""Overwrite method"""
		return plt.Circle((0.5, 0.5), 0.5)

	def _gen_axes_spines(self):
		"""Overwrite method"""
		return PolarAxes._gen_axes_spines(self)
