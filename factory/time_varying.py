import pyhrv
import warnings
import scipy as sp
import numpy as np
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyhrv import tools
import biosppy
import pyhrv.utils as utils


# TODO add dev_plot
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
	nn = tools.check_input(nni, rpeaks)

	# Check if parameter is on the list of invalid parameters (computational time of these parameters are too long or
	# the parameters are input parameters for PSD functions
	invalid_parameters = ['plot', 'tinn_m', 'tinn_n', 'fft_nfft', 'fft_window', 'fft_resampling_frequency',
						  'fft_interpolation', 'ar_nfft', 'ar_order', 'lomb_nfft', 'lomb_ma']

	# Check selected parameter
	if parameter is None:
		raise TypeError("No parameter set for 'parameter'")
	elif parameter in invalid_parameters:
		raise ValueError("Parameter '%s' is not supported by this function. Please select another one." % parameter)
	elif parameter not in utils.get_pyhrv_keys().keys():
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
					if 'plot' in str(a):
						# If functions has now plotting feature try regular function
						val = eval(func + '(nni=array)[\'%s\']' % parameter)
					else:
						raise TypeError(e)
		return val

	# Vars
	parameter_values = np.asarray([])

	# Get hrv_keys & the respective function
	# TODO change -1 to -2 when using the updated hrv_keys json
	hrv_keys = utils.get_pyhrv_keys()
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
	# TODO add available recommended minimum durations to the HRV keys json file
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
			args = (fig,)
			names = ("time_varying_%s" % parameter,)
			return biosppy.utils.ReturnTuple(args, names), t, parameter_values, parameter


if __name__ == '__main__':
	#
	nni = np.load('SampleNNISeries.npy')
	time_varying(nni=nni, parameter='sdnn', window='n30')


