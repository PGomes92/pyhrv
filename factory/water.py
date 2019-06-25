# Existing imports
from __future__ import division
import numpy as np
from pyhrv import tools
from pyhrv.frequency_domain import welch_psd, ar_psd, lomb_psd, _check_freq_bands

# New imports
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection, LineCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors


def waterfall_psd_plot(nni=None,
					   rpeaks=None,
					   segments=None,
					   method='ar',
					   fbands=None,
					   kwargs_method={},
					   segment_duration=300,
					   show=True,
					   legend=True):
	"""Creates 3D waterfall plot of PSD plots.

	Can be used to create a plot for comparison of multiple of PSD plots.

	Parameters
	----------
	nni : array
		NN-Intervals in [ms] or [s]
	rpeaks : array
		R-peak locations in [ms] or [s]
	segments : array of arrays
		Array with series of NNI from which the PSD plots should be computed.
	fbands : dict, optional
		Dictionary with frequency bands (2-element tuples or list)
		Value format:	(lower_freq_band_boundary, upper_freq_band_boundary)
		Keys:	'ulf'	Ultra low frequency		(default: none) optional
				'vlf'	Very low frequency		(default: (0.000Hz, 0.04Hz))
				'lf'	Low frequency			(default: (0.04Hz - 0.15Hz))
				'hf'	High frequency			(default: (0.15Hz - 0.4Hz))
	method : string, optional
		PSD computation method (default: 'ar'):
			'ar':		Autoregressive 				ar_psd()
			'welch':	Welch Method 				welch_psd()
			'lomb':		Lomb-Scargle periodogram 	lomb_psd()
	kwargs_method : dictionary, optional
		Dictionary with input parameters for the 'welch_psd()', 'ar_psd()' or 'lomb_psd()' method (default: {}).
	segment_duration : int, optional
		Duration of NNI segments from which the PSD are computed (default: 30)
	show : bool, optional
		If true, show PSD plot (default: True)
	colored : bool, optional
		If true, each segment in the waterfall PSD plot will get a different color, else it will be alternating black
		and gray (default: True)

	Returns (biosppy.utils.ReturnTuple Object)
	------------------------------------------
	[key : format]
		Description.
	X_waterfall_plot : matplotlib figure
		3D waterfall plot figure.
	X_waterfall_parameters : array of ReturnTuple objects
		Array of ReturnTuple objects containing the frequency domain parameters for each NNI segment
		(1st ReturnTuple object = parameters of the first NNI segment, 2nd ReturnTuple object = parameters of the
		second NNI segment, ...)

	Raises
	------
	ValueError


	Notes
	-----
	..	Only one type of input data is required (nni, rpeaks or segments)
	.. 	If both 'nni' and 'rpeaks' are provided, 'nni' will be chosen over the 'rpeaks'
	..  Segments will be chosen over 'nni' or 'rpeaks'
	..	NN and R-peak series provided in [s] format will be converted to [ms] format.

	"""
	# Check input values
	if segments is None:
		nn = tools.check_input(nni, rpeaks)
		y_unit_time = True
		segments, worked = tools.segmentation(nn, duration=segment_duration, full=False)
	else:
		y_unit_time = False
		worked = True

	# If segmentation is not possible, raise a warning and return the regular PSD plot
	if not worked:
		warnings.warn("NNI series is too short. Segmentation cannot be conducted using the defined segment duration. "
					  "Will return regular %s plot." % method)

		# Return Autorgressive plot
		if method == 'ar':
			# Unwrwap kwargs dictionary
			nfft = kwargs_method['nfft'] if 'nfft' in kwargs_method.keys() else 2 ** 8
			order = kwargs_method['order'] if 'order' in kwargs_method.keys() else 16

			nn = nn if nn is not None else segments
			return ar_psd(nni=nn, show=show, nfft=nfft, order=order)

		# Return Welch PSD plot
		elif method == 'welch':
			# Unwrap kwargs dictionary
			nfft = kwargs_method['nfft'] if 'nfft' in kwargs_method.keys() else 2 ** 12
			detrend = kwargs_method['detrend'] if 'detrend' in kwargs_method.keys() else True
			window = kwargs_method['window'] if 'window' in kwargs_method.keys() else 'hamming'

			# Compute PSD of current segment
			return welch_psd(nni=nn, nfft=nfft, detrend=detrend, window=window, show=show)

		# Return Lomb PSD plot
		elif method == 'lomb':
			# Unwrap kwargs dictionary
			nfft = kwargs_method['nfft'] if 'nfft' in kwargs_method.keys() else 2 ** 8
			ma_size = kwargs_method['ma_size'] if 'ma_size' in kwargs_method.keys() else None

			return lomb_psd(nni=nn, nfft=nfft, ma_size=ma_size, show=False)

	# If segmentation worked, proceed with plot
	else:
		# Verify or set default frequency bands
		fbands = _check_freq_bands(fbands)

		# Vars
		z_axis, plot_data, powers, freqs, segment_parameters = [], [], [], [], []

		# Create multiple Welch's plot from each segment
		for i, segment in enumerate(segments):
			# Get to PSD methods
			if method == 'ar':
				method_name = 'Autoregressive'

				# Unwrwap kwargs dictionary
				nfft = kwargs_method['nfft'] if 'nfft' in kwargs_method.keys() else 2 ** 12
				order = kwargs_method['order'] if 'order' in kwargs_method.keys() else 16
				a = segment
				# Compute PSD of current segment
				params, f, p = ar_psd(nni=segment, mode='dev', order=order, nfft=nfft, show=False, legend=False)

			elif method == 'welch':
				method_name = 'Welch\'s Method'

				# Unwrap kwargs dictionary
				nfft = kwargs_method['nfft'] if 'nfft' in kwargs_method.keys() else 2**12
				detrend = kwargs_method['detrend'] if 'detrend' in kwargs_method.keys() else True
				window = kwargs_method['window'] if 'window' in kwargs_method.keys() else 'hamming'

				# Compute PSD of current segment
				params, f, p = welch_psd(nni=segment, mode='dev', nfft=nfft, detrend=detrend, window=window, show=False)

			elif method == 'lomb':
				method_name = 'Lomb Scargle'

				# Unwrap kwargs dictionary
				nfft = kwargs_method['nfft'] if 'nfft' in kwargs_method.keys() else 2**8
				ma_size = kwargs_method['ma_size'] if 'ma_size' in kwargs_method.keys() else None

				params, f, p = lomb_psd(nni=segment, mode='dev', nfft=nfft, ma_size=ma_size, show=False)
				p[0] = 0
			else:
				raise ValueError("Unknown method '%s' selected. Please select a valid PSD estimation method (welch, "
								 "ar or lomb)." % method)

			# Store parameter results
			segment_parameters.append(params)

			# Get intervals between 0 and upper interval of the HF band
			f = [float(fval) for fval in f if fval <= fbands['hf'][1]]
			p = p[:len(f)]

			# Normalize data
			p = np.asarray([(p_val - np.min(p))/(np.max(p) - np.min(p)) for p_val in p])

			# Store values
			freqs.append(f)
			powers.append(p)

			# Zip data
			if y_unit_time:
				z_axis.append((i + 1) * float(segment_duration)/60)
			else:
				z_axis.append(i+1)

		# Create 3D figure object
		fig = plt.figure(figsize=(10, 4))
		wf_ax = Axes3D(fig)
		alpha = np.linspace(0.5, 1, len(freqs))
		colors = {'ulf': 'b', 'vlf': 'yellowgreen', 'lf': 'salmon', 'hf': 'lightskyblue'}

		# Plot segments
		for i, f_array in enumerate(freqs):
			f_array = np.asarray(f_array)

			# Get band specific frequencies
			_, vlf_i, lf_i, hf_i = _get_frequency_indices(f_array, fbands)

			# NOTE: Plottig frequency bands over each other reduces the amount of visualization artifacts
			fband_dict = {
				'vlf': np.append(vlf_i, np.append(lf_i, hf_i)),
				'lf': np.append(lf_i, hf_i),
				'hf': hf_i
			}

			for band in ['vlf', 'lf', 'hf']:
				# NOTE: do not try plot each plot each frequency band through loops, as this will cause visualization
				# artifacts!

				# Plot VLF band
				f = np.asarray(f_array[fband_dict[band]])
				p = np.asarray(powers[i][fband_dict[band]])
				p[0], p[-1] = 0, 0
				data = np.array([[f, p]]).transpose()
				poly = LineCollection(segments=data, facecolors=colors[band], linewidths=0.7, alpha=alpha[i])
				poly = LineCollection(data)
				wf_ax.add_collection3d(poly, zs=z_axis[i], zdir='y')

		# Axis settings & title
		if method == 'ar':
			wf_ax.set_title('3D Waterfall Plot - Autoregressive (Order: %i)' % order)
		else:
			wf_ax.set_title('3D Waterfall Plot - %s' % method_name)
		wf_ax.set_xlim([0, fbands['hf'][-1]])
		wf_ax.set_xlabel("Frequency [$Hz$]")
		wf_ax.set_ylim([0, z_axis[-1]])
		y_label = "Time [$min$]" if y_unit_time else "Segment [-]"
		wf_ax.set_ylabel(y_label)
		wf_ax.set_zlim([0, 1])
		wf_ax.set_zlabel("Power [$ms^2/Hz$]")
		wf_ax.view_init(elev=25, azim=-55)

		# Add legend
		if legend:
			legend = []
			legend.append(mpl.patches.Patch(facecolor=colors['vlf'], label='VLF'))
			legend.append(mpl.patches.Patch(facecolor=colors['lf'], label='LF'))
			legend.append(mpl.patches.Patch(facecolor=colors['hf'], label='HF'))
			wf_ax.legend(handles=legend, loc=6)

		# Show
		if show:
			plt.show()

		# Output
		return fig, segments

nni = np.load('SampleNNISeries.npy')
waterfall_psd_plot(nni, segment_duration=60, method='welch')
