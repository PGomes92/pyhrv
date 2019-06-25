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


def waterfall_psd_plot(nni=None, rpeaks=None, fbands=None, show=True, method='ar', duration=300, colored=True,
					   ar_order=32):
	# Check input values
	nn = tools.check_input(nni, rpeaks)

	# Verify or set default frequency bands
	fbands = _check_freq_bands(fbands)

	# Vars
	z_axis = []
	verts = []
	powers = []
	freqs = []

	# Segment signal
	segments, worked = tools.segmentation(nni, duration=duration, full=False)
	n_segments = len(segments)

	# Create 3D figure object
	fig = plt.figure()
	wf_ax = fig.gca(projection='3d')

	if worked:
		# Create multiple Welch's plot from each segment
		for i in range(n_segments):
			val = segments[i]

			# Plot psd
			if method == 'ar':
				_, f, p = ar_psd(val, mode='dev', order=ar_order)
			elif method == 'welch':
				_, f, p = welch_psd(val, mode='dev')
			elif method == 'lomb':
				pass
				#f, p = lomb_psd(val, mode='dev')
			else:
				raise ValueError("Unknown method '%s' selected. Please select a valid PSD estimation method (welch, "
								 "ar or lomb)" % method)

			# Store data
			freqs.append(np.asarray([float(fval) for fval in f if fval <= 0.4]))
			powers.append(np.asarray(p[:len(f)]))

		# Normalize data and zip it for the plot
		for i, power in enumerate(powers):
			# # Normalize data
			p = np.asarray([p_val/(power.max() - power.min()) for p_val in power])
			f = np.asarray(freqs[i])

			# Set values at 0Hz and 0.4Hz (upper limit of HF band) to 0 to prevent plotting visualization errors
			p[0] = 0
			p[-1] = 0

			# Add plot data
			z_axis.append((i + 1) * float(duration)/60)
			verts.append(list(zip(f, p)))

		# Get Colors
		if colored:
			facecolors = [plt.cm.jet(x * (1. / n_segments)) for x in range(n_segments)]
		else:
			facecolors = ['k']
			for i in segments[:-1]:
				color = 'k' if facecolors[-1] != 'k' else 'darkgrey'
				facecolors.append(color)
		# Plot
		poly = PolyCollection(verts, facecolors=facecolors, linewidths=0.7)
		poly.set_alpha(0.7)
		wf_ax.add_collection3d(poly, zs=z_axis, zdir='y')

		# Plot configurations
		wf_ax.set_xlim3d(0, 0.4)
		wf_ax.set_xlabel("Frequency")
		wf_ax.set_ylim3d(0, z_axis[-1])
		wf_ax.set_ylabel("t [min]")
		wf_ax.set_zlim3d(0, 1)
		wf_ax.set_zlabel("$ms^2 (Normalized)$")

		# Add fbands
		ys = [0, z_axis[-1]]
		_axis = ['']
		vlf = list(zip([fbands['vlf'][1], fbands['vlf'][1]], ys))
		lf = list(zip([fbands['lf'][1], fbands['lf'][1]], ys))
		hf = list(zip([fbands['hf'][1], fbands['hf'][1]], ys))
		lines = []
		lines.append(LineCollection((hf, hf), color='k', lw=0.5, linestyle='--'))
		lines.append(LineCollection((lf, lf), color='k', lw=0.5, linestyle='--'))
		lines.append(LineCollection((vlf, lf), color='k', lw=0.5, linestyle='--'))
		for line in lines:
			wf_ax.add_collection3d(line, zs=0, zdir='z')

		zs = [0, 1]
		vlf = list(zip([fbands['vlf'][1], fbands['vlf'][1]], zs))
		lf = list(zip([fbands['lf'][1], fbands['lf'][1]], zs))
		hf = list(zip([fbands['hf'][1], fbands['hf'][1]], zs))
		lines.append(LineCollection((hf, hf), color='k', lw=0.5, linestyle='--'))
		lines.append(LineCollection((lf, lf), color='k', lw=0.5, linestyle='--'))
		lines.append(LineCollection((vlf, lf), color='k', lw=0.5, linestyle='--'))
		for line in lines:
			wf_ax.add_collection3d(line, zs=55, zdir='y')

		# Add fbands
		axes = [[0, z_axis[-1]], [0, 1]]
		for _ax in axes:
			vlf = list(zip([fbands['vlf'][1], fbands['vlf'][1]], _ax))
			lf = list(zip([fbands['lf'][1], fbands['lf'][1]], zs))
			hf = list(zip([fbands['hf'][1], fbands['hf'][1]], zs))
			lines.append(LineCollection((hf, hf), color='k', lw=0.5, linestyle='--'))
			lines.append(LineCollection((lf, lf), color='k', lw=0.5, linestyle='--'))
			lines.append(LineCollection((vlf, lf), color='k', lw=0.5, linestyle='--'))

		# Show plot
		if show:
			plt.show()

	else:
		raise ValueError("The duration of the signal (%is) is too short to be segmented into %is segments."
			  % (np.sum(nn) / 1000, duration))

nni=np.load('SampleNNISeries.npy')
waterfall_psd_plot(nni, duration=60, method='welch')
