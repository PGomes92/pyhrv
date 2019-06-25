from __future__ import division

import numpy as np
import warnings
import json
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection

from pyhrv import tools
import pyhrv
from biosppy import utils


def radar_chart(nni=None,
                rpeaks=None,
                comparison_nni=None,
                comparison_rpeaks=None,
                parameters=None,
                reference_label='Reference',
                comparison_label='Comparison',
                show=True,
                legend=True):
    """

    Parameters
    ----------
    nni
    rpeaks
    comparison_nni
    comparison_rpeaks
    parameters
    reference_label
    comparison_label
    show
    legend

    Returns
    -------

    """
    # Check input data
    if nni is None and rpeaks is None:
        raise TypeError("No input data provided for baseline or reference NNI. Please specify the reference NNI series.")
    else:
        nn = tools.check_input(nni, rpeaks)

    if comparison_nni is not None and comparison_rpeaks is not None:
        raise TypeError("No input data provided for comparison NNI. Please specify the comarison NNI series.")
    else:
        comp_nn = tools.check_input(comparison_nni, comparison_rpeaks)

    if parameters is None:
        raise TypeError("No input list of parameters provided for 'reference'. Please specify a list of the parameters"
                        "to be computed and compared.")
    elif len(parameters) < 2:
        raise ValueError("Not enough parameters selected for a radar chart. Please specify at least 2 HRV parameters "
                         "listed in the 'hrv_keys.json' file.")

    # Register projection of custom RadarAxes class
    register_projection(RadarAxes)

    # Get HRV keys &
    # TODO change path to something that works in general
    para_func = json.load(open(os.path.join(os.path.split(__file__)[0], '../pyhrv/files/hrv_keys.json'), 'r'))
    unknown_parameters, ref_params, comp_params = [], {}, {}

    # Check if the provided input parameter exists in pyHRV (hrv_keys.json) & compute available parameters
    for p in parameters:
        p = p.lower()
        if p not in para_func.keys():
            # Save unknwon parameters
            unknown_parameters.append(p)
        else:
            # Compute available parameters
            ref_params[p] = eval(para_func[p][-1] + "(nni=nn)[p]")
            comp_params[p] = eval(para_func[p][-1] + "(nni=comp_nn)[p]")

    # Raise warning pointing out the unknwon parameters
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
    ax.set_varlabels([s.upper() for s in ref_params.keys()])
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
            _add_legend("%s:" % p.upper())

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
    args = (ref_params, comp_params, )
    names = ('ref_nni_segment', 'comparison_nni_segment', )
    return utils.ReturnTuple(args, names)


class RadarAxes(PolarAxes):
    """Child class of the PolarAxes class for the HRV parameter radar chart"""
    name = 'radar'

    def __init__(self, *args, **kwargs):
        super(RadarAxes, self).__init__(*args, **kwargs)
        self.set_theta_zero_location('N')
        self.theta = []

    def plot(self, *args, **kwargs):
        """Override plot so that line is closed by default"""
        lines = super(RadarAxes, self).plot(*args, **kwargs)
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

if __name__ == '__main__':
    ref_nni = np.load('../TestData/series_1.npy')
    comp_nni = np.load('../TestData/series_6.npy')

    params = ['nni_mean', 'nni_mean', 'sdnn', 'rmssd', 'sdsd', 'sdnn_index', 'sdann', 'nn50', 'nn20']
    radar_chart(ref_nni, comparison_nni=comp_nni, parameters=params)
    radar_chart(ref_nni, comparison_nni=comp_nni, parameters=params[:4])
    radar_chart(ref_nni, comparison_nni=comp_nni, parameters=params[-5:])

