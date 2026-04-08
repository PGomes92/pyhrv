# -*- coding: utf-8 -*-
"""
pyHRV - Parallel Processing Module
-----------------------------------

Provides multiprocessing utilities for parallel HRV computation.
The three HRV analysis domains (time, frequency, nonlinear) are independent
and can be computed concurrently using process pools.

Note: biosppy.utils.ReturnTuple does not pickle correctly, so worker
functions serialize results to plain dicts and the parent process
reconstructs ReturnTuple objects.

:copyright: (c) 2019 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.
"""

from __future__ import absolute_import

import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import biosppy.utils
import pyhrv.utils

# Force non-interactive backend in child processes
_PARALLEL_BACKEND = 'Agg'


def _rt_to_dict(rt):
    """Recursively convert a ReturnTuple to a pickle-safe dict.

    Handles nested ReturnTuples and drops matplotlib Figure objects
    (they can't be meaningfully transferred across processes).
    """
    result = {}
    for key in rt.keys():
        val = rt[key]
        if isinstance(val, biosppy.utils.ReturnTuple):
            result[key] = ('__rt__', _rt_to_dict(val))
        elif hasattr(val, 'savefig'):
            # matplotlib Figure — skip it, not needed for numeric results
            result[key] = ('__figure__', None)
        elif hasattr(val, 'tolist'):
            result[key] = val.tolist()
        else:
            result[key] = val
    return result


def _dict_to_rt(d):
    """Reconstruct a ReturnTuple from a dict produced by _rt_to_dict."""
    keys = []
    vals = []
    for k, v in d.items():
        keys.append(k)
        if isinstance(v, tuple) and len(v) == 2 and v[0] == '__rt__':
            vals.append(_dict_to_rt(v[1]))
        elif isinstance(v, tuple) and len(v) == 2 and v[0] == '__figure__':
            vals.append(None)
        else:
            vals.append(v)
    return biosppy.utils.ReturnTuple(tuple(vals), tuple(keys))


# ---------------------------------------------------------------------------
# Worker functions (run in child processes, return plain dicts)
# ---------------------------------------------------------------------------

def _worker_time_domain(nn_list, kwargs_time):
    """Compute time domain in a subprocess, return serializable dict."""
    import matplotlib
    matplotlib.use(_PARALLEL_BACKEND)
    import numpy as np
    import pyhrv.time_domain as td

    nn = np.array(nn_list)
    if kwargs_time is not None:
        threshold = kwargs_time.get('threshold', None)
        binsize = kwargs_time.get('binsize', 7.8125)
        plot = kwargs_time.get('plot', False)
        result = td.time_domain(nni=nn, show=False, threshold=threshold, plot=plot, binsize=binsize)
    else:
        result = td.time_domain(nni=nn, show=False, plot=False)
    return _rt_to_dict(result)


def _worker_frequency_domain(nn_list, fbands, kwargs_welch, kwargs_lomb, kwargs_ar):
    """Compute frequency domain in a subprocess, return serializable dict."""
    import matplotlib
    matplotlib.use(_PARALLEL_BACKEND)
    import numpy as np
    import pyhrv.frequency_domain as fd

    nn = np.array(nn_list)
    result = fd.frequency_domain(
        nni=nn, fbands=fbands,
        kwargs_welch=kwargs_welch, kwargs_lomb=kwargs_lomb, kwargs_ar=kwargs_ar,
        show=False
    )
    return _rt_to_dict(result)


def _worker_nonlinear(nn_list, kwargs_poincare, kwargs_sampen, kwargs_dfa):
    """Compute nonlinear parameters in a subprocess, return serializable dict."""
    import matplotlib
    matplotlib.use(_PARALLEL_BACKEND)
    import numpy as np
    import pyhrv.nonlinear as nl

    nn = np.array(nn_list)
    result = nl.nonlinear(
        nni=nn, show=False,
        kwargs_poincare=kwargs_poincare,
        kwargs_sampen=kwargs_sampen,
        kwargs_dfa=kwargs_dfa
    )
    return _rt_to_dict(result)


def _worker_welch(nn_list, fbands, kwargs_welch):
    """Compute Welch PSD in a subprocess."""
    import matplotlib
    matplotlib.use(_PARALLEL_BACKEND)
    import numpy as np
    import pyhrv.frequency_domain as fd

    nn = np.array(nn_list)
    if kwargs_welch is not None:
        detrend = kwargs_welch.get('detrend', True)
        window = kwargs_welch.get('window', 'hamming')
        nfft = kwargs_welch.get('nfft', 2**12)
        result = fd.welch_psd(nn, fbands=fbands, detrend=detrend, show=False, nfft=nfft, window=window)
    else:
        result = fd.welch_psd(nn, show=False, fbands=fbands)
    return _rt_to_dict(result)


def _worker_lomb(nn_list, fbands, kwargs_lomb):
    """Compute Lomb-Scargle PSD in a subprocess."""
    import matplotlib
    matplotlib.use(_PARALLEL_BACKEND)
    import numpy as np
    import pyhrv.frequency_domain as fd

    nn = np.array(nn_list)
    if kwargs_lomb is not None:
        nfft = kwargs_lomb.get('nfft', 2**8)
        ma_size = kwargs_lomb.get('ma_size', None)
        result = fd.lomb_psd(nn, fbands=fbands, ma_size=ma_size, show=False, nfft=nfft)
    else:
        result = fd.lomb_psd(nn, show=False, fbands=fbands)
    return _rt_to_dict(result)


def _worker_ar(nn_list, fbands, kwargs_ar):
    """Compute AR PSD in a subprocess."""
    import matplotlib
    matplotlib.use(_PARALLEL_BACKEND)
    import numpy as np
    import pyhrv.frequency_domain as fd

    nn = np.array(nn_list)
    if kwargs_ar is not None:
        nfft = kwargs_ar.get('nfft', 2**12)
        order = kwargs_ar.get('order', 16)
        result = fd.ar_psd(nn, fbands=fbands, order=order, show=False, nfft=nfft)
    else:
        result = fd.ar_psd(nn, show=False, fbands=fbands)
    return _rt_to_dict(result)


# ---------------------------------------------------------------------------
# Public parallel API
# ---------------------------------------------------------------------------

def parallel_hrv(nni, kwargs_time=None, fbands=None,
                 kwargs_welch=None, kwargs_lomb=None, kwargs_ar=None,
                 kwargs_nonlinear=None, max_workers=None):
    """Compute time, frequency, and nonlinear HRV domains in parallel.

    Parameters
    ----------
    nni : array_like
        NN intervals in [ms] or [s] (already validated).
    kwargs_time : dict, optional
        Kwargs for time_domain().
    fbands : dict, optional
        Frequency bands for frequency_domain().
    kwargs_welch, kwargs_lomb, kwargs_ar : dict, optional
        Method-specific kwargs for frequency domain.
    kwargs_nonlinear : dict, optional
        Kwargs for nonlinear(). Will be unpacked into kwargs_poincare,
        kwargs_sampen, kwargs_dfa.
    max_workers : int, optional
        Max number of worker processes (default: 3, one per domain).

    Returns
    -------
    t_results, f_results, n_results : tuple of ReturnTuple objects
    """
    if max_workers is None:
        max_workers = min(3, multiprocessing.cpu_count())

    # Unpack nonlinear kwargs
    kwargs_poincare = None
    kwargs_sampen = None
    kwargs_dfa = None

    if kwargs_nonlinear is not None:
        if not isinstance(kwargs_nonlinear, dict):
            raise TypeError("'kwargs_nonlinear' must be a dictionary.")
        kwargs_poincare = {
            'ellipse': kwargs_nonlinear.get('ellipse', True),
            'vectors': kwargs_nonlinear.get('vectors', True),
            'legend': kwargs_nonlinear.get('legend', True),
            'marker': kwargs_nonlinear.get('marker', 'o'),
        }
        kwargs_sampen = {
            'dim': kwargs_nonlinear.get('dim', 2),
            'tolerance': kwargs_nonlinear.get('tolerance', None),
        }
        kwargs_dfa = {
            'short': kwargs_nonlinear.get('short', None),
            'long': kwargs_nonlinear.get('long', None),
            'legend': kwargs_nonlinear.get('legend', True),
        }

    nn_list = list(nni)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        ft = executor.submit(_worker_time_domain, nn_list, kwargs_time)
        ff = executor.submit(_worker_frequency_domain, nn_list, fbands,
                             kwargs_welch, kwargs_lomb, kwargs_ar)
        fn = executor.submit(_worker_nonlinear, nn_list,
                             kwargs_poincare, kwargs_sampen, kwargs_dfa)

        t_results = _dict_to_rt(ft.result())
        f_results = _dict_to_rt(ff.result())
        n_results = _dict_to_rt(fn.result())

    return t_results, f_results, n_results


def parallel_time_domain(nni, threshold=None, plot=False, binsize=7.8125):
    """Compute time domain sub-functions in parallel.

    Each statistical function (sdnn, rmssd, etc.) is independent and
    submitted to a process pool concurrently.

    Parameters
    ----------
    nni : array_like
        NN intervals in [ms].
    threshold : int, optional
        Custom threshold for nnXX.
    plot : bool, optional
        If True, creates histogram plot.
    binsize : float, optional
        Histogram bin size.

    Returns
    -------
    results : biosppy.utils.ReturnTuple
    """
    import numpy as np
    import pyhrv.time_domain as td

    nn = np.array(nni)
    nn_list = list(nn)

    # These are fast functions — use ThreadPoolExecutor-like approach
    # but since they release the GIL via numpy, ProcessPool is fine
    with ProcessPoolExecutor(max_workers=min(5, multiprocessing.cpu_count())) as executor:
        futures = {
            'nni': executor.submit(_worker_td_func, nn_list, 'nni_parameters'),
            'hr': executor.submit(_worker_td_func, nn_list, 'hr_parameters'),
            'diff': executor.submit(_worker_td_func, nn_list, 'nni_differences_parameters'),
            'sdnn': executor.submit(_worker_td_func, nn_list, 'sdnn'),
            'sdnn_index': executor.submit(_worker_td_func, nn_list, 'sdnn_index'),
            'sdann': executor.submit(_worker_td_func, nn_list, 'sdann'),
            'rmssd': executor.submit(_worker_td_func, nn_list, 'rmssd'),
            'sdsd': executor.submit(_worker_td_func, nn_list, 'sdsd'),
            'nn50': executor.submit(_worker_td_func, nn_list, 'nn50'),
            'nn20': executor.submit(_worker_td_func, nn_list, 'nn20'),
        }

        # Collect results in order
        order = ['nni', 'hr', 'diff', 'sdnn', 'sdnn_index', 'sdann',
                 'rmssd', 'sdsd', 'nn50', 'nn20']
        results = _dict_to_rt(futures[order[0]].result())
        for key in order[1:]:
            results = pyhrv.utils.join_tuples(results, _dict_to_rt(futures[key].result()))

    if threshold is not None and threshold not in [50, 20]:
        results = pyhrv.utils.join_tuples(results, td.nnXX(nn, threshold=int(threshold)))

    results = pyhrv.utils.join_tuples(
        results, td.geometrical_parameters(nn, plot=plot, show=False, binsize=binsize)
    )
    return results


def parallel_frequency_domain(nni, fbands=None, show=False,
                              kwargs_welch=None, kwargs_lomb=None, kwargs_ar=None):
    """Compute Welch, Lomb-Scargle, and AR PSDs in parallel.

    Parameters
    ----------
    nni : array_like
        NN intervals in [ms].
    fbands : dict, optional
        Frequency bands.
    show : bool, optional
        If True, show plots after computation.
    kwargs_welch, kwargs_lomb, kwargs_ar : dict, optional
        Method-specific kwargs.

    Returns
    -------
    results : biosppy.utils.ReturnTuple
    """
    import numpy as np
    nn = np.array(nni)
    nn_list = list(nn)

    with ProcessPoolExecutor(max_workers=3) as executor:
        fw = executor.submit(_worker_welch, nn_list, fbands, kwargs_welch)
        fl = executor.submit(_worker_lomb, nn_list, fbands, kwargs_lomb)
        fa = executor.submit(_worker_ar, nn_list, fbands, kwargs_ar)

        welch_results = _dict_to_rt(fw.result())
        lomb_results = _dict_to_rt(fl.result())
        ar_results = _dict_to_rt(fa.result())

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return pyhrv.utils.join_tuples(welch_results, lomb_results, ar_results)


def _worker_td_func(nn_list, func_name):
    """Generic worker for individual time domain functions."""
    import matplotlib
    matplotlib.use(_PARALLEL_BACKEND)
    import numpy as np
    import pyhrv.time_domain as td

    nn = np.array(nn_list)
    func = getattr(td, func_name)
    result = func(nn)
    return _rt_to_dict(result)
