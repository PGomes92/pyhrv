.. _ref-parallelmodule:

Parallel Processing: parallel
==============================

Module: ``pyhrv.parallel``

The ``parallel`` module provides multiprocessing utilities for concurrent HRV computation. It uses Python's
``concurrent.futures.ProcessPoolExecutor`` to run independent analysis domains in separate processes.

No additional dependencies are required beyond the Python standard library.

.. _ref-parallel-hrv:

parallel_hrv()
--------------

.. py:function:: pyhrv.parallel.parallel_hrv(nni, kwargs_time=None, fbands=None, kwargs_welch=None, kwargs_lomb=None, kwargs_ar=None, kwargs_nonlinear=None, max_workers=None)

**Function Description**

Computes time domain, frequency domain, and nonlinear HRV parameters concurrently using 3 worker processes
(one per domain). This is the function called internally by ``hrv(parallel=True)``.

**Input Parameters**
   - ``nni`` (array): NN intervals in [ms] or [s] (already validated)
   - ``kwargs_time`` (dict, optional): Kwargs for ``time_domain()``
   - ``fbands`` (dict, optional): Frequency band specifications for ``frequency_domain()``
   - ``kwargs_welch`` (dict, optional): Kwargs for ``welch_psd()``
   - ``kwargs_lomb`` (dict, optional): Kwargs for ``lomb_psd()``
   - ``kwargs_ar`` (dict, optional): Kwargs for ``ar_psd()``
   - ``kwargs_nonlinear`` (dict, optional): Kwargs for ``nonlinear()``. Unpacked into ``kwargs_poincare``, ``kwargs_sampen``, ``kwargs_dfa``
   - ``max_workers`` (int, optional): Maximum number of worker processes (default: min(3, cpu_count))

**Returns**
   - ``t_results`` (ReturnTuple): Time domain results
   - ``f_results`` (ReturnTuple): Frequency domain results
   - ``n_results`` (ReturnTuple): Nonlinear results

**Example**

.. code-block:: python

   from pyhrv.parallel import parallel_hrv
   import pyhrv.utils

   nni = pyhrv.utils.load_sample_nni()
   t_results, f_results, n_results = parallel_hrv(nni=nni)

   print(t_results['sdnn'])
   print(f_results['fft_peak'])
   print(n_results['sd1'])

.. _ref-parallel-time:

parallel_time_domain()
----------------------

.. py:function:: pyhrv.parallel.parallel_time_domain(nni, threshold=None, plot=False, binsize=7.8125)

**Function Description**

Computes individual time domain functions (``nni_parameters``, ``sdnn``, ``rmssd``, etc.) concurrently using
up to 5 worker processes.

**Input Parameters**
   - ``nni`` (array): NN intervals in [ms]
   - ``threshold`` (int, optional): Custom threshold for ``nnXX`` parameter
   - ``plot`` (bool, optional): If True, creates histogram plot (default: False)
   - ``binsize`` (float, optional): Histogram bin size in [ms] (default: 7.8125)

**Returns**
   - ``results`` (ReturnTuple): All time domain results

**Example**

.. code-block:: python

   from pyhrv.parallel import parallel_time_domain
   import pyhrv.utils

   nni = pyhrv.utils.load_sample_nni()
   results = parallel_time_domain(nni=nni)
   print(results['sdnn'])
   print(results['rmssd'])

.. _ref-parallel-freq:

parallel_frequency_domain()
---------------------------

.. py:function:: pyhrv.parallel.parallel_frequency_domain(nni, fbands=None, show=False, kwargs_welch=None, kwargs_lomb=None, kwargs_ar=None)

**Function Description**

Computes Welch, Lomb-Scargle, and Autoregressive PSD methods concurrently using 3 worker processes.

**Input Parameters**
   - ``nni`` (array): NN intervals in [ms]
   - ``fbands`` (dict, optional): Frequency band specifications
   - ``show`` (bool, optional): If True, show plots after computation (default: False)
   - ``kwargs_welch`` (dict, optional): Kwargs for ``welch_psd()``
   - ``kwargs_lomb`` (dict, optional): Kwargs for ``lomb_psd()``
   - ``kwargs_ar`` (dict, optional): Kwargs for ``ar_psd()``

**Returns**
   - ``results`` (ReturnTuple): Combined frequency domain results from all three methods

**Example**

.. code-block:: python

   from pyhrv.parallel import parallel_frequency_domain
   import pyhrv.utils

   nni = pyhrv.utils.load_sample_nni()
   results = parallel_frequency_domain(nni=nni)
   print(results['fft_peak'])
   print(results['lomb_peak'])
   print(results['ar_peak'])

.. _ref-parallel-notes:

Notes
-----

.. note::

   - Parallel mode uses the ``Agg`` matplotlib backend in child processes. Plot figure objects (e.g. ``fft_plot``,
     ``poincare_plot``) are returned as ``None`` when using parallel functions.
   - ``biosppy.utils.ReturnTuple`` objects do not support pickling natively. The parallel module handles this by
     converting results to plain dicts in worker processes and reconstructing ``ReturnTuple`` objects in the parent.
   - For interactive plotting, use the sequential (non-parallel) functions instead.
