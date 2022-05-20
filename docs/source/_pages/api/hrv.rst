.. _ref-hrvfunc:

The HRV Function: hrv()
=======================

.. py:function:: pyhrv.hrv.hrv(nni=None, rpeaks=None, signal=None, sampling_rate=1000., interval=[0, 10], plot_ecg=True, plot_Tachogram=True, show=True, fbands=None, kwargs_ecg_plot=None, kwargs_tachogram=None, kwargs_)

**Function Description**

Computes all HRV parameters of the pyHRV toolkit (see list below).

.. seealso::

   * :ref:`ref-timemodule`
   * :ref:`ref-frequencymodule`
   * :ref:`ref-nonlinearmodule`

**Input Parameters**
   - ``nni`` (array): NN intervals in [ms] or [s]
   - ``rpeaks`` (array): R-peak times in [ms] or [s]
   - ``signal`` (array): ECG signal
   - ``sampling_rate`` (int, float, optional): Sampling rate in [Hz] used for the ECG acuqisition (default: 1000Hz)
   - ``interval`` (array, optional): Visualization interval of the Tachogram plot (default: None: [0s, 10s])
   - ``plot_ecg`` (bool, optional): If True, plots ECG signal with specified interval ('signal' must not be None)
   - ``plot_tachogram`` (bool, optional): If True, plots tachogram with specified interval
   - ``show`` (bool, optional): If True, shows the ECG plot figure (default: True)
   - ``fbands`` (dict, optional): Dictionary with frequency band specifications (default: None)
   - ``kwargs_ecg_plot`` (dict, optional): **kwargs for the plot_ecg() function (see 'tools.py' module)
   - ``kwargs_tachogram`` (dict, optional): **kwargs for the plot_tachogram() function (see 'tools.py' module)
   - ``kwargs_time`` (dict, optional): **kwargs for the time_domain() function (see 'time_domain()' function)
   - ``kwargs_welch`` (dict, optional): **kwargs for the 'welch_psd()' function (see 'frequency_domain.py' module)
   - ``kwargs_lomb`` (dict, optional): **kwargs for the 'lomb_psd()' function (see 'frequency_domain.py' module)
   - ``kwargs_ar`` (dict, optional): **kwargs for the 'ar_psd()' function (see 'frequency_domain.py' module)
   - ``kwargs_nonlinear`` (dict, optional): **kwargs for the 'nonlinear() function (see 'nonlinear.py' module)

.. note::

   If ``fbands`` is none, the default values for the frequency bands will be set.

      * VLF:   [0.00Hz - 0.04Hz]
      * LF:    [0.04Hz - 0.15Hz]
      * HF:    [0.15Hz - 0.40Hz]

   See **Application Notes** & **Examples & Tutorials** below for more information on how to define custom frequency bands.

   The ``show`` parameter is equally set for all plotting functions.

.. important::

   This function computes the Time Domain parameters using either the ``signal``, ``nni``, or ``rpeaks`` data. Provide
   only one type of data, as it is not required to pass all three types at once.

**Returns (ReturnTuple Object)**

The results of this function are returned in a ``biosppy.utils.ReturnTuple`` object. Use the following keys below (on the left) to index the results:

Time Domain:

   - ``nni_counter`` (int): Number of NNI (-)
   - ``nni_mean`` (float): Mean NNI [ms]
   - ``nni_min`` (int): Minimum NNI [ms]
   - ``nni_max`` (int): Maximum NNI [ms]
   - ``nni_diff_mean`` (float): Mean NNI difference [ms]
   - ``nni_diff_min`` (int): Minimum NNI difference [ms]
   - ``nni_diff_max`` (int): Maximum NNI difference [ms]
   - ``hr_mean`` (float): Mean heart rate [bpm]
   - ``hr_min`` (int): Minimum heart rate [bpm]
   - ``hr_max`` (int): Maximum heart rate [bpm]
   - ``hr_std`` (float): Standard deviation of the heart rate series [bpm]
   - ``sdnn`` (float): Standard deviation of NN intervals [ms]
   - ``sdnn_index`` (float): SDNN Index [ms]
   - ``sdann`` (float): SDANN [ms]
   - ``rmssd`` (float): Root mean of squared NNI differences [ms]
   - ``sdsd`` (float): Standard deviation of NNI differences [ms]
   - ``nnXX`` (int, optional): Number of NN interval differences greater than the specified threshold (-)
   - ``pnnXX`` (float, optional): Ratio between nnXX and total number of NN interval differences (-)
   - ``nn50`` (int): Number of NN interval differences greater 50ms
   - ``pnn50`` (float): Ratio between NN50 and total number of NN intervals [ms]
   - ``nn20`` (int): Number of NN interval differences greater 20ms
   - ``pnn20`` (float): Ratio between NN20 and total number of NN intervals [ms]
   - ``nn_histogram`` (matplotlib figure object): Histogram plot figure (only if input parameter ``plot`` is True
   - ``tinn_n`` (float): N value of the TINN computation (left corner of the interpolated triangle at (N, 0))
   - ``tinn_m`` (float): M value of the TINN computation (right corner of the interpolated triangle at (M, 0))
   - ``tinn`` (float): TINN (baseline width of the interpolated triangle) [ms]
   - ``tri_index`` (float): Triangular index [ms]

.. important::

   The ``XX`` in the ``nnXX`` and the ``pnnXX`` keys are substituted by the specified threshold.

   For instance, ``nnXX(nni, threshold=30)`` returns the custom ``nn30`` and ``pnn30`` parameters. Applying
   ``threshold=35`` as ``nnXX(nni, threshold=35)`` returns the custom ``nn35`` and ``pnn35`` parameters.

   These parameters are only returned if a custom threshold (``threshold``) has been defined in the input parameters.

Frequency Domain (X = one of the methods 'fft', 'ar', 'lomb'):

   - ``X_peak`` (tuple): Peak frequencies of all frequency bands [Hz]
   - ``X_abs`` (tuple): Absolute powers of all frequency bands [ms^2]
   - ``X_rel`` (tuple): Relative powers of all frequency bands [%]
   - ``X_log`` (tuple): Logarithmic powers of all frequency bands [log]
   - ``X_norm`` (tuple): Normalized powers of the LF and HF frequency bands [-]
   - ``X_ratio`` (float): LF/HF ratio [-]
   - ``X_total`` (float): Total power over all frequency bands [ms^2]
   - ``X_plot`` (matplotlib figure object): PSD plot figure object
   - ``fft_interpolation`` (str): Interpolation method used for NNI interpolation (hard-coded to 'cubic')
   - ``fft_resampling_frequency`` (int): Resampling frequency used for NNI interpolation [Hz] (hard-coded to 4Hz as recommended by the `HRV Guidelines <https://www.ahajournals.org/doi/full/10.1161/01.cir.93.5.1043>`_)
   - ``fft_window`` (str): Spectral window used for PSD estimation of the Welch's method
   - ``lomb_ma`` (int): Moving average window size
   - ``ar_interpolation`` (str): Interpolation method used for NNI interpolation (hard-coded to 'cubic')
   - ``ar_resampling_frequency`` (int): Resampling frequency used for NNI interpolation [Hz] (hard-coded to 4Hz as recommended by the `HRV Guidelines <https://www.ahajournals.org/doi/full/10.1161/01.cir.93.5.1043>`_)
   - ``ar_order`` (int): Autoregressive model order

Nonlinear:

   - ``poincare_plot`` (matploltib figure object): Poincaré plot figure
   - ``sd1`` (float): Standard deviation (SD1) of the major axis
   - ``sd2`` (float): Standard deviation (SD2) of the minor axis
   - ``sd_ratio`` (float): Ratio between SD1 and SD2 (SD2/SD1)
   - ``ellipse_area`` (float): Arrea S of the fitted ellipse
   - ``sample_entropy`` (float): Sample entropy of the NNI series
   - ``dfa_short`` (float): Alpha value of the short-term fluctuations (alpha1)
   - ``dfa_long`` (float): Alpha value of the long-term fluctuations (alpha2)

.. seealso::

   :ref:`ref-returntuple`

**Application Notes**

It is not necessary to provide input data for ``signal``, ``nni`` **and** ``rpeaks``. The parameter(s) of this
function will be computed with any of the input data provided (``signal``, ``nni`` **or** ``rpeaks``). The input data will be prioritized in the following order, in case multiple inputs are provided:

1. ``signal``, 2. ``nni``, 3. ``rpeaks``.

``nni`` or ``rpeaks`` data provided in seconds [s] will automatically be converted to ``nni`` data in  milliseconds [ms].

.. seealso::

   Section :ref:`ref-nnformat` for more information.

.. important::

   This function generates ``matplotlib`` plot figures which, depending on the backend you are using, can interrupt
   your code from being executed whenever plot figures are shown. Switching the backend and turning on the
   ``matplotlib`` interactive mode can solve this behavior.

   In case it does not - or if switching the backend is not possible - close all the plot figures to proceed with the
   execution of the rest your code after the ``plt.show()`` function.

   .. seealso::

      * :ref:`ref-matplotlib-workaround`
      * `More information about the matplotlib Interactive Mode <https://matplotlib.org/faq/usage_faq.html#what-is-interactive-mode>`_
      * `More information about matplotlib Backends <https://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`_

Incorrect frequency band specifications will be automatically corrected, if possible. For instance the following frequency bands contain overlapping frequency band limits which would cause issues when computing the frequency parameters:

.. code-block:: python

   fbands = {'vlf': (0.0, 0.25), 'lf': (0.2, 0.3), 'hf': (0.3, 0.4)}

Here, the upper band of the VLF band is greater than the lower band of the LF band. In this case, the overlapping frequency band limits will be switched:

.. code-block:: python

   fbands = {'vlf': (0.0, 0.2), 'lf': (0.25, 0.3), 'hf': (0.3, 0.4)}

.. warning::

   Corrections of frequency bands trigger ``warnings`` which are displayed in the Python console. It is recommended to watch out for these warnings and to correct the frequency bands given that the corrected bands might not be optimal.

Use the ``kwargs_ecg_plot`` dictionary to pass function specific parameters for the ``plot_ecg()`` function. The following keys are supported:

   - ``rpeaks`` (bool, optional): If True, marks R-peaks in ECG signal (default: True)
   - ``title`` (str, optional): Plot figure title (default: None)

Use the ``kwargs_tachogram`` dictionary to pass function specific parameters for the ``plot_tachogram()`` function. The following keys are supported:

   - ``hr`` (bool, optional): If True, plot HR seres in [bpm] on second axis (default: True)
   - ``title`` (str, optional): Optional plot figure title (default: None)

Use the ``kwargs_time`` dictionary to pass function specific parameters for the ``time_domain()`` function. The following keys are supported:

   - ``threshold`` (int, optional): Custom threshold in [ms] for the optional NNXX and pNNXX parameters (default: None)
   - ``plot`` (bool, optional): If True, creates histogram using matploltib, elss uses NumPy for histogram data only (geometrical parameters, default: True)
   - ``binsize`` (float, optional): Bin size in [ms] of the histogram bins - (geometrical params, default: 7.8125ms).

Use the ``kwargs_welch`` dictionary to pass function specific parameters for the ``welch_psd()`` method. The following keys are supported:

   - ``nfft`` (int, optional): Number of points computed for the FFT result (default: 2**12)
   - ``detrend`` (bool, optional): If True, detrend NNI series by subtracting the mean NNI (default: True)
   - ``window`` (scipy.window function, optional): Window function used for PSD estimation (default: 'hamming')

Use the ``lomb_psd`` dictionary to pass function specific parameters for the ``lombg_psd()`` method. The following keys are supported:

   - ``nfft`` (int, optional): Number of points computed for the Lomb-Scargle result (default: 2**8)
   - ``ma_order`` (int, optional): Order of the moving average filter (default: None; no filter applied)

Use the ``ar_psd`` dictionary to pass function specific parameters for the ``ar_psd()`` method. The following keys are supported:

   - ``nfft`` (int, optional): Number of points computed for the FFT result (default: 2**12)
   - ``order`` (int, optional): Autoregressive model order (default: 16)

Use the ``kwargs_nonlinear`` dictionary to pass function specific parameters for the ``nonlinear()`` function. The
following keys are supported:

   - ``ellipse`` (bool, optional): If True, shows fitted ellipse in plot (default: True)
   - ``vectors`` (bool, optional): If True, shows SD1 and SD2 vectors in plot (default: True)
   - ``legend`` (bool, optional): If True, adds legend to the Poincaré plot (default: True)
   - ``marker`` (str, optional): NNI marker in plot (must be compatible with the matplotlib markers (default: 'o')
   - ``dim`` (int, optional): Entropy embedding dimension (default: 2)
   - ``tolerance`` (int, float, optional): Tolerance distance for which the two vectors can be considered equal (default: std(NNI))
   - ``short`` (array, optional): Interval limits of the short-term fluctuations (default: None: [4, 16])
   - ``long`` (array, optional): Interval limits of the long-term fluctuations (default: None: [17, 64])

**Examples & Tutorials & Tutorials**

The following example codes demonstrate how to use the ``hrv()`` function.

You can choose either the ECG signal, the NNI series or the R-peaks as input data for the PSD estimation and
parameter computation:

.. code-block:: python

   # Import packages
   import biosppy
   import pyhrv.tools as tools
   from pyhrv.hrv import hrv

   # Load sample ECG signal
   signal = np.loadtxt('./files/SampleECG.txt')[:, -1]

   # Get R-peaks series using biosppy
   t, filtered_signal, rpeaks = biosppy.signals.ecg.ecg(signal)[:3]

   # Compute NNI series
   nni = tools.nn_intervals(t[rpeaks])

   # OPTION 1: Compute Time Domain parameters using the ECG signal
   signal_results = hrv(signal=filtered_signal)

   # OPTION 2: Compute Time Domain parameters using the R-peak series
   rpeaks_results = hrv(rpeaks=t[rpeaks])

   # OPTION 3: Compute Time Domain parameters using the NNI-series
   nni_results = hrv(nni=nni)

The output of of all three options above will be the same.

.. note::

   If an ECG signal is provided, the signal will be filtered and the R-peaks will be extracted using the
   ``biosppy.signals.ecg.ecg()`` function. Finally, the NNI series for the PSD estimation will be computed from the extracted
   R-peak series. The ECG plot is only generated if an ECG signal is provided.

.. seealso::

   `biosppy.signals.ecg.ecg() <https://biosppy.readthedocs.io/en/stable/biosppy.signals.html#biosppy.signals.ecg
   .ecg>`_

You can now access the parameters using the output parameter keys (works the same for the ``rpeaks_results`` and
``nni_results``):

.. code-block:: python

   # Print SDNN
   print(signal_results['sdnn'])

   # Print RMSSD
   print(signal_results['rmssd'])

Use the `kwargs` input dictionaries to provide custom input parameters.

.. code-block:: python

   # Define custom input parameters using the kwargs dictionaries
   kwargs_time = {'threshold': 35}
   kwargs_nonlinear = {'vectors': False}
   kwargs_welch = {'nfft': 2**8}
   kwargs_lomb = {'nfft': 2**16}
   kwargs_ar = {'nfft': 2**8}
   kwargs_tachogram = {'hr': False}
   kwargs_ecg_plot = {'title': 'My ECG Signal'}

   # Compute HRV parameters
   hrv(nni=nni, kwargs_time=kwargs_time, kwargs_nonlinear=kwargs_nonlinear, kwargs_ar=kwargs_ar,
      kwargs_lomb=kwargs_lomb, kwargs_welch=kwargs_welch, kwargs_tachogram=kwargs_tachogram)

pyHRV is robust against invalid parameter keys. For example, if an invalid input parameter such as `nfft` is
provided with the `kwargs_time` dictionary, this parameter will be ignored and a warning message will
be issued.

.. code-block:: python

   # Define custom input parameters using the kwargs dictionaries
   kwargs_time = {
      'threshold': 35,     # Valid key, will be used
      'nfft': 2**8         # Invalid key for the time domain, will be ignored
   }

   # Compute HRV parameters
   hrv(nni=nni, kwargs_time=kwargs_time)

This will trigger the following warning message.

.. warning::

   `Unknown kwargs for 'time_domain()': nfft. These kwargs have no effect.`
