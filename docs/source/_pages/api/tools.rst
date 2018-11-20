Tools Module
============

The *Tools Module* contains general purpose functions and key functionalities (e.g. computation of NNI series) for the entire HRV package, among other useful features for HRV analysis workflow (e.g., HRV report, HRV export/import).

.. seealso::

   `pyHRV Tools Module source code <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/tools.py>`_

.. contents:: Module Contents

.. important::

   Some of the examples below use an ECG signal recorded with the OpenSignals (r)evolution software and are loaded with the ``opensignalsreader`` package.

   These examples do work with any other ECG signal independently of the acquisition software, of course.

   The sample NNI series used in the some examples below were taken from the NNI samples that come with the pyHRV
   package.

.. seealso::

   Useful links:

   * `OpenSignals (r)evolution software <http://bitalino.com/en/software>`_
   * `Sample ECG file acquired with the OpenSignals software <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/samples/SampleECG.txt>`_
   * :ref:`ref-samples` (docs)
   * `Sample NNI Series on GitHub <https://github.com/PGomes92/pyhrv/tree/master/pyhrv/samples>`_
   * `series_1.npy (file used in the examples below) <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/samples/series_1.npy>`_

NN Intervals: nn_intervals()
############################

.. py:function:: pyhrv.tools.nn_intervals(rpeaks=None)

**Function Description**

Computes the series of NN intervals [ms] from a series of successive R-peak locations.

**Input Parameters**
   - ``rpeaks`` (array): R-peak times in [ms] or [s].

**Returns**
   - ``nni`` (array): Series of NN intervals in [ms].

**Computation**

The NN interval series is computed from the R-peak series as follows:

.. math::

   NNI_{j} = R_{j+1} - R_{j}

for :math:`0 <= j <= (n - 1)`

with:

   * :math:`NNI_j`: NN interval j
   * :math:`R_j`: Current R-peak j
   * :math:`R_{j+1}`: Successive R-peak j+1
   * :math:`n`: Number of R-peaks

**Application Notes**

The ``nni`` series will be returned in [ms] format, even if the ``rpeaks`` are provided in [s] format.

.. seealso::

   :ref:`ref-nnformat` for more information about the [s] to [ms] conversion.

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import biosppy
   import pyhrv.tools as tools
   from opensignalsreader import OpenSignalsReader

   # Load sample ECG signal stored in an OpenSignals file
   signal = OpenSignalsReader('SampleECG.txt').signal('ECG')

   # Get R-peak locations (and hide the ECG plots)
   rpeaks = biosppy.signals.ecg.ecg(signal, show=False)[2]

   # Compute NNI parameters
   nni = tools.nn_intervals(rpeaks)

.. _ref-nnformat:

NN Format: nn_format()
######################

.. py:function:: pyhrv.tools.nn_format(nni=None)

**Function Description**

Checks the format of the NNI series and converts data in [s] to [ms] format. Additionally, it ensures that the data will be returned in the ``NumPy`` array format.

**Input Parameters**
   - ``nni`` (array): NNI series [ms] or [s].

**Returns**
   - ``nni`` (array): NNI series [ms] and NumPy array format.

**Computation**

The automatic [s] to [ms] conversion occurs on a threshold based identification whether the data is in [s] or [ms] format: if the maximum value of the input array is < 10, then the data is assumed to be in [s] format.

This conversion process is based on the following two assumptions:

   * Any interval data in [s] format ranges between 0.2s (:math:`\hat{=}300bpm`) and 1.5s (:math:`\hat{=}40bpm`). Any interval greater 1.5s is highly unlikely to occur, and even if it does, it does still not reach the specified maximum interval limit of 10s (:math:`\hat{=}6bpm`)
   * The provided NNI series has been filtered from NNI outliers caused by signal artifacts (e.g. ECG signal loss)

.. note::

   It is important to filter the NNI series from the intervals caused by signal artifacts first, otherwise the returned series will be influenced by these NNI and distort all HRV parameter results.

**Application Notes**

The ``nni`` series will be returned in [ms] format, even if the ``rpeaks`` are provided in [s] format.

.. seealso::

   :ref:`ref-nnformat` for more information about the [s] to [ms] conversion.

**Example**

The following example code demonstrates how to use this function:

.. note::

   This functions is intended to be used by the parameter functions of ``pyHRV``, an external use might not be appropriate.

.. code-block:: python

   # Import packages
   import biosppy
   import pyhrv.tools as tools
   from opensignalsreader import OpenSignalsReader

   # Load sample ECG signal stored in an OpenSignals file
   signal = OpenSignalsReader('./samples/SampleECG.npy').signal('ECG')

   # Get R-peak locations
   rpeaks = biosppy.signals.ecg.ecg(signal)[2]

   # Compute NNI parameters
   nni = tools.nn_intervals(rpeaks)

   # Confirm [ms] format
   nni_in_ms = tools.nn_format(nni)


NN Interval Differences: nn_diff()
##################################

.. py:function:: pyhrv.tools.nn_diff(nn=None)

**Function Description**

Computes the series of NN interval differences [ms] from a series of successive NN intervals.

**Input Parameters**
   - ``nni`` (array): NNI series in [ms] or [s].

**Returns**
   - ``nn_diff_`` (array): Series of NN interval differences in [ms].

**Computation**

The NN interval series is computed from the R-peak series as follows:

.. math::

   \Delta NNI_j = NNI_{j+1} - NNI_j

for :math:`0 <= j <= (n - 1)`

with:

   * :math:`\Delta NNI_j`: NN interval j
   * :math:`NNI_j`: Current NNI j
   * :math:`NNI_{j+1}`: Successive NNI j+1
   * :math:`n`: Number of NNI

**Application Notes**

The ``nn_diff_`` series will be returned in [ms] format, even if the ``nni`` are provided in [s] format.

.. seealso::

   :ref:`ref-nnformat` for more information about the [s] to [ms] conversion.

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import biosppy
   import pyhrv.tools as tools
   from opensignalsreader import OpenSignalsReader

   # Load sample ECG signal stored in an OpenSignals file
   signal = OpenSignalsReader('./samples/SampleECG.npy').signal('ECG')

   # Get R-peak locations
   rpeaks = biosppy.signals.ecg.ecg(signal)[2]

   # Compute NNI parameters
   nni = tools.nn_intervals(rpeaks)

   # Compute NNI differences
   delta_nni = tools.nn_diff(nni)

.. _ref-hr:

Heart Rate: heart_rate()
########################

.. py:function:: pyhrv.tools.heart_rate(nni=None, rpeaks=None)

**Function Description**

Computes a series of Heart Rate values in [bpm] from a series of NN intervals or R-peaks in [ms] or [s] or the HR from a single NNI.

**Input Parameters**
   - ``nni`` (int, float, array): NN interval series in [ms] or [s]
   - ``rpeaks`` (array): R-peak locations in [ms] or [s]

**Returns**
   - ``hr`` (array): Series of NN intervals in [ms].

**Computation**

The Heart Rate series is computed as follows:

.. math::

   HR_j = \frac{60000}{NNI_j}

for :math:`0 <= j <= n`

with:

   * :math:`HR_j`: Heart rate j (in [bpm])
   * :math:`NNI_j`: NN interval j (in [ms])
   * :math:`n`: Number of NN intervals

**Application Notes**

The input ``nni`` series will be converted to [ms], even if the ``rpeaks`` are provided in [s] format.

.. seealso::

   :ref:`ref-nnformat` for more information about the [s] to [ms] conversion.

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.tools as tools

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute Heart Rate series
   hr = tools.heart_rate(nn)

It is also possible to compute the HR from a single NNI:

.. code-block:: python

   # Compute Heart Rate from a single NNI
   hr = tools.heart_rate(800)
   # here: hr = 75 [bpm]

.. Attention::

   In this case, the input NNI must be provided in [ms] as the [s] to [ms] conversion is only conducted for series of NN Intervals.

Plot ECG: plot_ecg()
####################

.. py:function:: pyhrv.tools.plot_ecg(signal=None, t=None, samplin_rate=1000., interval=None, rpeaks=True, figsize=None, title=None, show=True)

**Function Description**

Plots ECG signal on a medical grade ECG paper-like figure layout.

An example of an ECG plot generated by this function can be seen here:

.. figure:: /_static/ecg10.png

The x-Division does automatically adapt to the visualized interval (e.g., 10s interval -> 1s, 20s interval -> 2s, ...).

**Input Parameters**
   - ``signal`` (array): ECG signal (filtered or unfiltered)
   - ``t`` (array, optional): Time vector for the ECG signal (default: None)
   - ``sampling_rate`` (int, float, optional): Sampling rate of the acquired signal in [Hz] (default: 1000Hz)
   - ``interval`` (array, optional): Visualization interval of the ECG signal plot (default: [0s, 10s])
   - ``rpeaks`` (bool, optional): If True, marks R-peaks in ECG signal (default: True)
   - ``figsize`` (array, optional): Matplotlib figure size (width, height) (default: None: (12, 4)
   - ``title`` (str, optional): Plot figure title (default: None)
   - ``show`` (bool, optional): If True, shows the ECG plot figure (default: True)

**Returns**
   - ``fig_ecg`` (matplotlib figure object): Matplotlibe figure of the ECG plot

**Application Notes**

The input ``nni`` series will be converted to [ms], even if ``nni`` are provided in [s] format.

.. seealso::

   :ref:`ref-nnformat` for more information about the [s] to [ms] conversion.

This functions marks, by default, the detected R-peaks. Use the ``rpeaks`` input parameter to turn on (``rpeaks=True``) or to turn of (``rpeaks=False``) the visualization of these markers.

.. important::

   This parameter will have no effect if the number of R-peaks within the visualization interval is greater than 50. In this case, for reasons of plot clarity, no R-peak markers will be added to the plot.

The time axis scaling will change depending on the duration of the visualized interval:

   * t in [s] if visualized duration <= 60s
   * t in [mm:ss] (minutes:seconds) if 60s < visualized duration <= 1h
   * t in [hh:mm:ss] (hours:minutes:seconds) if visualized duration > 1h

**Example**

.. code-block:: python

   # Import
   import pyhrv.tools as tools
   from opensignalsreader import OpenSignalsReader

   # Load ECG data
   signal = OpenSignalsReader('SampleECG.txt').signal('ECG')

   # Plot ECG
   tools.plot_ecg(signal)

The plot of this example should look like the following plot:

.. figure:: /_static/ecg10.png
   :align: center

   Default visualization interval of the ``plot_ecg()`` function.

Use the ``interval`` input parameter to change the visualization interval using a 2-element array (``[lower_interval_limit, upper_interval_limit]``; default: 0s to 10s). Additionally, use the ``rpeaks`` parameter to toggle the R-peak markers.

The following code sets the visualization interval from 0s to 20s and hides the R-peak markers:

.. code-block:: python

   # Plot ECG
   tools.plot_ecg(signal, interval=[0, 20], rpeaks=False)

The plot of this example should look like the following plot:

.. figure:: /_static/ecg20.png
   :align: center

   Visualizing the first 20 seconds of the ECG signal without R-peak markers.

Use the ``title`` input parameter to add titles to the ECG plot:

.. code-block:: python

   # Plot ECG
   tools.plot_ecg(signal, interval=[0, 20], rpeaks=False, title='This is a Title')

.. figure:: /_static/ecg20title.png
   :align: center

   ECG plot with custom title.

Tachogram: tachogram()
######################

.. py:function:: pyhrv.tools.tachogram(signal=None, nn=None,rpeaks=None, sampling_rate=1000., hr=True, interval=None, title=None, figsize=None, show=True)

**Function Description**

Plots Tachogram (NNI & HR) of an ECG signal, NNI or R-peak series.

An example of a Tachogram plot generated by this function can be seen here:

.. figure:: /_static/tachogram10.png

**Input Parameters**
   - ``signal`` (array): ECG signal (filtered or unfiltered)
   - ``nni`` (array): NN interval series in [ms] or [s]
   - ``rpeaks`` (array): R-peak locations in [ms] or [s]   - ``t`` (array, optional): Time vector for the ECG signal (default: None)
   - ``sampling_rate`` (int, optional): Sampling rate in [hz] of the ECG signal (default: 1000Hz)
   - ``hr`` (bool, optional): If True, plot HR seres in [bpm] on second axis (default: True)
   - ``interval`` (array, optional): Visualization interval of the Tachogram plot (default: None: [0s, 10s])
   - ``title`` (str, optional): Optional plot figure title (default: None)
   - ``figsize`` (array, optional): Matplotlib figure size (width, height) (default: None: (12, 4))
   - ``show`` (bool, optional): If True, shows the ECG plot figure (default: True)

**Returns**
   - ``fig`` (matplotlib figure object): Matplotlib figure of the Tachogram plot.

**Application Notes**

The input ``nni`` series will be converted to [ms], even if the ``rpeaks`` or ``nni`` are provided in [s] format.

.. seealso::

   :ref:`ref-nnformat` for more information about the [s] to [ms] conversion.

**Example**

The following example demonstrates how to load an ECG signal recorded with the OpenSignals (r)evolution and loaded with the `opensignalsreader` package.

.. code-block:: python

   # Import
   import pyhrv.tools as tools
   from opensignalsreader import OpenSignalsReader

   # Load ECG data
   signal = OpenSignalsReader('SampleECG.txt').signal('ECG')

   # Plot ECG
   tools.plot_ecg(signal)

Alternatively, use R-peak data to plot the histogram...

.. code-block:: python

   # Import
   import biosppy
   import pyhrv.tools as tools
   from opensignalsreader import OpenSignalsReader

   # Load ECG data
   signal = OpenSignalsReader('SampleECG.txt').signal('ECG')

   # Extract R-peaks
   rpeaks = biosppy.signals.ecg.ecg(signal)[2]

   # Plot ECG
   tools.tachogram(rpeaks=rpeaks)

... or directly the NNI series...

.. code-block:: python

   # Compute NNI intervals from the R-peaks
   nni = tools.nn_intervals(rpeaks)

   # Plot ECG
   tools.tachogram(nni=nni)

The plots generated by the examples above should look like the plot below:

.. figure:: /_static/tachogram10.png
   :align: center

   Tachogram with default visualization interval.

Use the ``interval`` input parameter to change the visualization interval (default: 0s to 10s; here: 0s to 20s):

.. code-block:: python

   # Plot ECG
   tools.plot_ecg(signal, interval=[0, 20])

The plot of this example should look like the following plot:

.. figure:: /_static/tachogram20.png
   :align: center

   Tachogram with custom visualization interval.

.. note::

   Interval limits which are out of bounce will automatically be corrected.

   Example:
      * lower limit < 0 -> lower limit = 0
      * upper limit > maximum ECG signal duration -> upper limit = maximum ECG signal duration


Set the ``hr`` parameter to ``False`` in case only the NNI Tachogram is needed:

.. code-block:: python

   # Plot ECG
   tools.plot_ecg(signal, interval=[0, 20], hr=False)

.. figure:: /_static/tachogramNoHR.png
   :align: center

   Tachogram of the NNI series only.

The time axis scaling will change depending on the duration of the visualized interval:

   * t in [s] if visualized duration <= 60s
   * t in [mm:ss] (minutes:seconds) if 60s < visualized duration <= 1h
   * t in [hh:mm:ss] (hours:minutes:seconds) if visualized duration > 1h

.. figure:: /_static/tachogramlong.png
   :align: center

   Tachogram of an ~1h NNI series.

Check Interval: check_interval()
################################

.. py:function:: pyhrv.tools.check_interval(interval=None, limits=None, default=None)

**Function Description**

General purpose function that checks and verifies correctness of interval limits within optionally defined valid interval specifications and/or default values if no interval is specified.

This function can be used to set visualization intervals, check overlapping frequency bands, or for other similar purposes, and is intended to automatically catch possible error sources due to invalid intervals boundaries.

**Input Parameters**
   - ``interval`` (array): Input interval [min, max] (default: None)
   - ``limits`` (array): Minimum and maximum allowed interval limits (default: None)
   - ``default`` (array): Specified default interval (e.g. if ``interval`` is None) (default: None)

**Returns**
   - ``interval`` (array): Interval with correct(ed)/valid interval limits.

**Raises**
   - ``TypeError`` If no input data is specified.
   - ``ValueError`` If the input interval[s] have equal lower and upper limits.

**Computation**

The input data is provided as ``interval = [int(lower_limit), int(upper_limit)]``. Depending on the limits, the following conditions should be met:

   * If ``lower_limit > upper_limit``: the interval limits will be switched to ``interval = [upper_limit, lower_limit]``
   * If ``lower_limit == upper_limit``: raises ``ValueError``

If minimum and maximum intervals are specified, i.e. ``limit = [int(minimum), int(maximum)]``, the following additional actions may occur:

   * If ``lower_limit < minimum``: the lower limit will be set to the minimum allowed limit ``lower_limit = minimum``
   * If ``upper_limit > maximum``: the upper limit will be set to the maximum allowed limit ``upper_limit = maximum``

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv.tools as tools

   # Check valid interval limits; returns interval without modifications
   interval = [0, 10]
   res = tools.check_interval(interval)

   # Check invalid interval limits; returns corrected interval limits
   interval = [10, 0]
   res = tools.check_interval(interval)
   # here: res = [0, 10]

You can specify valid minimum and maximum values for the interval limits. If an interval with limits outside the valid
region are provided, the limits will be set to the specified valid minimum and maximum values:

.. code-block:: python

   # Specify minimum and maximum valid values (here: [2, 8]); interval is out of valid interval
   interval = [0, 10]
   limits = [2, 8]
   res = tools.check_interval(interval, limits)
   # here: res = [2, 8]

You can specify default values for this function. These can be used if no interval is specified by the user and default values should apply (e.g. when integrating this function in custom functions with dynamic intervals).

.. code-block:: python

   # Don't specify intervals or limits, but set a default values (here: [0, 10])
   res = tools.check(interval=None, limits=None, default=[0, 10])

.. _ref-segmentation:

Segmentation: segmentation()
############################

.. py:function:: pyhrv.tools.segmentation(nn=None,rpeaks=None, overlap=False, duration=300)

**Function Description**

Segmentation of NNI series into individual segments of specified duration (e.g. splitting R-peak locations into 5 minute segments for computation of the SDNN index).

.. note::

   The segmentation of the NNI series can only be conducted if the sum of the NNI series (i.e. the maximum duration) is greater than the specified segment duration (``segment``).

   .. seealso::

      **Application Notes** below for more information.

**Input Parameters**
   - ``nni`` (array): NN interval series in [ms] or [s]
   - ``full`` (bool, optional): If True, returns last segment, even if the last segment is singificantly shorter than the specified duration (default: True)
   - ``overlap`` (bool, optional): If True, allow to return NNI that go from the interval of one segment to the successive segment (default: False)
   - ``duration`` (int, optional): Segment duration in [s] (default: 300s)

**Returns**
   - ``segments`` (array of arrays): Array with the segmented NNI series.
   - ``control`` (bool): If True, segmentation was possible.

.. seealso::

   **Application Notes** below for more information about the returned segmentation results.

**Raises**
   - ``TypeError``: If ``nni`` input data is not specified

**Application Notes**

The function returns the results in an array of arrays if a segmentation of the signal was possible. This requires the sum of the provided NNI series (i.e. the maximum duration) to be greater than the specified segment duration (``segment``). In this case, a segmentation can be conducted and the segments with the respective NNIs will be returned along with the control variable set to ``True``.

If a segmentation cannot be conducted, i.e. the maximum duration of the NNI series is shorter than the specified segment duration, the input unmodified NNI series will be returned along with the control variable set to ``False``.

You can use the control variable to test whether the segmentation could be conducted successfully or not.

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.tools as tools

   # Load Sample NNI series (~5min)
   nni = np.load('series_1.npy')

   # Segment NNI series with a segment duration of [60s]
   segments, control = tools.segmentation(nn=nni, duration=60)

This will return 5 segments and the control variable will be ``True``. Use the code below to see the exact results:

.. code-block:: python

   # Print control variable
   print("Segmentation?", control)

   # Print segments
   for i, segment in enumerate(segments):
      print("Segment %i" % i)
      print(segment)

HRV Reports: hrv_report()
#########################

.. py:function:: pyhrv.tools.hrv_report(results=None, path=None, rfile=None, nn=None, info={}, file_format='txt', delimiter=';', hide=False, plots=False)

**Function Description**

Generates HRV report (in .txt or .csv format) of the provided HRV results. You can find a sample report generated with this function `here <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/files/SampleReport.txt>`_.

**Input Parameters**
   - ``results`` (dict, ReturnTuple object): Computed HRV parameter results
   - ``path`` (str): Absolute path of the output directory
   - ``rfile`` (str): Output file name
   - ``nni`` (array, optional): NN interval series in [ms] or [s]
   - ``info`` (dict, optional): Dictionary with HRV metadata
   - ``file_format`` (str, optional): Output file format, select 'txt' or 'csv' (default: 'txt')
   - ``delimiter`` (str, optional): Delimiter separating the columns in the report (default: ';')
   - ``hide`` (bool, optional): Hide parameters in report that have not been computed
   - ``plots`` (bool, optional): If True, save plot figures in .png format

.. note::

   The ``info`` dictionary can contain the following metadata:

      * key: ``file`` - Name of the signal acquisition file
      * key: ``device`` - ECG acquisition device
      * key: ``identifier`` - ECG acquisition device identifier (e.g. MAC address)
      * key: ``fs`` - Sampling rate used during ECG acquisition
      * key: ``resolution`` - Resolution used during acquisition

   Any other key will be ignored.

.. important::

   It is recommended to use absolute file paths when using the ``path`` parameter to ensure the correct functionality of this function.

**Raises**
   - ``TypeError``: If no HRV results are provided
   - ``TypeError``: If no file or directory path is provided
   - ``TypeError``: If the specified selected file format is not supported
   - ``IOError``: If the selected output file or directory does not exist

**Application Notes**

This function uses the weak ``_check_fname()`` function found in this module to prevent the (accidental) overwriting of existing HRV reports. If a file with the file name ``rfile`` does exist in the specified ``path``, then the file name will be incremented.

For instance, if a report file with the name  *SampleReport.txt* exists, this file will not be overwritten, instead, the file name of the new report will be incremented to *SampleReport_1.txt*.

If the file with the file name *SampleReport_1.txt* exists, the file name of the new report will be incremented to *SampleReport_2.txt*, and so on...

.. important::

   The maximum supported number of file name increments is limited to 999 files, i.e., using the example above, the
   implemented file protection mechanisms will go up to *SampleReport_999.txt*.

If no file name is provided, an automatic file name with a time stamp will be generated for the generated report
(*hrv_report_YYYY-MM-DD_hh-mm-ss.txt*  or *hrv_report_YYYY-MM-DD_hh-mm-ss.txt*).

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv
   import numpy as np
   import pyhrv.tools as tools

   # Load Sample NNI series (~5min)
   nni = np.load('series_1.npy')

   # Compute HRV results
   results = pyhrv.hrv(nn=nni)

   # Create HRV Report
   tools.hrv_report(results, rfile='SampleReport', path='/my/favorite/path/')


This generates a report looking like the one below:

.. figure:: /_static/samplereport.png
   :scale: 50%

.. seealso::

   * `Sample report in .txt format <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/files/SampleReport.txt>`_
   * `Sample report in .csv format <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/files/SampleReport.csv>`_

.. _ref-hrvexport:

HRV Export: hrv_export()
########################

.. py:function:: pyhrv.tools.hrv_export(results=None, path=None, efile=None, comment=None, plots=False)

**Function Description**

Exports HRV results into a JSON file. You can find a sample export generated with this function `here <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/files/SampleExport.json>`_.

**Input Parameters**
   - ``results`` (dict, ReturnTuple object): Computed HRV parameter results
   - ``path`` (str): Absolute path of the output directory
   - ``efile`` (str): Output file name
   - ``comment`` (str, optional): Optional comment
   - ``plots`` (bool, optional): If True, save figures of the results in .png format

.. important::

   It is recommended to use absolute file paths when using the ``path`` parameter to ensure the correct operation of this function.

**Returns**
   - ``efile`` (str): Absolute path of the output report file (may vary from the input data)

**Raises**
   - ``TypeError``: If no HRV results are provided
   - ``TypeError``: If no file or directory path is provided
   - ``TypeError``: If specified selected file format is not supported
   - ``IOError``: If the selected output file or directory does not exist

**Application Notes**

This function uses the weak ``_check_fname()`` function found in this module to prevent the (accidental) overwriting of existing HRV exports. If a file with the file name ``efile`` exists in the specified ``path``, then the file name will be incremented.

For instance, if an export file with the name  *SampleExport.json* exists, this file will not be overwritten, instead,
the file name of the new export file will be incremented to *SampleExport_1.json*.

If the file with the file name *SampleExport_1.json* exists, the file name of the new export will be incremented to
*SampleExport_2.json*, and so on.

.. important::

   The maximum supported number of file name increments is limited to 999 files, i.e., using the example above, the
   implemented file protection mechanisms will go up to *SampleExport_999.json*.

If no file name is provided, an automatic file name with a time stamp will be generated for the generated report
(*hrv_export_YYYY-MM-DD_hh-mm-ss.json*).

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv
   import numpy as np
   import pyhrv.tools as tools

   # Load Sample NNI series (~5min)
   nni = np.load('series_1.npy')

   # Compute HRV results
   results = pyhrv.hrv(nn=nni)

   # Export HRV results
   tools.hrv_export(results, efile='SampleExport', path='/my/favorite/path/')


.. seealso::

   * `Sample HRV export <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/files/SampleExport.json>`_

HRV Import: hrv_import()
########################

.. py:function:: pyhrv.tools.hrv_import(hrv_file=None)

**Function Description**

Imports HRV results stored in JSON files generated with the 'hrv_export()'.

.. seealso::

   * :ref:`ref-hrvexport` function
   * `Sample HRV export <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/files/SampleExport.json>`_

**Input Parameters**
   - ``hrv_file`` (str, file handler): File handler or absolute string path of the HRV JSON file

**Returns**
   - ``output`` (ReturnTuple object): All HRV parameters stored in a ``biosppy.utils.ReturnTuple`` object

**Raises**
   - ``TypeError``: If no file path or handler is provided

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv.tools as tools

   # Import HRV results
   hrv_results = tools.hrv_import('/path/to/my/HRVResults.json')

.. seealso::

   `HRV keys file <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/files/hrv_keys.json>`_ for a full list of HRV parameters and their respective keys.

Join Tuples: join_tuples()
##########################

.. py:function:: pyhrv.tools.join_tuples(*args)

**Function Description**

Joins multiple biosppy.utils.ReturnTuple objects into one biosppy.utils.ReturnTuple object.

.. seealso::

   :ref:`ref-returntuple`

**Input Parameters**
   - ``*args`` (biosppy.utils.ReturnTuple): Multiple biosppy.utils.ReturnTuple objects (can also be stored in an array)

**Returns**
   - ``output`` (ReturnTuple object): biosppy.utils.ReturnTuple object with the content of all input tuples/objects merged together.

**Raises**
   - ``TypeError``: If no input data is provided
   - ``TypeError``: If input data contains non-biosppy.utils.ReturnTuple objects

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv.tools as tools

   # Join multiple ReturnTuple objects
   tuples = tools.join_tuples(return_tuple1, return_tuple2, return_tuple3)

Standard Deviation: std()
#########################

.. py:function:: pyhrv.tools.std(array=None, dof=1)

**Function Description**

Computes the standard deviation of a data series.

**Input Parameters**
   - ``array`` (array): Data series
   - ``dof`` (int, optional): Degrees of freedom (default: 1)

**Returns**
   - ``result`` (float): Standard deviation of the input data series

**Raises**
   - ``TypeError``: If no input array is provided

**Computation**

The standard deviation is computed according to the following formula:

.. math::

   SD = \sqrt{\frac{1}{n-dof} \sum_{i=1}^{n} (NNI_j - \overline{NNI})^2}

with:
   * :math:`SD`: Standard Deviation
   * :math:`n`: Number of NN Intervals
   * :math:`dof`: Degrees of Freedom
   * :math:`NNI_j`: NN Interval j
   * :math:`\overline{NNI}`: Mean NN Interval

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv.tools as tools

   # Sample array
   data = [600, 650, 800, 550, 900, 1000, 750]

   # Compute standard deviation
   sd = tools.std(data)
   # sd = 163.2993161855452

Time Vector: time_vector()
##########################

.. py:function:: pyhrv.tools.time_vector(signal=None, sampling_rate=1000.)

**Function Description**

Computes time vector based on the sampling rate of the provided input signal.

**Input Parameters**
   - ``signal`` (array): ECG signal (or any other sensor signal)
   - ``sampling_rate`` (int, float, optional): Sampling rate of the input signal in [Hz] (default: 1000Hz)

**Returns**
   - ``time_vector`` (array): Time vector for the input signal sampled at the input ``sampling_rate``

**Raises**
   - ``TypeError``: If no input array is provided

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv.tools as tools
   from opensignalsreader import OpenSignalsReader

   # Load sample ECG signal stored in an OpenSignals file
   acq = OpenSignalsReader('./samples/SampleECG.npy')
   signal = acq.signal('ECG')
   sampling_rate = acq.sampling_rate

   # Compute time vector
   t = tools.time_vector(signal, sampling_rate)

.. _ref-checkinput:

Check Input: check_input()
##########################

.. py:function:: pyhrv.tools.check_input(nn=None, rpeaks=None)

**Function Description**

Checks if input series of NN intervals or R-peaks are provided and, if yes, returns a NN interval series in [ms] format.

**Input Parameters**
   - ``nni`` (array): NN interval series in [ms] or [s] (default: None)
   - ``rpeaks`` (array): R-peak locations in [ms] or [s] (default: None)

**Returns**
   - ``nni`` (array): NN interval series in [s] (default: None)

**Raises**
   - ``TypeError``: If no R-peak data or NN intervals provided

**Application Notes**

This function is mainly used by the parameter computation functions of the ``time_domain.py``, the ``frequency_domain.py``, and the ``nonlinear.py`` modules.

The ``nni`` series will be returned in [ms] format, even if the input ``rpeaks`` or ``nni`` are provided in [s] format.

.. seealso::

   :ref:`ref-nnformat` for more information about the [s] to [ms] conversion.
