Utils Module
============

The *Utils Module* contains general purpose functions (utilities) to support the features of the pyHRV toolbox incl. loading NNI sample data, check input data of HRV functions, segment arrays, check for duplicate file names to avoid accidental overwriting of files, and other utilities. These functions do not compute any HRV parameters nor provide any parameter-specific features (e.g. comparison plots).

.. seealso::

   `pyHRV Utils Module source code <https://github.com/PGomes92/pyhrv/blob/master/pyhrv/utils.py>`_

.. contents:: Module Contents

.. _ref-nni-sample:

Loading NNI Sample Data: load_sample_nni()
##########################################

.. py:function:: pyhrv.utils.load_sample_nni(series='short')

**Function Description**

Returns a short-term (5min) or long-term (60min) series of sample NNI found in the pyhrv/files/ directory.

These sample series were extracted from the `MIT-BIH NSRDB Database from physionet.org <https://physionet.org/physiobank/database/nsrdb/>`_ and can be found in the `samples <https://github.com/PGomes92/pyhrv/tree/master/samples>`_ folder.

**Input Parameters**
   - ``series`` (str): If 'long', returns a 60min NNI series, if 'short' returns a 5min NNI series (default: 'short').

**Returns**
   - ``nni`` (array): Series of NN intervals in [ms].

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv

   # Load short sample series (5min)
   nni = pyhrv.utils.load_sample_nni()

   # Load long sample series (60min)
   nni = pyhrv.utils.load_sample_nni(series='long')

Loading the pyHRV Keys: load_hrv_keys_json()
############################################

.. py:function:: pyhrv.utils.load_hrv_keys_json()

**Function Description**

Loads the content of the 'hrv_keys.json' file found in the 'pyhrv/files/' directory.

.. note::

   The *hrv_keys.json* file contains all the keys and additional information related to the computed HRV parameters. These keys are used throughout pyHRV when exporting and importing data.

   For example, the :code:`pyhrv.time_domain.sdnn()` function returns a *ReturnTuple* object from the BioSPPy package in which the result of the function can be accessed using the :code:`sdnn` key.

**Input Parameters**
   - none

**Returns**
   - ``hrv_keys`` (dict): Content of the pyhrv/files/hrv_keys.json file in a dictionary

**Example**

The following example code demonstrates how to use this function:

.. code-block:: python

   # Import packages
   import pyhrv

   # Load content of the hrv_keys.json file
   hrv_keys = pyhrv.utils.load_hrv_keys_json()

.. _ref-checkinput:

Check Input: check_input()
##########################

.. py:function:: pyhrv.utils.check_input(nn=None, rpeaks=None)

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

.. _ref-nnformat:

NN Format: nn_format()
######################

.. py:function:: pyhrv.utils.nn_format(nni=None)

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
   import pyhrv
   from opensignalsreader import OpenSignalsReader

   # Load sample ECG signal stored in an OpenSignals file
   signal = OpenSignalsReader('./samples/SampleECG.npy').signal('ECG')

   # Get R-peak locations
   rpeaks = biosppy.signals.ecg.ecg(signal)[2]

   # Compute NNI parameters
   nni = pyhrv.utils.nn_intervals(rpeaks)

   # Confirm [ms] format
   nni_in_ms = pyhrv.utils.nn_format(nni)


Check Interval: check_interval()
################################

.. py:function:: pyhrv.utils.check_interval(interval=None, limits=None, default=None)

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
   import pyhrv

   # Check valid interval limits; returns interval without modifications
   interval = [0, 10]
   res = pyhrv.utils.check_interval(interval)

   # Check invalid interval limits; returns corrected interval limits
   interval = [10, 0]
   res = pyhrv.utils.check_interval(interval)
   # here: res = [0, 10]

You can specify valid minimum and maximum values for the interval limits. If an interval with limits outside the valid
region are provided, the limits will be set to the specified valid minimum and maximum values:

.. code-block:: python

   # Specify minimum and maximum valid values (here: [2, 8]); interval is out of valid interval
   interval = [0, 10]
   limits = [2, 8]
   res = pyhrv.utils.check_interval(interval, limits)
   # here: res = [2, 8]

You can specify default values for this function. These can be used if no interval is specified by the user and default values should apply (e.g. when integrating this function in custom functions with dynamic intervals).

.. code-block:: python

   # Don't specify intervals or limits, but set a default values (here: [0, 10])
   res = pyhrv.utils.check_interval(interval=None, limits=None, default=[0, 10])

.. _ref-segmentation:

Segmentation: segmentation()
############################

.. py:function:: pyhrv.utils.segmentation(nn=None,rpeaks=None, overlap=False, duration=300)

**Function Description**

Segmentation of NNI series into individual segments of specified duration (e.g. splitting R-peak locations into 5 minute segments for computation of the SDNN index).

.. note::

   The segmentation of the NNI series can only be conducted if the sum of the NNI series (i.e. the maximum duration) is greater than the specified segment duration (``segment``).

   .. seealso::

      **Application Notes** below for more information.

**Input Parameters**
   - ``nni`` (array): NN interval series in [ms] or [s]
   - ``full`` (bool, optional): If True, returns last segment, even if the last segment is singificantly shorter than the specified duration (default: True)
   - ``duration`` (int, optional): Segment duration in [s] (default: 300s)
   - ``warn`` (bool, optional): If True, raise a warning message if a segmentation could not be conducted (duration > NNI series duration)

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
   import pyhrv

   # Load Sample NNI series (~5min)
   nni = pyhrv.utils.load_sample_nni()

   # Segment NNI series with a segment duration of [60s]
   segments, control = pyhrv.utils.segmentation(nn=nni, duration=60)

This will return 5 segments and the control variable will be ``True``. Use the code below to see the exact results:

.. code-block:: python

   # Print control variable
   print("Segmentation?", control)

   # Print segments
   for i, segment in enumerate(segments):
      print("Segment %i" % i)
      print(segment)

Join Tuples: join_tuples()
##########################

.. py:function:: pyhrv.utils.join_tuples(*args)

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
   import pyhrv

   # Join multiple ReturnTuple objects
   tuples = pyhrv.utils.join_tuples(return_tuple1, return_tuple2, return_tuple3)

Standard Deviation: std()
#########################

.. py:function:: pyhrv.utils.std(array=None, dof=1)

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
   import pyhrv

   # Sample array
   data = [600, 650, 800, 550, 900, 1000, 750]

   # Compute standard deviation
   sd = pyhrv.utils.std(data)
   # sd = 163.2993161855452

Time Vector: time_vector()
##########################

.. py:function:: pyhrv.utils.time_vector(signal=None, sampling_rate=1000.)

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
   import pyhrv
   from opensignalsreader import OpenSignalsReader

   # Load sample ECG signal stored in an OpenSignals file
   acq = OpenSignalsReader('./samples/SampleECG.npy')
   signal = acq.signal('ECG')
   sampling_rate = acq.sampling_rate

   # Compute time vector
   t = pyhrv.utils.time_vector(signal, sampling_rate)
