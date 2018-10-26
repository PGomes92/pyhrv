Time Domain Module
==================
The ``time_domain.py`` module contains all the functions to compute the HRV time domain parameters. The individual parameter level functions are listed in the table below.

.. contents:: Module Contents


NNI Parameters
##############

.. py:function:: pyhrv.time_domain.nn_parameters(nn=None, rpeaks=None)

**Function Description:**

Computes basic statistical parameters from a series of NN intervals (# of intervals, mean, min, max).

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object):**
   - ``nn_counter`` (int): Number of NNI (-) [key: 'nn_counter']
   - ``nn_mean`` (float): Mean NNI (ms) [key: 'nn_mean']
   - ``nn_min`` (int): Minimum NNI (ms) [key: 'nn_min']
   - ``nn_max`` (int): Maximum NNI (ms) [key: 'nn_max']

**Notes**: If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**:
The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute NNI parameters
   results = td.nn_parameters(nn)

   # Get & print minimum NNI
   print(results['nn_min'])


∆NNI Parameters
###############

.. py:function:: pyhrv.time_domain.nn_differences_parameters(nn=None, rpeaks=None)

**Function Description:**

Computes basic statistical parameters from a series of NN interval differences (# of intervals, mean, min, max).

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object):**
   - ``nn_diff_mean`` (float): Mean NNI difference (ms) [key: 'nn_diff_mean']
   - ``nn_diff_min`` (int): Minimum NNI difference (ms) [key: 'nn_diff_min']
   - ``nn_diff_max`` (int): Maximum NNI difference (ms) [key: 'nn_diff_max']

**Notes**:If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**:
The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute NNI differences parameters
   results = td.nn_differences_parameters(nn)

   # Get & print minimum NNI difference
   print(results['nn_diff_min'])

Heart Rate Parameters
#####################

.. py:function:: pyhrv.time_domain.hr_parameters(nn=None, rpeaks=None)

**Function Description:**

Computes basic statistical parameters from a series of heart rate (HR) data (mean, min, max, standard deviation)

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object):**
   - ``hr_mean`` (float): Mean heart rate (bpm) [key: 'hr_mean']
   - ``hr_min`` (int): Minimum heart rate (bpm) [key: 'hr_min']
   - ``hr_max`` (int): Maximum heart rate (bpm) [key: 'hr_max']
   - ``hr_std`` (float): Standard deviation of the heart rate series (bpm) [key: 'hr_std']

**Notes**:If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**:
The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute HR parameters
   results = td.hr_parameters(nn)

   # Get & print minimum HR difference
   print(results['hr_min'])


.. _ref-sdnn:

SDNN
####
.. py:function:: pyhrv.time_domain.sdnn(nn=None, rpeaks=None)

**Function Description:**

Computes the Standard Deviation of an NN interval series (SDNN).

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object):**
   - ``sdnn`` (float): Standard deviation of NN intervals (ms) [key: 'sdnn']

**Parameter Computation**

The SDNN parameter is computed according to the following formula:

.. math::

   SDNN = \sqrt{\frac{1}{n - 1} \sum_{j=1}^{n} (NNI_j - \overline{NNI})^2}

with:
   * n: Number of NNI
   * NNI_j: NNI j
   * :math:`\overline{NNI}`: Mean of NNI series

**Notes**: If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**:
The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute SDNN parameter
   results = td.sdnn(nn)

   # Get & print SDNN
   print(results['sdnn'])

SDNN Index
##########
.. py:function:: pyhrv.time_domain.sdnn_index(nn=None, rpeaks=None, full=False, duration=300)

**Function Description**

Computes the SDNN Index of an NNI series with a specified segmentation duration of ``duration`` (300s=5min by default).

**Input Parameters**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).
   - ``full`` (bool): If true, returns the last segment even if its duration is significantly shorter than ``duration`` (default: ``False``).
   - ``duration`` (int): Maximum duration per segment in (s) (default: 300s)

   .. note:: ``full`` is ``False`` by default which causes the last segment to be dropped. For instance, if processing an NNI series of 12.5min and the default segment duration of 5min, the segmentation function would split this series into 3 segments of 5min, 5min and 2.5min in duration. In this case, the last segment greatly alter the SDNN Index. Set the ``full`` parameter to ``False`` to drop the last segment or to ``True`` to compute the SDNN Index even with shorter segments.

**Returns (ReturnTuple Object)**
   - ``sdnn_index`` (float): SDNN Index (ms) [key: 'sdnn_index']

**Parameter Computation**

The SDNN Index is computed using the ``pyhrv.time_domain.sdnn()`` (:ref:`ref-sdnn`) and the ``pyhrv.tools.segmentation()`` functions. First, the input NNI series is segmented into segments of ~5min in duration. Second, the SDNN parameter of each segment is computed. Finally, the mean value of all computed SDNN values is computed.

These steps are presented in the flow chart below.

.. figure:: /_static/sdnn_index.png
   :align: center

**Application Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute SDNN Index parameter
   results = td.sdnn_index(nn)

   # Get & print SDNN Index
   print(results['sdnn_index'])

SDANN
#####
.. py:function:: pyhrv.time_domain.sdann(nn=None, rpeaks=None, full=False, duration=300)

**Function Description**

Computes the SDANN of an NNI series with a specified segmentation duration of ``duration`` (300s=5min by default).

**Input Parameters**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).
   - ``full`` (bool): If true, returns the last segment even if it's significantly shorter than ``duration`` (default: :code:`False`).
   - ``duration`` (int): Maximum duration per segment in (s) (default: 300s)

   .. note:: ``full`` is ``False`` by default which causes the last segment to be dropped. For instance, if processing an NNI series of 12.5min and the default segment duration of 5min, the segmentation function would split this series into 3 segments of 5min, 5min and 2.5min in duration. In this case, the last segment greatly alter the SDNN Index. Set the ``full`` parameter to ``False`` to drop the last segment or to ``True`` to compute the SDNN Index even with shorter segments.

**Returns (ReturnTuple Object)**
   - ``sdann`` (float): SDANN (ms) [key: 'sdann']

**Parameter Computation**

The SDANN is computed using the ``pyhrv.time_domain.sdnn()`` and the ``pyhrv.tools.segmentation()`` functions. First, the input NNI series is segmented into segments of ~5min in duration. Second, the mean of each segment is computed. Finally, the SDNN all computed mean values is computed.

These steps are presented in the flow chart below.

.. figure:: /_static/sdann.png
   :align: center

**Application Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute SDANN parameter
   results = td.sdann(nn)

   # Get & print SDANN
   print(results['sdann'])

RMSSD
#####
.. py:function:: pyhrv.time_domain.rmssd(nn=None, rpeaks=None)

**Function Description**

Computes the root mean of squared NNI differences.

**Input Parameters**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object)**
   - ``rmssd`` (float): Root mean of squared NNI differences (ms) [key: 'rmssd']

**Parameter Computation**

The RMSSD parameter is computed according to the following formula:

.. math::

   RMSSD = \sqrt{\frac{1}{n - 1} \sum_{j=1}^{n} \Delta {NNI_j}^2}

with:
   * :math:`n`: Number of NNI
   * :math:`\Delta NNI_j`: NNI differences

**Notes**
If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute RMSSD parameter
   results = td.rmssd(nn)

   # Get & print RMSSD
   print(results['rmssd'])

.. _ref-sdsd:

SDSD
#####
.. py:function:: pyhrv.time_domain.rmssd(nn=None, rpeaks=None)

**Function Description**

Standard deviation of NNI differences.

**Input Parameters**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object)**
   - key: ``sdsd`` (float): Standard deviation of NNI differences (ms) [key: 'sdsd']

**Parameter Computation**

The SDSD parameter is computed according to the following formula:

.. math::

   SDSD = \sqrt{\frac{1}{n - 1} \sum_{j=1}^{n} (\Delta {NNI_j} - \overline{\Delta NNI}^2}

with:
   * :math:`n`: Number of NNI
   * :math:`\Delta NNI_j`: NNI differences
   * :math:`\overline{NNI}`: Mean NNI

**Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute SDSD parameter
   results = td.sdsd(nn)

   # Get & print SDSD
   print(results['sdsd'])

NNXX
####

.. py:function:: pyhrv.time_domain.nnXX(nn=None, rpeaks=None, threshold=None)

**Function Description**

Finds number of NN interval differences greater than a specified threshold and ratio between number of intervals > threshold and total number of NN interval differences.

**Input Parameters**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).
   - ``threshold`` (int): hreshold for nnXX values in (ms).

**Returns (ReturnTuple Object)**
   - ``nnXX`` (int): Number of NN interval differences greater than the specified threshold (-) [key: 'nnXX' with XX = specified threshold]
   - ``pnnXX`` (float): Ratio between nnXX and total number of NN interval differences (-) [key: 'pnnXX' with XX = specified threshold]

**Exceptions:**
   - ``TypeError``: If no threshold is specified.
   - ``ValueError``: Threshold <= 0.

**Parameter Computation**

[TEXT]

**Application Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute nnXX parameters with a threshold of 42ms
   results = td.nnXX(nn, threshold=42)

   # Get & print nnXX (here: nn42) value
   print(results['nn42'])

NN50
####

.. py:function:: pyhrv.time_domain.nn50(nn=None, rpeaks=None)

**Function Description:**

Find number of NN interval differences which are greater 50ms (NN50) and ratio between NN50 and total amount of NN intervals.

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object):**
   - ``nn50`` (int): Number of NN interval differences greater 50ms [key: 'nn50']
   - ``pnn50`` (float): Ratio between NN50 and total number of NN intervals (ms) [key: 'pnn50']

**Parameter Computation**

[TEXT]

**Application Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute nn50 parameters
   results = td.nn50(nn)

   # Get & print nn50
   print(results['nn50'])

NN20
####

.. py:function:: pyhrv.time_domain.nn20(nn=None, rpeaks=None)

**Function Description:**

Find number of NN interval differences which are greater 20ms (NN20) and the ratio between NN20 and total amount of NN intervals.

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object):**
   - ``nn20`` (int): Number of NN interval differences greater 20ms [key: 'nn20']
   - ``pnn20`` (float): Ratio between NN20 and total number of NN intervals (ms) [key: 'pnn20']

**Parameter Computation**

The ``nn20`` parameter is the number of any ∆NNI > 20ms. The ``pnn20`` parameter is the ratio between the ``nn20`` parameter and the total number of ∆NNI.

**Notes**: If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

**Example**:
The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute nn20 parameters
   results = td.nn20(nn)

   # Get & print nn20
   print(results['nn20'])


Geometrical Parameters
######################

The geometrical parameters are computed based on the NNI histogram distribution. The TINN and Triangular Index are, in the most cases, provided together. However, ``pyHRV`` provides individual functions to individually compute the TINN (``pyhrv.time_domain.tinn()``) and Triangular Index (``pyhrv.time_domain.triangular_index()``) parameters. Additionally, the ``pyhrv.time_domain.geometrical_parameters()`` function allows you to compute all geometrical parameters and to join them in a single NNI histogram with using only a single function.

.. _ref-tinn:

TINN
----

.. py:function:: pyhrv.time_domain.tinn(nn=None, rpeaks=None, binsize=7.8125, plot=True, show=True, figsize=None, legend=True)

**Function Description:**

This function fits an interpolated triangle into the NNI histogram and computes its baseline width. See *Parameter
Computation* below for detailed information about the computation. As result, an NNI histogram (plot) as shown below is
computed.

.. figure:: /_static/tinn.png
   :align: center
   :scale: 40%

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).
   - ``binsize`` (int, float, optional): Bin size of the histogram bins (default: 7.8125ms).
   - ``plot`` (bool, optional): If True, create the histogram plot figure using ``matplotlib``. If False, the histogram data is computed using ``numpy`` with generating a histogram plot figure (default: True).
   - ``show`` (bool, optional): If True, shows the histogram plot figure (default: True).
   - ``figsize`` (array, optional): 2-element array with the ``matplotlib`` figure size ``figsize``. Format: ``figsize=(width, height)`` (default: will be set to (6, 6) if input is None).
   - ``legend`` (bool, optional): If True, adds legend to the histogram plot figure (default: True).

.. note::

   The ``binsize`` is pre-defined at 7.8125ms and is determined from the minimum suitable sampling frequency for ECG signals of 128Hz as recommended by the `HRV Guidelines <https://www.ahajournals.org/doi/full/10.1161/01.cir.93.5.1043>`_. At this sampling frequency, the temporal resolution of the signal used to derive NNI series is limited at 7.8125ms (= 1/128Hz).

**Returns (ReturnTuple Object)**

The results of this function are returned in a ``biosppy.utils.ReturnTuple`` object (see also :ref:`ref-returntuple`. Use the following keys below (on the left) to index the results.

   - ``tinn_histogram`` (matplotlib figure object): Histogram plot figure (only if input parameter ``plot`` is True
   - ``tinn_n`` (float): N value of the TINN computation (left corner of the interpolated triangle at (N, 0))
   - ``tinn_m`` (float): M value of the TINN computation (right corner of the interpolated triangle at (M, 0))
   - ``tinn`` (float): TINN (baseline width of the interpolated triangle) [ms]

**Parameter Computation**

The TINN parameters are computed based on the interpolation of a triangle into the NNI distribution. The positioning of the triangle's edges are determined by the following procedure: The first edge is positioned at the point *(D(X), X)* with *D(X)* being the histogram's maximum and *X* the bin containing the maximum. The other two edges are positioned at the points *(N, 0)* and *(M, 0)*. Finally, *N* and *M* are determined by finding the interpolated triangle with the best fit to the NNI histogram using the least squares method, as presented by the following formula:

.. math::

   E(n, N, M) = min{\sum (D(X) - q(n, N, M))^2}

with:
   * :math:`E(n)`: Error of the triangular interpolation with the best fit to the distribution
   * :math:`D(X)`: NNI distribution
   * :math:`q(n, N, m)`: Triangular interpolation function
   * :math:`n`: Bin
   * :math:`N`: N value determining the left corner of the interpolated triangle
   * :math:`M`: M value determining the right corner of the interpolated triangle

The main flow of this function is presented in the following flowchart:

.. figure:: /_static/tinn_flowchart.png
   :align: center

**Application Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

Use the ``legend`` input parameter do show or hide the legend in the histogram figure.

The ``show`` parameter only has effect if ``plot`` is set to True. If ``plot`` is False, no plot figure will be
generated, therefore, no figure can be shown using the ``show`` parameter.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute TINN parameters
   results = td.tinn(nn)

   # Get & print TINN and the N value
   print(results['tinn'])
   print(results['tinn_n'])

.. _ref-triindex:

Triangular Index
----------------

.. py:function:: pyhrv.time_domain.triangular_index(nn=None, rpeaks=None, binsize=7.8125, plot=True, show=True, figsize=None,legend=True)

**Function Description:**

Computes the triangular index based on the NN interval histogram.

.. figure:: /_static/trindex.png
   :align: center
   :scale: 40%

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).
   - ``binsize`` (int, float, optional): Bin size of the histogram bins (default: 7.8125ms).
   - ``plot`` (bool, optional): If True, create the histogram plot figure using ``matplotlib``. If False, the histogram data is computed using ``numpy`` with generating a histogram plot figure (default: True).
   - ``show`` (bool, optional): If True, shows the histogram plot figure (default: True).
   - ``figsize`` (array, optional): 2-element array with the ``matplotlib`` figure size ``figsize``. Format: ``figsize=(width, height)`` (default: will be set to (6, 6) if input is None).
   - ``legend`` (bool, optional): If True, adds legend to the histogram plot figure (default: True).

.. note::

   The ``binsize`` is pre-defined at 7.8125ms and is determined from the minimum suitable sampling frequency for ECG signals of 128Hz as recommended by the `HRV Guidelines <https://www.ahajournals.org/doi/full/10.1161/01.cir.93.5.1043>`_. At this sampling frequency, the temporal resolution of the signal used to derive NNI series is limited at 7.8125ms (= 1/128Hz).

**Returns (ReturnTuple Object):**

The results of this function are returned in a ``biosppy.utils.ReturnTuple`` object (see also :ref:`ref-returntuple`. Use the following keys below (on the left) to index the results.

   - ``tinn_histogram`` (matplotlib figure object): Histogram plot figure (only if input parameter ``plot`` is True
   - ``tinn_n`` (float): N value of the TINN computation (left corner of the interpolated triangle at (N, 0))
   - ``tinn_m`` (float): M value of the TINN computation (right corner of the interpolated triangle at (M, 0))
   - ``tinn`` (float): TINN (baseline width of the interpolated triangle) [ms]

**Parameter Computation**

The Triangular Index is computed as the ratio between the total number of NNIs and the maximum of the NNI histogram
distribution (D(x)).

.. math::

   Tri = \frac{Number of NNIs}{D(X)}

**Application Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

Use the ``legend`` input parameter do show or hide the legend in the histogram figure.

The ``show`` parameter only has effect if ``plot`` is set to True. If ``plot`` is False, no plot figure will be
generated, therefore, no figure can be shown using the ``show`` parameter.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute Triangular Index parameters
   results = td.triangular_index(nn)

   # Get & print the Triangular Index
   print(results['tri_index'])

Geometrical Parameters Function
-------------------------------

.. py:function:: pyhrv.time_domain.geometrical_parameters(nn=None, rpeaks=None, binsize=7.8125, plot=True, show=True, figsize=None, legend=True)

**Function Description:**

Computes all the geometrical parameters based on the NNI histogram (Triangular Index, TINN, N, M) and returns them in a single histogram plot figure.

.. figure:: /_static/geometrical.png
   :align: center
   :scale: 40%

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).
   - ``binsize`` (int, float, optional): Bin size of the histogram bins (default: 7.8125ms).
   - ``plot`` (bool, optional): If True, create the histogram plot figure using ``matplotlib``. If False, the histogram data is computed using ``numpy`` with generating a histogram plot figure (default: True).
   - ``show`` (bool, optional): If True, shows the histogram plot figure (default: True).
   - ``figsize`` (array, optional): 2-element array with the ``matplotlib`` figure size ``figsize``. Format: ``figsize=(width, height)`` (default: will be set to (6, 6) if input is None).
   - ``legend`` (bool, optional): If True, adds legend to the histogram plot figure (default: True).

.. note::

   The ``binsize`` is pre-defined at 7.8125ms and is determined from the minimum suitable sampling frequency for ECG signals of 128Hz as recommended by the `HRV Guidelines <https://www.ahajournals.org/doi/full/10.1161/01.cir.93.5.1043>`_. At this sampling frequency, the temporal resolution of the signal used to derive NNI series is limited at 7.8125ms (= 1/128Hz).

**Returns (ReturnTuple Object):**

The results of this function are returned in a ``biosppy.utils.ReturnTuple`` object (see also :ref:`ref-returntuple`. Use the following keys below (on the left) to index the results.

   - ``nn_histogram`` (matplotlib figure object): Histogram plot figure (only if input parameter ``plot`` is True
   - ``tinn_n`` (float): N value of the TINN computation (left corner of the interpolated triangle at (N, 0))
   - ``tinn_m`` (float): M value of the TINN computation (right corner of the interpolated triangle at (M, 0))
   - ``tinn`` (float): TINN (baseline width of the interpolated triangle) [ms]
   - ``tri_index`` (float): Triangular index [ms]

**Parameter Computation**

See :ref:`ref-tinn` and :ref:`ref-triindex` for detailed information.

**Application Notes**

If both ``nn`` and ``rpeaks`` are provided, the ``nn`` will be chosen over the ``rpeaks`` to avoid additional computational costs.

``nn`` data provided in seconds (s) will automatically converted to milli seconds (ms). See section :ref:`ref-nnformat` for more information.

Use the ``legend`` input parameter do show or hide the legend in the histogram figure.

The ``show`` parameter only has effect if ``plot`` is set to True. If ``plot`` is False, no plot figure will be
generated, therefore, no figure can be shown using the ``show`` parameter.

**Example**

The following example code demonstrates how to use this function and how access the results stored in the ``biosppy.utils.ReturnTuple`` object. This example uses a NNI series from the ``pyhrv/samples/`` folder (see :ref:`ref-samples` for more information).

.. code-block:: python

   # Import packages
   import numpy as np
   import pyhrv.time_domain as td

   # Load sample data
   nn = np.load('./samples/series_1.npy')

   # Compute Triangular Index parameters
   results = td.geometrical_parameters(nn)

   # Get & print the Triangular Index & the TINN value
   print(results['tri_index'])
   print(results['tinn'])


.. _ref-timedomain:


