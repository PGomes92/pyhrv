pyHRV Function Levels
=====================

The HRV feature extraction functions of the pyHRV toolbox have been implemented and categorized into three levels, which are intended to facilitate the usage and increase the usability of the toolbox according to the needs of the user or programming background. This multilevel-architecture is illustrated in the Figure below.

.. figure:: /_static/function_levels.png
   :align: center

   Multilevel architecture of the pyHRV toolbox.

**Level 1 - HRV Level:**

This level consists of a single function that allows you to compute the full range of HRV parameters using only a single line of code by calling the `hrv()` function found in the `hrv.py`. This function calls all underlying HRV parameter functions and returns the bundled results in a `biosppy.utils.ReturnTuple()` object. Custom settings or parameters for the computation of specific parameters can be passed to the hrv() function if needed using the `kwargs` input dictionaries. Otherwise, the default values will be used.

.. seealso::

  :ref:`ref-hrvfunc`


**Level 2 - Domain Level:**

This level of module/domain functions intended for users that want to compute all the parameters of a specific domain. Using the example of the Time Domain, the `time domain()` function calls all the parameter function of this domain and returns the results bundled in a `biosppy.utils.ReturnTuple()` object. As in the previous level, user-specific settings and parameters can be passed to the domain functions using the available `kwargs` dictionaries. The module level function can be found in the respective domain modules.

.. seealso::

   * :ref:`ref-timedomain` (time_domain.py)
   * :ref:`ref-frequencydomain` (frequency_domain.py)
   * :ref:`ref-nonlinear` (nonlinear.py)

**Level 3 - Parameter Level:**

This level parameter-specific functions for the case that only want to compute single, specific parameters
(individually) (e.g. `sdnn()` returns the SDNN parameter). This allows you to select only the parameters or features
required for your specific application. User-specific settings and parameters can be made directly when overloading
the function.