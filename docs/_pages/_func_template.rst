FUNCTION NAME
#############

.. py:function:: pyhrv.MODULE.FUNCTION(nn=None, rpeaks=None)

**Function Description:**

[DESCRIPTION]

**Input Parameters:**
   - ``nn`` (array): NN intervals in (ms) or (s).
   - ``rpeaks`` (array): R-peak times in (ms) or (s).

**Returns (ReturnTuple Object):**
   - ``nn20`` (int): Number of NN interval differences greater 20ms [key: 'nn20']
   - ``PARAM1`` (TYPE): DESCRIPTION [key: 'KEY']

**Parameter Computation**

[DESCRIPTION]

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
   print(results['nn20'])a