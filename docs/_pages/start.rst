Getting Started
===============

Installation
############
Use the ``pip`` tool to download ``pyHRV`` from the `Python Package Index (PyPi)
<https://pypi.org>`_. If you are using macOS or a Linux-based operating system, open the *Terminal* and type in the
command below. Windows users may use the *Command Prompt*.

.. code:: bash

   $ pip install pyhrv


``pyHRV`` depends on the following third-party packages which will automatically be installed using the
installation command, if these are not already on your machine:

* `biosppy <https://github.com/PIA-Group/BioSPPy>`_
* `numpy <http://www.numpy.org>`_
* `scipy <http://scipy.org>`_
* `matplotlib <https://matplotlib.org>`_
* `nolds <https://github.com/CSchoel/nolds>`_
* `spectrum <https://github.com/withspectrum/spectrum>`_

.. note::

   Alternatively, it is recommended to install the `Anaconda <https://www.anaconda.com>`_ software which comes with a compatible Python 2.7 distribution and all the necessary (and more) third-party packages for scientific computing.


R-Peak Extraction with BioSPPy
##############################
``BioSppy`` is a toolbox for biosignal processing and comes with built-in ECG processing and R-peak detection
algorithms. These can be used to compute the NNI series upon which the HRV parameters can be computed.

An example of this procedure is demonstrated below using ECG data acquired with the `BITalino (r)evolution <www
.bitalino.com>`_ hardware and the `OpenSignals (r)evolution <http://bitalino.com/en/software>`_ software. The ECG
signals are imported and converted to mV using the `opensignalsreader <https://github
.com/PGomes92/opensignalsreader>`_ package.

.. code:: python

      import biosppy
      import numpy as np
      import pyhrv.tools as tools
      from opensignalsreader import OpenSignalsReader

      # Load sample ECG signal & extract R-peaks using BioSppy
      signal = OpenSignalsReader('./samples/SampleECG.txt').signal('ECG')
      signal, rpeaks = biosppy.signals.ecg.ecg(signal, show=False)[:2]

      # Compute NNI
      nni = tools.nn_intervals(rpeaks)

.. note::

   ``pyHRV`` can of course be used with any ECG - or ECG Lead I like - signal and is not limited to signals acquired with specific devices or software. The instructions above are merely an example.

.. _ref-samples:

Sample NNI Series
#################
`pyHRV` comes with 50 NNI samples which can be used if you have no ECG data or NNI series available. These sample series where extracted from the `MIT-BIH NSRDB Database from physionet.org <https://physionet.org/physiobank/database/nsrdb/>`_ and can be found in the `./samples/ <https://github.com/PGomes92/pyhrv/tree/master/pyhrv/samples>`_ folder.

.. _ref-returntuple:

The biosppy.utils.ReturnTuple Object
########################################
The results of the ``pyHRV`` parameter functions wrapped and returned in ``biosppy.utils.ReturnTuple`` objects. This package-specific class combines the advantages of Python dictionaries (indexing using keywords) and Python tuples (immutable). Parameter values stored
in the ``ReturnTuple`` object can be accessed as follows:

.. code-block:: python

   from biosppy import utils

   # Store sample data in a ReturnTuple object
   args = (500, 600, )
   names = ('parameter1', 'parameter2', )
   results = utils.ReturnTuple(args, names)

   # Get and print 'parameter1'
   print(results['parameter1'])


.. seealso::

   - `BioSPPy API Reference - ReturnTuple <https://biosppy.readthedocs.io/en/stable/biosppy.html#biosppy.utils.ReturnTuple>`_
   - `Note on ReturnTuple objects <https://biosppy.readthedocs.io/en/stable/tutorial.html#a-note-on-return-objects>`_
