.. pyHRV - OpenSource Python Toolbox for Heart Rate Variability documentation master file, created by
   sphinx-quickstart on Wed Oct 17 02:39:51 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/pyhrv.png
    :align: center

``pyHRV`` is a toolbox for Heart Rate Variability (HRV) written in Python.
The toolbox bundles a selection of functions to compute Time Domain, Frequency Domain, and nonlinear HRV parameters,
along with other additional features designed to support your HRV research.

.. toctree::
   :maxdepth: 2
   :numbered:
   :glob:
   :caption: Contents:

   _pages/start
   _pages/api
   _pages/tutorials
   _pages/license

Highlights
----------
This Python package computes...

-  ... fundamental HRV data series (NNI, ∆NNI, HR)
-  ... HRV Time Domain parameters
-  ... HRV Frequency Domain parameters
-  ... nonlinear HRV parameters

... and comes with variety of additional HRV tools, such as...

-  ... ECG and Tachogram plotting features
-  ... export and import of HRV results to .JSON files
-  ... HRV report generation in .TXT or .CSV formats
-  ... and many other useful features to support your HRV research!

Installation
------------

This package can be installed using the `pip` tool:

.. code:: bash

    $ pip install pyhrv

.. note::

   This has been primarily developed for the Python 2.7 programming language. Running the pip command above may cause
   errors when trying to install the package using Python 3.

   In this case, try to install the pyHRV dependencies first:

   .. code:: bash

       $ pip install biosppy
       $ pip install matplotlib
       $ pip install numpy
       $ pip install scipy
       $ pip install nolds
       $ pip install spectrum


Context of this Work
--------------------
This package is part of the master thesis "Development of an Open-Source Python Toolbox for Heart Rate Variability
(HRV)" developed at the University of Applied Sciences Hamburg, Germany, in collaboration with PLUX Wireless
Biosignals, S.A., and the IT - Instituto de Telecomunicações.
