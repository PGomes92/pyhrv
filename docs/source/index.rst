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
-  ... HRV report generation in .TXT, .CSV and .PDF formats
-  ... and many other useful features to support your HRV research!

Installation
------------

This package can be installed using the `pip` tool:

.. code:: bash

    $ pip install pyhrv

Disclaimer & Context
--------------------
This program is distributed in the hope it will be useful and provided to you "as is", but WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. This program is NOT intended for medical diagnosis. We expressly disclaim any liability whatsoever for any direct, indirect, consequential, incidental or special damages, including, without limitation, lost revenues, lost profits, losses resulting from business interruption or loss of data, regardless of the form of action or legal theory under which the liability may be asserted, even if advised of the possibility of such damages.

This package has initially (up to version 0.3) been developed within the scope of my master thesis "Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)" at the University of Applied Sciences Hamburg, Germany (Faculty Life Sciences, Department of Biomedical Engineering) and PLUX wireless biosignals, S.A., Lisbon, Portugal.
