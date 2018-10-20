# -*- coding: utf-8 -*-
"""
Heart Rate Variability
----------------------

This package provides tools and functions to compute HRV parameters.

Notes
-----
..  This package is part of the master thesis
	"Development of an Open-Source Python Toolbox for Heart Rate Variability (HRV)".
..	This package is a contribution to the open-source biosignal processing toolbox 'BioSppy':
	https://github.com/PIA-Group/BioSPPy
..	See 'references.py' for the full list of references

Author
------
..  Pedro Gomes, Master Student, University of Applied Sciences Hamburg

Thesis Supervisors
------------------
..  Hugo Silva, PhD, Instituto de Telecomunicacoes, PLUX wireless biosignals S.A.
..  Prof. Dr. Petra Margaritoff, University of Applied Sciences Hamburg

Last Update
-----------
13-09-2018

:copyright: (c) 2018 by Pedro Gomes (HAW Hamburg)
:license: BSD 3-clause, see LICENSE for more details.
"""
from __future__ import absolute_import
from pyhrv.__version__ import __version__
from pyhrv.hrv import hrv

# Metadata
__author__ = "Pedro Gomes"
__email__ = "pedro.gomes@haw-hamburg.de"
__maintainer__ = "Pedro Gomes"
__status__ = "Development"
__license__ = "BSD 3-Clause License"
name = "pyhrv"
description = "Python toolbox for Heart Rate Variability."
