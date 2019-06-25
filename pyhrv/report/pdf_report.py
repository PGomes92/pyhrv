#!/usr/bin/env python -W ignore::FutureWarning
# -*- coding: utf-8 -*-
"""
pyHRV - PDF Report Generator
----------------------------

This modules contains the PDFReport class which can be used to create HRV PDF reports for HRV results computed with
pyHRV. The PDF module uses LaTeX and pdflatex to generate the report. The LaTeX files can be edited to create custom
PDF reports if needed (see 'Notes' below for more information).

Notes
-----
..  The use of this module requires a installed LaTex distribution (e.g. MikTex) and the following LaTeX libraries
    .. amsmath          .. array            .. colorbl          .. datetime         .. fancyhdr         .. geometry
    .. graphicx         .. helvet           .. hyperref         .. layouts          .. titlesec         .. tikz
    .. xcolor
..  All parameter results for the PDF report are stored as new LaTeX commands in the parameter.tex file
..  Feel free to create your own PDF report format by editing the LaTeX templates and adding the HRV results using the
    commands found in the parameters.tex file
..  Supporting pyHRV when creating custom PDF reports using your own LaTeX templates is not mandatory but would be highly
    and greatly appreciated! You can do this by mentioning pyHRV and the link to the GitHub repository to the PDF footer:

    Example Footer Text:    HRV results & report generated with pyHRV (v.X.X.X)
                            https://github.com/PGomes92/pyhrv

Author
------
..  Pedro Gomes, pgomes92@gmail.com

Docs
----
..	You can find the documentation for this module here:
    [LINK]

Last Update
-----------
11-05-2019

:copyright: (c) 2019 by Pedro Gomes
:license: BSD 3-clause, see LICENSE for more details.

"""
# Compatibility
from __future__ import print_function

# Imports
import os
import sys
import shutil
import warnings
import subprocess
import numpy as np
import datetime as dt
from distutils.spawn import find_executable

# BioSPPy imports
from biosppy.signals.ecg import ecg

# Local imports / pyHRV imports
import pyhrv
import pyhrv.time_domain as td
import pyhrv.frequency_domain as fd
import pyhrv.nonlinear as nl
from pyhrv import tools
from pyhrv import utils

class PDFReport:

    def __init__(self, signal=None, nni=None, rpeaks=None, results=None, sampling_rate=None):
        """

        Parameters
        ----------
        signal : array_like
            ECG signal
        nni : array_like
            NN intervals in (ms) or (s)
        rpeaks : array_like
            R-peak times in (ms) or (s)
        results : dict, biosppy.utils.ReturnTuple object
            pyHRV parameter results
        sampling_rate : int, float
            Sampling rate used for the ECG acquisition in (Hz)

        """
        # Check LaTeX requirements
        # self._check_latex_requirements()

        # Todo fix file paths
        self._temp_path = os.path.split(__file__)[0]
        self._figure_path = os.path.join(self._temp_path, 'figures/')
        self._path = os.path.split(__file__)[0]
        self._main = os.path.join(self._path, 'main.tex')
        self._build = os.path.join(self._path, 'build')
        self._file = None

        # Initialize general info
        self._parameter_output = ''
        self._general_info = {'subject': 'n/a', 'experiment': 'n/a', 'age': 'n/a', 'gender': 'n/a', 'comment': 'n/a'}
        self._general_info_latex = {}
        self.set_general_info()
        self.figsizes = (12, 4)

        # Reset sections that should be shown in the report
        self.sections = self._reset_sections()
        self.results = results if results is not None else {}
        self.signal = signal
        self.sampling_rate = sampling_rate if sampling_rate is not None else 1000.

        # Get the list of available pyHRV parameters, keys, and other information (see ./files/hrv_keys.json)
        self.hrv_keys = utils.get_pyhrv_keys()

        # Get the NNI series
        if self.signal is not None:
            rpeaks = ecg(self.signal, self.sampling_rate, show=False)[2]
            self.nni = utils.check_input(None, rpeaks)
        else:
            # Check input data series
            if nni is not None or rpeaks is not None:
                self.nni = utils.check_input(nni, rpeaks)

        # Clear all the data and files from the working directory
        self.clear()


    def _check_latex_requirements(self):
        # Check if LaTeX is installed
        if not find_executable('latex'):
            raise IOError('The PDF Report generation is based on LaTex, which does not seem to be installed on your'
                          'machine. Please ensure that LaTeX is installed in order to use this feature.')

        # TODO add package path for Windows and Linux
        # Set default MiKTeX package path
        if sys.platform == 'darwin':
            package_path = '/Library/Application Support/MiKTeX/texmfs/install/tpm/packages/'
        elif sys.platform == 'linux':
            package_path = ''
        elif sys.platform == 'win32':
            package_path = ''

        # Check if the default package path exists
        if not os.path.exists(package_path):
            raise IOError("The MiKTeX package path '%s' could not be found. Please install MiKTeX or correct the path"
                          "where your LaTeX packages are stored." % package_path)

        # Get the list of required LaTeX packages from the settings.tex file
        required_packages = []
        with open(os.path.join(os.path.split(__file__)[0], 'settings.tex'), 'r') as f:
            for line in f.readlines():
                if 'usepackage' in line:
                    required_packages.append(line.split('{')[-1].replace('}', '').replace('\n', ''))

        # Check if the LaTeX packages are available
        missing_packages = []
        for package in required_packages:
            if not os.path.exists(os.path.join(package_path, '%s.tpm' % package)):
                missing_packages.append(package)

        if missing_packages != []:
            raise IOError("Could not find the following LaTeX packages installed on your machine. Please make sure to"
                          "install them using the MiKTeX console before proceeding: %s" % missing_packages)

    def set_general_info(self, subject=None, experiment=None, age=None, gender=None, comment=None):
        """Set information that will be listed in the General Info table on page 1 of the PDF report

        Parameters
        ----------
        subject : str
            Name of the subject.
        experiment : str
            Information about the experiment or experiment identifier.
        age : int
            Age of the subject.
        gender : str
            Gender of the subject.
        comment : str
            Comments.

        """
        self._general_info_latex['version'] = self._new_command('version', str(pyhrv.__version__), 'r')

        for var in self._general_info.keys():
            self._general_info[var] = eval(var)
            if eval(var) is None:
                self._general_info_latex[var] = self._new_command(var, 'n/a', 'r')
            else:
                self._general_info_latex[var] = self._new_command(var, eval(var), 'r')

    def create_report(self, fname=None, output_path=None, terminal_output=False):
        """Creates the PDF report by generating all figures, setting the LaTeX file containing the parameter results and
        and runs pdflatex to create the PDF report.

        Parameters
        ----------
        fname : str, optional
            Report file name. If None (i.e. not provided), the default file name will created with time stamp:
            pyhrv_report_YYYY-MM-DD_HH-MM-SS.pdf
        output_path : str
            Output path where the PDF report will be saved (default: 'build' folder in the pdfreport submodule of the
            pyHRV package
        terminal_output : boolean
            If True, outputs pdflatex output to Python terminal, else surpresses output

        """
        # Create figures for the PDF report
        self._prepare_figures()
        self._create_parameter_listing()

        # Get all necessary paths
        self._path = os.path.split(__file__)[0]
        self._main = os.path.join(self._path, 'main.tex')
        self._build = os.path.join(self._path, 'build')
        if output_path is not None:
            if os.path.isdir(output_path):
                self._path = output_path
            else:
                raise ValueError("Invalid output path '%s'. Please check the provided output path." % output_path)

        # Set the PDF file name
        if fname is None:
            fname = 'pyhrv_report%s.pdf' % dt.datetime.now().strftime('_%Y-%m-%d_%H-%M-%S')
        else:
            error = True
            while error:
                try:
                    fname = tools._check_fname(os.path.join(self._build, '%s.pdf' % fname), 'pdf')[1][1]
                    error = False
                except:
                    error = True

        # Copy the main.tex template to a new tex file with the PDF file name
        self._file = os.path.join(self._path, '%s.tex' % fname.split('.')[0])
        with open(self._main, 'r') as r:
            with open(self._file, 'w') as w:
                for line in r.readlines():
                    w.write(str(line))

        # Run PDF LaTeX twice to get all the references
        if terminal_output:
            print("Running PDFLaTeX (Step 1 of 2)")
            subprocess.call(['pdflatex', self._file], cwd=self._path, shell=False)
            print("Running PDFLaTeX (Step 2 of 2)")
            subprocess.call(['pdflatex', self._file], cwd=self._path, shell=False)
        else:
            subprocess.call(['pdflatex', self._file], cwd=self._path, shell=False, stdout=open(os.devnull, 'wb'))
            subprocess.call(['pdflatex', self._file], cwd=self._path, shell=False, stdout=open(os.devnull, 'wb'))

        # Move the PDF file and clean the working directory
        self._file = self._file.split('\\')[-1].split('.')[0]
        shutil.move(os.path.join(self._path, '%s.pdf' % self._file), os.path.join(os.path.join(self._build, '%s.pdf' % self._file)))
        self.clear()

    def clear(self):
        """Clears working files created during the creation of the PDF report in the following folders:
            > ./pdffactory/
            > ./pdffactory/figures folder

        ...leaving only the following files available in ./pdffactory/
            > __init__.py
            > main.tex
            > parameters.tex    <-- also resets this file (all values to 'n/a')
            > pyhrv.png
            > pdf_report.py
            > settings.tex

        """
        # Delete all available figures
        for f in os.listdir(self._figure_path):
            os.remove(self._figure_path + f)

        # Remove all working files (e.g. main.aux)
        for file in os.listdir(self._path):
            if file not in ['__init__.py', 'main.tex', 'parameters.tex', 'pyhrv.png', 'pdf_report.py', 'settings.tex']:
                file = os.path.join(self._path, file)
                if os.path.isfile(file):
                    os.remove(file)

        # Reset parameters and the visibility of the section
        with open(os.path.join(self._temp_path, 'parameters.tex'), 'w') as f:
            f.write('')

        # Reset the visibility of all sections
        self.sections = self._reset_sections()

    def _reset_sections(self):
        """Sets the visibility of all sections of the PDF report to False (hidden)"""
        self.sections = {
            'ecgsignal': False,
            'tachogram': False,
            'hrheatplot': False,
            'tdomain': False,
            'fdomain': False,
            'welch': False,
            'ar': False,
            'lomb': False,
            'ndomain': False,
            'poincare': False,
            'dfa': False,
            'sampen': False
        }
        return self.sections

    def _create_parameter_listing(self, reset=False):
        """Create .tex file containing all the parameter values as LaTex variables/commands.

        Parameters
        ----------
        reset : True
            Sets all parameter values to 'n/a' (default: False)

        """
        # Get available keys & replace keys to work with latex
        exceptions = {
            'sd1': 'sdone',
            'sd2': 'sdtwo',
            'nn50': 'nnfifty',
            'pnn50': 'pnnfifty',
            'nn20': 'nntwenty',
            'pnn20': 'pnntwenty',
            'dfa_alpha1': 'dfaalphaone',
            'dfa_alpha2': 'dfaalphatwo',
            'dfa_alpha1_beats': 'dfaalphaonebeats',
            'dfa_alpha2_beats': 'dfaalphatwobeats'
        }

        # Check if unknown keys are provided
        unknown_parameter_keys = []
        for key in self.results.keys():
            if key not in self.hrv_keys and key not in ['dfa_plot', 'ar_resampling_frequency', 'ar_interpolation']:
                unknown_parameter_keys.append(key)

        # Throw warning to visualize all the unknown parameters
        if unknown_parameter_keys:
            warnings.warn("\nUnknown HRV parameter or parameter keys provided which will not be included "
                          "in the PDF report: %s" % unknown_parameter_keys, stacklevel=2)

        self._parameter_output += "% General Information\n"
        if 'comment' in results.keys():
            self._general_info_latex['comment'] = results['comment']

        for key in self._general_info_latex.keys():
            self._parameter_output += self._general_info_latex[key]

        # Prepare parameter commands
        self._parameter_output += '\n% Report parameter values'
        fbands = {}
        for p in self.hrv_keys.keys():
            # Get a LaTeX compatible parameter key
            _p = p if p not in exceptions.keys() else exceptions[p]

            # Sort the HRV parameter results
            if p in self.results.keys() and not reset:
                self._set_section(p)
                if 'fft_' not in p and 'ar_' not in p and 'lomb_' not in p:
                    if 'alpha' and 'beats' in p:
                        self._new_command(_p, '%i-%i beats' % (self.results[p][0], self.results[p][-1]))
                    else:
                        self._new_command(_p, self.results[p])
                else:
                    _p = p.replace('fft_', '').replace('ar_', '').replace('lomb_', '')
                    if 'plot' not in p:
                        if 'bands' in p:
                            for band in ['vlf', 'ulf', 'lf', 'hf']:
                                _low = self.results[p][band][0] if self.results[p][band] is not None else 0
                                _high = self.results[p][band][1] if self.results[p][band] is not None else 0
                                fbands[band] = (_low, _high)
                        else:
                            method = ''
                            if 'fft_' in p:
                                method = 'fft'
                            elif 'ar_' in p:
                                method = 'ar'
                            elif 'lomb_' in p:
                                method = 'lomb'

                            if isinstance(self.results[p], tuple):
                                vals = self.results[p]
                                if len(vals) == 4:
                                    for i, band in enumerate(['ulf', 'vlf', 'lf', 'hf']):
                                        self._new_command('%s%s%s' % (method, band, _p), vals[i])
                                elif len(vals) == 3:
                                    self._new_command('%sulf%s' % (method, _p), 'n/a')
                                    for i, band in enumerate(['vlf', 'lf', 'hf']):
                                        self._new_command('%s%s%s' % (method, band, _p), vals[i])
                                elif len(vals) == 2 and 'norm' in _p:
                                    for i, band in enumerate(['lf', 'hf']):
                                        self._new_command('%s%s%s' % (method, band, _p), vals[i])
                            else:
                                self._new_command('%s%s' % (method, _p), self.results[p])
            else:
                self._new_command(_p, 'n/a')

        for band in fbands.keys():
            self._new_command('%slow' % band, fbands[band][0])
            self._new_command('%shigh' % band, fbands[band][1])

        self._parameter_output += '\n% Report section visibility'
        for sec in self.sections.keys():
            self._parameter_output += "\n\\%s%s" % (sec, str(self.sections[sec]).lower())

        # Reset parameters and the visibility of the section
        with open(os.path.join(self._temp_path, 'parameters.tex'), 'w') as f:
            f.write(self._parameter_output)

    def _prepare_figures(self):
        """Re-creates plot figures based on the results for the PDF report"""
        plots = {}

        # Tachogram
        try:
            fig = tools.tachogram(nni=self.nni,
                                  show=False,
                                  figsize=self.figsizes,
                                  interval='complete')['tachogram_plot']
            plots['tachogram'] = fig
            self._set_section('tachogram')
        except Exception as e:
            self._set_section('tachogram', False)
            warnings.warn("\nAn error occurred while trying to create the Tachogram figure for "
                          "the PDF report: \n%s'" % str(e), stacklevel=2)

        # ECG signal plot
        try:
            if self.signal is not None:
                plots['ecg_plot'] = tools.plot_ecg(signal=self.signal,
                                                   show=False,
                                                   interval='complete',
                                                   figsize=self.figsizes)['ecg_plot']
                self._set_section('ecg_plot')
        except Exception as e:
            self._set_section('ecg_plot', False)
            warnings.warn("\nAn error occurred while trying to create the ECG plot figure for "
                          "the PDF report: \n%s'" % str(e), stacklevel=2)

        # Histogram
        try:
            plots['histogram'] = td.triangular_index(nni=self.nni,
                                                     show=False,
                                                     legend=False,
                                                     figsize=self.figsizes)['tri_histogram']
            self._set_section('histogram')
        except Exception as e:
            self._set_section('histogram', False)
            warnings.warn("\nAn error occurred while trying to create the NNI histogram figure for "
                          "the PDF report: \n%s'" % str(e), stacklevel=2)

        # HR Heat Plot
        if isinstance(self._general_info['age'], int) and isinstance(self._general_info['gender'], str):
            try:
                plots['hr_heatplot'] = utils.heart_rate_heatplot(nni=self.nni,
                                                                 show=False,
                                                                 age=self._general_info['age'],
                                                                 gender=self._general_info['gender'],
                                                                 figsize=self.figsizes)['hr_heatplot']
                self._set_section('hrheatplot')
            except Exception as e:
                self._set_section('hrheatplot', False)
                warnings.warn("\nAn error occurred while trying to create the HR heatplotfor "
                              "the PDF report: \n%s'" % str(e), stacklevel=2)

        # Welch's Plot
        if 'fft_plot' in self.results.keys():
            try:
                plots['fft_plot'] = fd.welch_psd(nni=self.nni,
                                                 fbands=self.results['fft_bands'],
                                                 window=self.results['fft_window'],
                                                 show=False,
                                                 show_param=False,
                                                 figsize=self.figsizes)['fft_plot']
                self._set_section('fft_plot')
            except Exception as e:
                self._set_section('fft_plot', False)
                warnings.warn("\nAn error occurred while trying to create the FFT/Welch figure for "
                              "the PDF report: \n%s'" % str(e), stacklevel=2)

        # Welch's Plot
        if 'ar_plot' in self.results.keys():
            try:
                plots['ar_plot'] = fd.ar_psd(nni=self.nni,
                                             fbands=self.results['ar_bands'],
                                             show=False,
                                             show_param=False,
                                             figsize=self.figsizes)['ar_plot']
                self._set_section('ar_plot')
            except Exception as e:
                self._set_section('ar_plot', False)
                warnings.warn("\nAn error occurred while trying to create the AR PSD figure for "
                              "the PDF report: \n%s'" % str(e), stacklevel=2)

        # Welch's Plot
        if 'lomb_plot' in self.results.keys():
            try:
                plots['lomb_plot'] = fd.lomb_psd(nni=self.nni,
                                                 fbands=self.results['lomb_bands'],
                                                 ma_size=self.results['lomb_ma'],
                                                 show=False,
                                                 show_param=False,
                                                 figsize=self.figsizes)['lomb_plot']
                self._set_section('lomb_plot')
            except Exception as e:
                self._set_section('lomb_plot', False)
                warnings.warn("\nAn error occurred while trying to create the AR PSD figure for "
                              "the PDF report: \n%s'" % str(e), stacklevel=2)

        # Poincare
        if 'poincare_plot' in self.results.keys():
            try:
                plots['poincare_plot'] = nl.poincare(nni=self.nni, show=False, legend=False)['poincare_plot']
                self._set_section('poincare_plot')
            except Exception as e:
                self._set_section('poincare_plot', False)
                warnings.warn("\nAn error occurred while trying to create the Poincare plot figure for "
                              "the PDF report: \n%s'" % str(e), stacklevel=2)

        # DFA
        if 'dfa_plot' in self.results.keys():
            try:
                fig = nl.dfa(nn=self.nni, show=False, legend=False)['dfa_plot']
                plots['dfa_plot'] = fig
                self._set_section('dfa')
            except Exception as e:
                self._set_section('dfa', False)
                warnings.warn("\nAn error occurred while trying to create the DFA plot figure for "
                              "the PDF report: %s'" % str(e), stacklevel=2)

        # Save all plot figures
        for f in plots.keys():
            plots[f].savefig("%s%s.png" % (self._figure_path, f.replace('_', '')), dpi=300, bbox_inches='tight')

    def _set_section(self, key, state=True):
        """Activates visibility of the sections in the PDF report depending on the available parameters
        (e.g. if time domain parameters available in 'results' = show time domain table, else hide time domain table
        in report)

        Parameters
        ----------
        key : str
            pyHRV parameter key available in the hrv_keys.json file
        state : boolean
            Value of the section with True = visible and False = hidden in the PDF report (default: True)

        """
        if key in self.hrv_keys.keys():
            # Check the available parameters and set the visibility of sections and tables...
            if self.hrv_keys[key][0] == 'time':
                self.sections['tdomain'] = state
            elif self.hrv_keys[key][0] == 'frequency_fft':
                self.sections['fdomain'] = state
                self.sections['welch'] = state
            elif self.hrv_keys[key][0] == 'frequency_ar':
                self.sections['fdomain'] = state
                self.sections['ar'] = state
            elif self.hrv_keys[key][0] == 'frequency_lomb':
                self.sections['fdomain'] = state
                self.sections['lomb'] = state
            elif key in ['sd1', 'sd2', 'sd_ratio', 'ellipse_area']:
                self.sections['ndomain'] = state
                self.sections['poincare'] = state
            elif 'dfa' in key:
                self.sections['ndomain'] = state
                self.sections['dfa'] = state
            elif key == 'sampen':
                self.sections['ndomain'] = state
                self.sections['sampen'] = state

        # ...and of the plot figures
        if key == 'tachogram':
            self.sections['tachogram'] = state
        if key == 'ecg_plot':
            self.sections['ecgsignal'] = state
        if key == 'histogram':
            self.sections['histogram'] = state
        if key == 'fft_plot':
            self.sections['fdomain'] = state
            self.sections['welch'] = state
        if key == 'ar_plot':
            self.sections['fdomain'] = state
            self.sections['ar'] = state
        if key == 'lomb_plot':
            self.sections['fdomain'] = state
            self.sections['lomb'] = state
        if key == 'poincare_plot':
            self.sections['ndomain'] = state
            self.sections['poincare'] = state
        if key == 'dfa_plot':
            self.sections['dfa'] = state
        if key == 'hrheatplot':
            self.sections['hrheatplot'] = state

    def _new_command(self, var, val, mode='w'):
        """Returns a new commands in LaTeX compatible \newcommand strings

        Parameters
        ----------
        var : string
            Name of the new command
        val : str, int, float
            Value of the new command
        mode : str, char
            If 'w', writes the resulting newcommand string to self._parameter_output string, else returns the
            resulting newcommand string

        Raises
        ------
        ValueError
            If an unknown mode other than 'r' or 'w' is provided

        """
        if mode == 'w':
            if isinstance(val, str):
                self._parameter_output += ("\\newcommand{\\%s}{%s}\n" % (var, val)).replace('_', '')
            elif isinstance(val, float):
                self._parameter_output += ("\\newcommand{\\%s}{%.3f}\n" % (var, val)).replace('_', '')
            elif isinstance(val, int):
                self._parameter_output += ("\\newcommand{\\%s}{%i}\n" % (var, val)).replace('_', '')
            else:
                self._parameter_output += ("\\newcommand{\\%s}{%s}\n" % (var, 'n/a')).replace('_', '')
        elif mode == 'r':
            if isinstance(val, str):
                return ("\\newcommand{\\%s}{%s}\n" % (var, val)).replace('_', '')
            elif isinstance(val, float):
                return ("\\newcommand{\\%s}{%.3f}\n" % (var, val)).replace('_', '')
            elif isinstance(val, int):
                return ("\\newcommand{\\%s}{%i}\n" % (var, val)).replace('_', '')
            else:
                return ("\\newcommand{\\%s}{%s}\n" % (var, 'n/a')).replace('_', '')
        else:
            raise ValueError("Unknown mode '%s'. Set 'w' for writing into the parameters.tex file or 'r' to return the"
                             "string containing the new LaTeX command.")


if __name__ == '__main__':
    # Load 2 minute Sample Data
    ecg_raw = np.loadtxt('SampleECG_Long.txt')[:1000*60*2, -1]

    # Get rpeaks and NNI series
    rpeaks = ecg(signal=ecg_raw, show=False)[2]
    nni = tools.nn_intervals(rpeaks)
    #
    # # Compute a series of HRV results
    # results = tools.join_tuples(td.time_domain(nni, show=False),
    #                             nl.nonlinear(nni, show=False),
    #                             fd.frequency_domain(nni, show=False))

    # Step 1: Create a PDFReport object and pass the ECG signal, NNI, or R-Peaks series and the results
    report = PDFReport(signal=ecg_raw, nni=nni, results={})

    # Step 2: Set general information about the acquisition
    report.set_general_info(subject='Jon Doe',
                            experiment='Sample Report',
                            age=27,
                            gender='male',
                            comment='This is a sample comment in a sample report')

    # Step 3: Create the PDF report
    report.create_report(terminal_output=True)
