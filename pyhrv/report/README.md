# pyHRV Reports

pyHRV comes with built-in features to create reports from the computed HRV results in .TXT, .CSV. and .PDF format (new in v.0.4). Find below more information about the different output formats.

## .TXT Report
Follow the code example below to generate a .TXT report. An example report can be found here. [Click here for an example TXT report](../files/SampleReport.txt)

```python
    # Imports
    import pyhrv

    # Load 5 minute sample NNI series
    nni = pyhrv.utils.load_sample_nni()

    # Compute HRV parameters (& hide
    results = pyhrv.hrv(nni, show=False)

    # Create HRV TXT report
    pyhrv.report.hrv_report(results, path='./files/', rfile='SampleReport', file_format='txt')
```

## .CSV Report
Follow the code example below to generate a .CSV report. An example report can be found here. [Click here for an example CSV report](../files/SampleReport.csv)

```python
    # Imports
    import pyhrv

    # Load 5 minute sample NNI series
    nni = pyhrv.utils.load_sample_nni()

    # Compute HRV parameters (& hide
    results = pyhrv.hrv(nni, show=False)

    # Create HRV CSV report
    pyhrv.report.hrv_report(results, path='./files/', rfile='SampleReport', file_format='csv')
```

## .PDF Report
#### Requirements & Dependencies
The pyHRV PDF report generator is a LaTex powered feature which allows you generate high-quality reports in PDF format. [Click here for an example PDF report](../files/SampleReport.pdf).

A LaTex distribution is required to be installed on your computer along with the following LaTeX packages in order for this feature to be used:

* amsmath
* array
* colorbl
* datetime
* fancyhdr
* geometry
* graphics
* helvet
* hyperref
* layouts
* titlesec
* tikz
* xcolor

Visit the [The LaTeX Project Website](https://www.latex-project.org/get/) to download a LaTex distribution for your operating system and for more information about how to install the required packages.

#### Creating a PDF report
See the following example demonstrating how to compute the HRV parameters and to generate a PDF report.

```python
    # Imports
    import pyhrv

    # Load 5 minute sample NNI series
    nni = pyhrv.utils.load_sample_nni()

    # Compute HRV parameters (& hide
    results = pyhrv.hrv(nni, show=False)

    # Step 1: Create a PDFReport object and pass the ECG signal, NNI, or R-Peaks series and the results
    report = pyhrv.report.PDFReport(nni=nni, results=results)

    # Step 2: Set general information about the acquisition
    report.set_general_info(subject='Jon Doe',
                            experiment='Sample Report',
                            age=27,
                            gender='male',
                            comment='This is a sample comment in a sample report')

    # Step 3: Create the PDF report
    report.create_report(terminal_output=True)
```

#### Create your own PDF report template
The LaTeX template files for the PDF report are available in th the [pyhrv/report/templates](./templates) folder. Open these files in the LaTeX Editor of your preference to create your own report template. 

Note that the HRV results are stored as new commends in the [parameters.tex](parameters.tex) file which are then used throughout the report template. It is recommended to change only the template while leaving the parameters.tex file as is.

If possible, help pyHRV grow by leaving a reference on your custom reports. For this, use either the citation format suggested on this [repositories' main README](https://github.com/PGomes92/pyhrv) or add a note to the footer with the repositories URL:

```latex
   \lfoot{HRV results \& report generated with pyHRV (v.\version)\\\url{https://github.com/PGomes92/pyhrv}}
```
