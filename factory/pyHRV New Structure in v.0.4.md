
    pyhrv                           # Toolbox
    ├── quickstart                  # Quickstart Guide
    ├── files                       # pyHRV support files
    |   ├── hr_heatplot.json        # HR reference normal values for the hr_heatplot() function
    |   ├── hrv_keys.json           # HRV keys to access the parameter results stored in
    |   ├── references.txt          # Publications references on which pyHRV functions are based
    |   ├── SampleECG.txt           # Sample ECG signal (for testing purposes)
    |   ├── SampleExport.json       # Sample HRV results export (for demonstration purposes)
    |   ├── SampleNNISeries.npy     # Sample NNI series (for testing purposes)
    |   └── SampleECG.txt           # Sample ECG signal (for testing purposes)
    |      
    ├── report                      # Subpackage for PDF, TXT, and CSV report generation
    |   ├── build                   # Default path for generated PDF, TXT, and CSV reports 
    |   ├── figure                  # Working directory where figures for the PDF report are temporarily stored
    |   ├── templates               # LaTeX templates for the PDF reports (feel free to customize)
    |   ├── __init__.py             # Report init file
    |   ├── main.tex                # Main LaTeX file
    |   ├── parameters.tex          # File containing all the variables in which the HRV parameters will be stored
    |   ├── pdf_report.py           # PDF report module
    |   ├── pyhrv.png               # pyHRV logo for the PDF header
    |   ├── README.md               # Report README
    |   ├── SampleReport.csv        # Sample report in CSV format
    |   ├── SampleReport.pdf        # Sample report in PDF format
    |   ├── SampleReport.txt        # Sample report in TXT format
    |   └── settings.tex            # Settings file for the LaTeX-based project
    |   
    ├── __init__.py                 # pyHRV init file
    ├── __version__.py              # pyHRV version
    ├── hrv.py                      # HRV function
    ├── frequency_domain_.py        # Frequency Domain functions
    ├── nonlinear.py                # Nonlinear Parameters functions
    ├── time_domain.py              # Time Domain functions
    ├── tools.py                    # HRV tools functions (Tachogram, ECG, radar chart, HR heatplot, and other 
    |                               # Comparison functions
    └── utils.py                    # General data verification, conversion, file handling, etc.


Tools.py
########

- [x] nn_intervals()
- [x] nn_diff()
- [x] plot_ecg()
- [x] tachogram()
- [x] heart_rate()
- [x] heart_rate_heatplot()
- [x] radar_chart()
- [x] time_varying()
- [x] hrv_export()
- [x] hrv_import()

**Restructure - utils.py**
- [x] check_interval()
- [x] check_limits()
- [x] segmentation()
- [x] join_tuples()
- [x] check_fname()
- [x] std()
- [x] time_vector()
- [x] check_input()
- [x] get_pyhrv_keys()
- [x] load_sample_nni()
- [x] pyHRVRadarAxes(PolarAxes)
