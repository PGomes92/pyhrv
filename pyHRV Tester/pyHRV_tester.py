import sys
import os
import numpy as np
import datetime as dt
import pyhrv
from pyhrv.report import hrv_report

PRINTER = True

nni_data = pyhrv.utils.load_sample_nni("long")
rpeaks_data = np.cumsum(nni_data)

tfile = open('pyhrv_test_report_python_%s-%s-%s.txt' % sys.version_info[:3], 'w')

print(pyhrv.__version__)


def output(info, output=''):
	msg = "[%s] %*s %s\n" % (info, 8-len(info), '', output)
	tfile.write(msg)
	if PRINTER:
		print(msg)


def output_success(func_name, data_type):
	output('OK', '%s() passed with %s input data.' % (func_name, data_type))


def output_error(func_name, error, data_type):
	output('ERROR', '%s() encountered the following error with %s input data and did not pass:\n\t\t\t--> %s' % (func_name, data_type, error))


def output_header(header=''):
	output_seperator()
	output('%s' % header.upper())


def output_seperator():
	output('%s' % '#' * 120)


def test_func(f, nni=True):
		input_ = 'nni' if nni else 'rpeaks'
		data_ = 'nni_data' if nni else 'rpeaks_data'
		try:
			# Try to pass the show and mode argument to to suppress PSD plots
			eval(f + '(%s=%s, mode=\'dev\')' % (input_, data_))
		except TypeError as e:
			if 'mode' in str(e):
				try:
					# If functions has now mode feature but 'mode' argument, but a plotting feature
					eval(f + '(%s=%s, show=False)' % (input_, data_))
				except TypeError as a:
					if 'show' in str(a):
						# If functions has now plotting feature try regular function
						eval(f + '(%s=%s)' % (input_, data_))
					else:
						raise TypeError(e)


pyhrv_funcs = pyhrv.utils.load_hrv_keys_json()


def test_module(funcs):
	for func in funcs:
		# Test with NNI as input data
		try:
			test_func(func)
			output_success(func, 'NNI')
		except Exception as e:
			output_error(func, e, 'NNI')

		# test with rpeaks as input data
		try:
			test_func(func, False)
			output_success(func, 'R-Peak')
		except Exception as e:
			output_error(func, e, 'R-Peak')


#######################
# GENERAL INFORMATION #
#######################
output_seperator()
output('Python', 'Version %s' % sys.version)
output('pyHRV', 'Version %s' % pyhrv.__version__)
output('Date', dt.datetime.now())

#####################
# TIME DOMAIN TESTS #
#####################
output_header('Time Domain')
time_funcs = set([x[-1] for x in pyhrv.utils.load_hrv_keys_json().values() if x[0] == 'time'])
test_module(time_funcs)
test_module(["pyhrv.time_domain.time_domain"])

##########################
# FREQUENCY DOMAIN TESTS #
##########################
output_header('Frequency Domain')
frequency_funcs = set([x[-1] for x in pyhrv.utils.load_hrv_keys_json().values() if 'frequency' in x[-1]])
test_module(frequency_funcs)
test_module(["pyhrv.frequency_domain.frequency_domain"])


##########################
# NONLINEAR DOMAIN TESTS #
##########################
output_header('Nonlinear Parameters')
nonlinear_funcs = set([x[-1] for x in pyhrv.utils.load_hrv_keys_json().values() if 'nonlinear' in x[-1]])
test_module(nonlinear_funcs)
test_module(["pyhrv.nonlinear.nonlinear"])

######################
# HRV FUNCTION TESTS #
######################
output_header('HRV Function')
test_module(["pyhrv.hrv"])

###############
# HRV Reports #
###############
output_header('Testing HRV TXT & CSV Reports')
results = pyhrv.hrv(nni_data, show=False)
for ftype in ['txt', 'csv']:
	try:
		hrv_report(results, path=os.getcwd(), rfile='SampleReport', file_format=ftype)
		output_success('pyhrv.report.hrv_report() %s report' % ftype, 'NNI')
	except Exception as e:
		output_error('pyhrv.report.hrv_report() %s report' % ftype, e, 'NNI')


#########
# TOOLS #
#########
output_header('Testing functions of the tools module')
tools1 = ['pyhrv.tools.%s' % x for x in ['tachogram', 'heart_rate', 'heart_rate_heatplot', 'time_varying']]
test_module(tools1)

try:
	pyhrv.tools.radar_chart(nni=nni_data, comparison_nni=nni_data, parameters=['sdnn', 'sdnn'], show=False)
	output_success('pyhrv.tools.radar_chart', 'NNI')
except Exception as e:
	output_error('pyhrv.tools.radar_chart', e, 'NNI')

try:
	pyhrv.tools.radar_chart(rpeaks=nni_data, comparison_rpeaks=nni_data, parameters=['sdnn', 'sdnn'], show=False)
	output_success('pyhrv.tools.radar_chart', 'NNI')
except Exception as e:
	output_error('pyhrv.tools.radar_chart', e, 'NNI')