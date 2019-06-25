# HRV SETUP SCRIPT

import setuptools
from pyhrv import __author__, __version__, __email__, name, description
# import matplotlib
# matplotlib.use('TkAgg')
with open("README.md", "r") as fh:
	long_description = fh.read()

# Create setup
setuptools.setup(
	name=name,
	version=__version__,
	author=__author__,
	author_email=__email__,
	description=description,
	long_description=long_description,
	long_description_content_type="text/markdown",
	python_requires='>=2.7',
	url="https://github.com/PGomes92/pyhrv",
	keywords=['Heart Rate Variability', 'HRV'],
	setup_requires=[
		'numpy',
		'scipy',
		'biosppy',
		'matplotlib',
		'nolds',
		'spectrum',
	],
	install_requires=[
		'biosppy',
		'matplotlib',
		'numpy',
		'scipy',
		'nolds',
		'spectrum',
	],

	packages=setuptools.find_packages(),
	package_data={
		"pyhrv": ["files/*", "README.md", "references.txt", "files/quickstart/*"]},
	include_package_data=True,
	classifiers=[
		'Intended Audience :: Developers',
		'Intended Audience :: Education',
		'Intended Audience :: Science/Research',
		'Natural Language :: English',
		'License :: OSI Approved :: BSD License',
		'Programming Language :: Python',
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3.4',
		'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7',
		'Operating System :: OS Independent',
	],
)
