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
	python_requires='>=3.6',
	url="https://github.com/PGomes92/pyhrv",
	keywords=['Heart Rate Variability', 'HRV'],
	install_requires=[
		'biosppy>=2.0.0',
		'matplotlib',
		'numpy',
		'scipy',
		'nolds',
		'spectrum>=0.8.1',
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
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7',
		'Programming Language :: Python :: 3.8',
		'Programming Language :: Python :: 3.9',
		'Programming Language :: Python :: 3.10',
		'Programming Language :: Python :: 3.11',
		'Programming Language :: Python :: 3.12',
		'Programming Language :: Python :: 3.13',
		'Operating System :: OS Independent',
	],
)
