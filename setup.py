import setuptools
from calkulate import __version__

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name         = 'Calkulate',
    version      = __version__,
    author       = 'Matthew P. Humphreys',
    author_email = 'm.p.humphreys@cantab.net',
    description  = 'Seawater total alkalinity from titration data',
    url          = 'https://github.com/mvdh7/calkulate',
    packages     = setuptools.find_packages(),
    install_requires = [
        'numpy>=1.15',
        'scipy>=1.1'],
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    classifiers = (
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'),)
