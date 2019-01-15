import setuptools

with open('README.md','r') as fh:
    long_description = fh.read()

setuptools.setup(
    name         = 'calkulate',
    version      = '2.0.21',
    author       = 'Matthew P. Humphreys',
    author_email = 'm.p.humphreys@cantab.net',
    description  = 'Seawater total alkalinity from titration data',
    url          = 'https://github.com/mvdh7/calkulate',
    packages     = setuptools.find_packages(),
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    classifiers = (
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',),)
