import setuptools
import calkulate as calk

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Calkulate",
    version=calk.__version__,
    author=calk.__author__,
    author_email="m.p.humphreys@icloud.com",
    description="Seawater total alkalinity from titration data",
    url="https://github.com/mvdh7/calkulate",
    packages=setuptools.find_packages(),
    install_requires=["numpy>=1.17", "PyCO2SYS==1.4.0", "scipy>=1.4", "pandas>=1"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
