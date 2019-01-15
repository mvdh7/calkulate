# Introduction

# Software set-up

## Installation

Installation instructions for Python (steps 1 to 6) and MATLAB (1 to end):

  1. [Download and install Anaconda](https://www.anaconda.com/download/). Choose the Python 3.X version (although the 2.X should work too, if you already have it)

  2. Open the Anaconda Prompt (Windows) or Terminal (Mac/Linux)

  3. Create a new Python 3.6 environment by entering the following:

```
conda create -n calkenv python=3.6 numpy scipy
```

  4. Activate the new environment

  Mac/Linux:

```
source activate calkenv
```

  Windows:

```
activate calkenv
```

  You should now see the environment's name (i.e. `calkenv`) appear in brackets at the start of each line in the Anaconda Prompt/Terminal

  5. Install the Calkulate package into the environment using pip:

```
pip install calkulate
```

  6. You should now be able to use Calkulate in this Python environment. If you wish to also use Calkulate in MATLAB, continue to step 7 onwards

  7. Still within the Anaconda Prompt/Terminal (making sure that the calkenv environment is active), run Python:

```
python
```

  8. Find the location of this environment's Python executable by entering the following 2 lines:

```python
from sys import executable
print(executable)
```

  9. Copy the string that appears. It should look something like (on Windows):

```
C:\Users\username\anaconda\Anaconda3\envs\calkenv\python.exe
```

  This string is the value for the `python_exe` variable that goes into the MATLAB function `calk_initpy()`

  10. Exit python

```python
exit()
```

  11. Download the MATLAB function wrappers (LINK). These are a set a functions that make it easier for you to use some parts of Calkulate within MATLAB, although they are just for convenience -- it's possible to use the entire program without them.

  12. Move the downloaded folder to a sensible location, and add it (plus all subfolders) to your MATLAB search path

  13. Before you can execute the MATLAB functions you must first run the `callk_initpy()` function at least once (per MATLAB session), with the input `callk_initpy()` string obtained in step 9

## Updates

### Python

Update instructions for Python:

  1. Open the Anaconda Prompt

  2. Activate the calkenv environment (installation instructions, step 4)

  3. Upgrade the Calkulate package using pip:

```
pip install calkulate --upgrade --no-cache-dir
```

### MATLAB

Update instructions for MATLAB:

  1. Delete your original Calkulate scripts

  2. Replace them in full with the new versions

## Testing

### MATLAB

First, check that MATLAB and Python are talking to each other properly by entering the following into the MATLAB Command Window (with `python_exe` set to the value of the string you identified using `sys.executable` during the installation process):

```matlab
calk_initpy(python_exe)
pyversion
```

The `executable` that then appears in the Command Window should match the `python_exe` string that you used.

Next, test that you can import Python's `numpy` and `scipy` packages:

```matlab
np = py.importlib.import_module('numpy');
sp = py.importlib.import_module('scipy');
```

If these are working, then new variables `np` and `sp` of type *module* should appear in your MATLAB Workspace. If not, or if you get an error, then there is something wrong with your Python installation.

If all is well, try the same for Calkulate:

```matlab
calk = py.importlib.import_module('calkulate');
```

Finally, test out the MATLAB functions as follows:

```matlab
[Macid,pH,Tk,Msamp,Cacid,S,XT,KX] = calk_Dickson1981;
```

This should import the simulated titration data from Table 1 of Dickson (1981). A plot of the Free scale pH (`pH`) against the acid mass (`Macid`) should appear as follows.

# References

Dickson, A. G. (1981). An exact definition of total alkalinity and a procedure for the estimation of alkalinity and total inorganic carbon from titration data. Deep-Sea Res. Pt A 28, 609–623. <a href="https://doi.org/10.1016/0198-0149(81)90121-7">doi:10.1016/0198-0149(81)90121-7</a>.
