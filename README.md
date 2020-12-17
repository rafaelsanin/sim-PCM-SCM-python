# PCMSCM toolbox (Python version)

Python version of Gregory Plett's PCMSCM toolbox battery Parallel-cell modules (PCM) and Series-cell modules (SCM) charge/discharge simulator. The original code is written in Matlab which is available in the PCMSCM toolbox at
[mocha-java.uccs.edu/BMS2/](http://mocha-java.uccs.edu/BMS2/).

## PCM simulator

The Parallel-cell modules (PCM) simulator is the `simPCM.py` file, it loads `A123modeldyn.pickle`, a 1RC Enhanced-self Correcting (ESC) model generated with the Gavin Wiggins Python version of Gregory Plett's enhanced self-correcting (ESC) model, available at [github.com/batterysim/esctoolbox-python](https://github.com/batterysim/esctoolbox-python). The `simPCM.py` simulator also imports the ModelDyn model object from `models.py` and the `OCVfromSOCtemp.py` function as modules.

See the comments in each file for more information.

## SCM model

Pending

## Installation

Requires Python 3.6, Matplotlib, NumPy, and Pandas. The preferred method to
install Python 3 and associated packages is with the Anaconda or Miniconda
distribtion available at
[continuum.io/downloads](https://www.continuum.io/downloads).

## Usage

Clone or download the files to your local machine. Start iPython and run the `simPCM.py` file.
