# matstanlib
**matstanlib** is a MATLAB library for visualization, processing, and analysis of output from Bayesian models fit with [Stan](https://mc-stan.org/).  

**matstanlib** places a heavy emphasis on diagnostic analysis of model output.  the library includes a full set of automated computational checks, consistent with current best practice for Bayesian modeling, as well as a wide variety of diagnostic plots.  

**matstanlib** also supports the use of model output as the basis for inference.  the libarary includes many functions to plot posterior densities and estimates, in addition to analysis functions for parameter estimation and model comparison.  

**matstanlib** is being actively developed and maintained.


## Dependencies
The library is designed to be used in conjunction with Stan and one of the MATLAB interfaces to Stan.  

As such, to use **matstanlib** you first must have a working installation of:
* MATLAB (R2020a or later) + MATLAB's Statistics and Machine Learning Toolbox (aka "stats toolbox")
* the command-line version of Stan: [CmdStan](https://mc-stan.org/users/interfaces/cmdstan.html)
* a MATLAB interface to Stan (such as [MatlabStan](https://github.com/brian-lau/MatlabStan) or [Trinity](https://github.com/joachimvandekerckhove/trinity))

The current test & development environment for **matstanlib** is MATLAB R2021a, CmdStan v2.26.1 and v2.27.0, and a fork of MatlabStan.  


## Installation
Simply clone (or download and unzip) the **matstanlib** repository on your machine, and add the absolute path to the **matstanlib** folder, including all subfolders, to MATLAB's search path:
```
addpath(genpath('/path/to/matstanlib'))
```
To ensure **matstanlib** is always available, add the above line (with a real path, of course) to your `startup.m` file.

## Documentation
**matstanlib** is heavily internally documented.  To see the most recent & complete documentation for a given function, run `help` followed by the name of a function at the MATLAB command line. For example:
```
help extractsamples
```
```
help tracedensity
```
```
help computebfmi
```

The [**matstanlib** wiki](https://github.com/baribault/matstanlib/wiki) offers guidance on how to use **matstanlib**, gives an overview of all **matstanlib** functions, and more.  (The wiki is currently a work in progress!) 

## Getting started
If you are starting from scratch, try running one of the scripts in the `example` folder to see how **matstanlib** functions support a Bayesian workflow in MATLAB.  The `example_correlation.m` script is an ideal place to start! 

If you already have Stan model output in the workspace, start with the functions in the `core` folder:
1. `extractsamples.m` converts your output to a more standard format that is compatibile with all **matstanlib** functions
2. `mcmctable.m` uses the posterior samples to compute summary statistics and convergence diagnostics for each parameter
3. `interpretdiagnostics.m` prints a full diagnostic report and warns when any computational checks fail

To start writing your own Bayesian modeling scripts with Stan & **matstanlib**, check out :skull: `skeleton.m`.
