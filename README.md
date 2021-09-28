# matstanlib
matstanlib is a library of MATLAB functions for diagnostic visualization, processing, and analysis of output from Bayesian models fit with [Stan](https://mc-stan.org/).  

On the diagnostic side, matstanlib offers a full set automated diagnostic checks, consistent with current best practice for Bayesian modeling, as well as a wide variety of diagnostic plots.  

On the analysis side, matstanlib offers even more plots to facilitate visual exploration of model results, in addition to analysis functions for parameter estimation and model comparison.  

matstanlib is being actively developed and maintained.

## Prerequisites
The library is designed to be used in conjunction with Stan and one of the MATLAB interfaces to Stan.  

As such, to use matstanlib you first must have a valid standalone installation of Stan:
* CmdStan: https://mc-stan.org/users/interfaces/cmdstan.html

Please follow the installation instructions specific to your operating system, which are included in CmdStan's User's Guide.  You can verify your installation is working by also following the User's Guide instructions for the `bernoulli` example.  

You should also have a recent version of MATLAB, including "stats toolbox":
* MATLAB R2020a or later
* MATLAB's Statistics and Machine Learning Toolbox

To check if stats toolbox is installed, at the MATLAB command line run: `help ksdenisty`.  If MATLAB returns `ksdensity not found.`, then stats toolbox is not yet installed. 

Finally, you will need one of the following MATLAB interfaces to Stan:
* MatlabStan: https://github.com/brian-lau/MatlabStan 
** Note: MatlabStan is not actively maintained and its code requires editing to work correctly with recent releases of Stan.  
* Trinity: https://github.com/joachimvandekerckhove/trinity
** Note: Trinity is not actively maintained, but works without error recent Stan versions.

Please follow the installation instructions provided with your chosen interface. 

matstanlib is currently being actively tested on MATLAB R2021a, CmdStan v2.27.0, and a fork of MATLABStan.  


## Installation
Simply clone or download the matstanlib repository to your machine, and add the absolute path to the matstanlib folder to MATLAB's search path.  
```
addpath(genpath('/path/to/matstanlib'))
```
To ensure matstanlib is always available, I recommend that you add the above line (with a real path, of course) your `startup.m` file. 


## Getting started
Check out the scripts in the `example` folder to see how matstanlib may be used in a variety of context.  You can also open up `skeleton.m' to get started writing your own modeling scripts.  
