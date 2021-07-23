--------------------------------------------------------------------------
--- MATSTANLIB v0.6 ------------------------------------------------------
--- a library of helper functions for MATLAB interfaces to Stan ----------
--------------------------------------------------------------------------

this library contains functions for processing, analyzing, and visualizing 
MCMC-generated posterior samples and diagnostic quantities returned by any 
MATLAB interfaces to Stan (mc-stan.org), such as:
        MATLABStan  >>  https://github.com/brian-lau/MatlabStan
        Trinity     >>  https://github.com/joachimvandekerckhove/trinity

the library also contains a few scripts with examples of proper usage.  

please note that in this library of functions, it is generally true that:
    SAMPLES is a structure containing samples for many parameters, where 
        each fieldname is a parameter, and each field's value is a matrix 
        of HMC samples of size: 
        [nIterations nChains <arrayDimensions> <parameterDimensions>]. 
        this is format is intentionally similar to the format used by 
        other interfaces to Stan, such as PyStan. 
    CHAINS is a [nIterations nChains]-size matrix of HMC samples for 
        a scalar parameter or a single instance of a nonscalar parameter.
    a "parameter name" is a string(/char-type array) that matches a 
        monitored parameter (of any dimension) in the Stan model 
        specification (e.g., 'mu', 'alpha_group', 'gamma'). 
    a "parameter instance name" is a string(/char-type array) that denotes 
        a specific instance of a monitored nonscalar parameter.  such a 
        string will take the form of the parameter names + indices, where 
        the index or comma-separated indices are appended in square 
        brackets (e.g., 'mu[3,1,4]', 'alpha_group[2]'). 
        (but please not that a scalar parameter also qualifies as a 
        parameter instance in most contexts in matstanlib.)
    PARAMETER is a parameter name string. 
    PARAMETERS or PARAMETERNAMES is a cell containing one or more parameter 
        name strings. 
    INSTANCES or INSTANCENAMES is a cell of parameter instance name strings
        and/or parameter name strings for scalar parameters. 
    PARAMETERREQUEST may be either a string or cell of strings, where each 
        string is either a parameter name or parameter instance name.  
        whether instances names may be included depends on the function. 
and that:
    wherever I say string, I actually mean character array. 
    (i.e., if string-type is specified, char-type is what is actually 
    intended and syntactically required (@ischar == true).)

I have tried to ensure that this package has good documentation.
typing 'help FUNCTIONNAME' at the command line will give information on how 
to call each function.  scripts are heavily commented. 

however, please also note that this library is a work in progress.  while most
functions and scripts are complete, a few are still being actively worked on, 
and so are not guaranteed to be free of errors (either thrown or silent). 
code may be updated or refactored at any time.  

once the library is in a satisfactory state, I will make an official release.

this code currently requires that MATLAB's Statistics Toolbox be installed.

--------------------------------------------------------------------------

ALL INCLUDED SCRIPTS AND FUNCTIONS ARE BY BETH BARIBAULT (c) 2019, 2020, 2021 ---

--------------------------------------------------------------------------

PLOTTING
---COLORS
    > colorlibrary.m                returns a library of named colors for plotting
    > getcolors.m                   returns RGB values for color or style name
    > makecolormap.m                returns a colormap matrix
---CHAINS
    > tracedensity.m                plots chain traces, w/ optional diagnostic overlays
    > plotautocorr.m                plots the autocorrelation function, by chain
---POSTERIORS
    > plotdensity.m                 plots a univariate density, with optional overlays
    > overlaydensity.m              overlays one or more univariate densities
    > jointdensity.m                plots a bivariate denisty & marginals
    > multidensity.m                plots a matrix of bivariate densities
    > horzdensity.m                 plots smoothed densities for multiple parameters
    > violindensity.m               plots violin densities for multiple parameters
---PARAMETER ESTIMATES
    > plotrecovery.m                plots recovery (true vs. estimated parameter values)
    > plotintervals.m               plots credible intervals for parameter instances
---PREDICTIVES
    > postpredhist.m                ...
---DIAGNOSTICS
    > tracedensity.m                plots chain traces
    > rankplots.m                plot one rank plot per chain
    > plotlp.m                      plots lp__ and other diagnostic quantites
    > parcoordivergent.m            plots joint posterior samples across parameters
    > plotdivergences.m             plots divergent transitions, one rug plot per chain
---MISC
    > hyperpriortester.m            plots simulated priors for given hyperpriors

--------------------------------------------------------------------------

STAN OUTPUT MANIPULATION
    > extractsamples.m              reformat output from a MATLAB interface to Stan
    > getparaminstances.m           generates a list of (indexed) parmeter instances
    > removechain.m                 toss a chain
    > removeiters.m                 toss some iterations
    > transformparameterization.m   apply known reparameterizations post-hoc

--------------------------------------------------------------------------

STAN OUTPUT DIAGNOSTICS
    > computeess.m                  computes effecive sample size
    > computerhat.m                 computes Rhat
    > mcmctable.m                   generates a table of MCMC diagnostics & summary stats
    > interpretdiagnostics.m        interprets MCMC diagnostics and prints a report
    > rhattable.m                   ...

--------------------------------------------------------------------------

STAN BASIC ANALYSIS
    > smoothdensity.m               estimates a posterior density with kernel smoothing
    > getsamplestats.m              collects summary statistics for posterior samples
    > computecredint.m              computes a credible interval
    > computeWAIC.m                 computes WAIC from log_lik samples
    > savagedickey.m                computes a Savage-Dickey Bayes factor

--------------------------------------------------------------------------

STAN BASIC MODELS & SCRIPTS
    > skeleton.m                    bare bones outline of a script
    > example_funnel.m              Neal's funnel; generates divergent transitions
    > example_funnel_bg.m           Betancourt & Giorolami's funnel
    > example_funnel_ncp.m          hierarchical funnel; non-centered parameterization
    > example_correlation.m         Pearson correlation; plots & a Bayes factor
    > example_RL.m                  reinforcement learning; a hierarchical cognitive model
    > example_eightschools.m        ???
    > example_figures.m             *** recreates the figures in [matstanlib paper] ***

--------------------------------------------------------------------------

MISC MATSTANLIB HELPERS
    > writestanfile.m               writes a cell of Stan code strings to a .stan file
    > collecttruevalues.m           generates a structure of true parameter values
    > str2ind.m                     (called by many matstanlib functions)
    > reindexvector.m               converts a vector of values to sequential indices

--------------------------------------------------------------------------

MISC
    > getdeps.m                     prints a dependency report for an .m file

--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------
