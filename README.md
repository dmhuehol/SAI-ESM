# SAI-ESM for Hueholt et al. "Climate speeds help frame relative ecological risk in future climate change and stratospheric aerosol injection scenarios"
Code for analyzing model output from stratospheric aerosol injection (SAI) experiments in Earth system models, particularly the Community Earth System model (CESM). In particular, this code is designed to calculate and visualize the climate speed of 2m temperature as a metric for high-level perturbations to ecology.  

The version of this code available in `reproducibility-climate-speeds-help-frame` specifically accompanies:  
**Hueholt, D.M.**, E.A. Barnes, J.W. Hurrell, A.L. Morrison. "Climate speeds help frame relative ecological risk in future climate change and stratospheric aerosol injection scenarios", submitted Oct. 2023.  

I spun this code off of the main branch of [SAI-ESM](https://github.com/dmhuehol/SAI-ESM) written for [Hueholt et al. 2023 "Assessing Outcomes in Stratospheric Aerosol Injection Scenarios Shortly After Deployment"](https://doi.org/10.1029/2023EF003488), expecting this would be an effective way to reuse a code base I had already spent significant time developing. In reality, the purposes of these projects quickly diverged, and this fundamental choice resulted in irretrievably snarled code. This was a valuable learning experience--but it does, unfortunately, mean that reproducing the figures from "Climate speeds help frame ecological risk..." is not as smooth as for "Assessing Outcomes..."!

## Replicating Hueholt et al.
The method to make each figure is documented in the docstring and comments of each script beginning with `wrap_`. Note that most figures generate without features such as axes, annotations, or titles, as I prefer to add these manually in Keynote for flexibility.

## Get data
Data will be permanently archived at the Open Science Foundation before publication.

## Brief description of code within package
All files written in Python except for the README itself.
* `CustomExceptions`: Custom exceptions
* `fun_calc_var`: Functions for calculations more complicated than unit conversions, but not so complicated they need a new file
* `fun_convert_unit`: Functions for simple in-line unit conversions and calculations
* `fun_derive_data`: Functions to derive and save data from a base dataset to a new netCDF file
* `fun_ens_plot`: Functions to plot data as timeseries with ensemble visualizations
* `fun_plot_slice_globe`: Functions to plot time-slice maps (as for climate speed)
* `fun_plot_tools`: Plotting functions
* `fun_process_data`: Functions for operations used across many scripts like extracting metadata from filenames, opening files in, managing area inputs and realizations, etc.
* `fun_robustness`: Functions to calculate robustness of trends in data (see Hueholt et al. 2023 "Assessing Outcomes..." Supplemental for description of robustness)
* `fun_special_plot`: Special plotting functions (e.g., Fig. 3, plotting colorbars alone)
* `README`: This README file
* `region_library`: Contains regions that can be called from plotting functions
* `wrap_areaexposed_script`: Wrap functions to plot area exposed to a given value of climate speed (Supplementary Fig. 2)
* `wrap_ensplots_script`: Wrap plotting functions to make timeseries with ensemble visualizations (Supplementary Fig. 1)
* `wrap_plot_slice_globe`: Wrap plotting functions to make time-slice maps (Fig. 1, Fig. 2b-d, Supplementary Fig. 3-6)
* `wrap_rangeplot_script`: Wrap plotting functions to plot range of ensemble variability (Fig. 2a, Supplementary Fig. 7)
* `wrap_wrae_script`: Plot warming rate vs. area exposed to threshold value (Fig. 3, Supplementary Fig. 8)

## Questions and answers
### What's with all the variable names referencing "control" and "feedback"?
"Control" refers to the runs following a climate change trajectory WITHOUT SAI ("no-SAI" in Hueholt et al.), while "Feedback" refers to runs where SAI is also deployed in the model ("SAI" in Hueholt et al.). "Control" and "Feedback" were terms that come out of the development of GLENS and ARISE. Choosing better variable names would certainly be part of a hypothetical rewrite of this code base!

### What do the "fun_", "run_", and "wrap_" in the filenames mean?
Generally, code in this repository follows the naming schema:
*    "fun_" contains individual functions
*    "wrap_" wraps these functions to be run from the command line

## Sources and Credit
Unless specified otherwise, all code and documentation was written by [Daniel Hueholt](https://www.hueholt.earth) as a Graduate Research Assistant advised by Prof. [Elizabeth Barnes](https://barnes.atmos.colostate.edu/) and Prof. [James Hurrell](https://sites.google.com/rams.colostate.edu/hurrellgroup/home) at [Colorado State University](https://www.colostate.edu/).  

This work was supported by the National Science Foundation (NSF) Graduate Research Fellowship (Grant 006784) (DMH) and the Defense Advanced Research Projects Agency (DARPA, Grant HR00112290071) (JWH, EAB, ALM). The views, opinions, and/or findings expressed are those of the authors and should not be interpreted as representing the official views or policies of the Department of Defense or the U.S. government.

* Code in this project is licensed under the Open Software License 3.0, included with this repository as LICENSE.txt
* Figures and text associated with this project are licensed under Creative Commons Attribution Share Alike 4.0 International
