# SAI-ESM for Hueholt et al. "Speed of environmental change frames relative ecological risk in climate change and climate intervention scenarios"
Code for analyzing model output from stratospheric aerosol injection (SAI) experiments in Earth system models, particularly the Community Earth System model (CESM). In particular, this code is designed to calculate and visualize the climate speed of 2m temperature as a metric for high-level perturbations to ecology. This code includes all required content for the Nature Research Code and Software Submission Checklist.  

The version of this code available in `reproducibility-climate-speeds-help-frame` specifically accompanies:  
**Hueholt, D.M.**, E.A. Barnes, J.W. Hurrell, A.L. Morrison. "Speed of environmental change frames relative ecological risk in climate change and climate intervention scenarios", accepted in principle Mar. 2024.  

I spun this code off of the main branch of [SAI-ESM](https://github.com/dmhuehol/SAI-ESM) written for [Hueholt et al. 2023 "Assessing Outcomes in Stratospheric Aerosol Injection Scenarios Shortly After Deployment"](https://doi.org/10.1029/2023EF003488), expecting this would be an effective way to reuse a code base I had already spent significant time developing. In reality, these projects quickly diverged, and this fundamental choice resulted in irretrievably snarled code. This was a valuable learning experience--but it does, unfortunately, mean reproducing the figures from "Speed of environmental change frames relative..." is not as smooth as for "Assessing Outcomes..."!

## System requirements
* Code written in Python 3.11.7 on a 2020 MacBook Pro running Mac OS Sonoma 14.2.1. The computational complexity is not high and it should function on other platforms.
* At least 1.3GB of free space is required for the accompanying data
* Dependencies are listed in `requirements.txt`
* No non-standard hardware is required

## Install software
**Instructions**: Clone the repository using the method of your choice. 
**Typical install time**: A few seconds for download.

## Get data
Data is permanently archived at the Open Science Framework: [doi.org/10.17605/OSF.IO/Z37ES](https://doi.org/10.17605/OSF.IO/Z37ES). This repository includes a datasheet (Gebru et al. 2021, [Connolly & Hueholt in prep](https://github.com/dmhuehol/Datasheets-for-Earth-Science-Datasets)) for comprehensive documentation of the data. Download this data to the path of your choice, and use this as the ``dataPath`` when running code.

## Replicating Hueholt et al.
The production of each figure within the Main Text and Supplementary Information of Hueholt et al. "Climate speeds help frame..." is documented in the docstring and comments of each script beginning with `wrap_`. Note that most figures generate by default without features such as axes, annotations, or titles, as I add these manually in Keynote for flexibility. Run time for each figure ranges from a few seconds (for ``wrap_plot_slice_globe``) to up to a couple of minutes (for ``wrap_wrae_script``).

## Brief description of code within package
All files written in Python except when specified otherwise.
* `CustomExceptions`: Custom exceptions
* `do_cdo_annualmean.sh`: Shell script example of using CDO to merge CESM files (CDO commands for non-CESM files derived from this manually)
* `fun_calc_var`: Functions for calculations more complicated than unit conversions, but not so complicated they need a new file
* `fun_convert_unit`: Functions for simple in-line unit conversions and calculations
* `fun_derive_data`: Functions to derive and save data from a base dataset to a new netCDF file
* `fun_ens_plot`: Functions to plot data as timeseries with ensemble visualizations
* `fun_plot_slice_globe`: Functions to plot time-slice maps (as for climate speed)
* `fun_plot_tools`: Plotting functions
* `fun_process_data`: Functions for operations used across many scripts like extracting metadata from filenames, opening files in, managing area inputs and realizations, etc.
* `fun_robustness`: Functions to calculate robustness of trends in data (see Hueholt et al. 2023 "Assessing Outcomes..." Supplemental for description of robustness)
* `fun_special_plot`: Special plotting functions (e.g., Fig. 3, plotting colorbars alone)
* `get_periods`: Obtain time periods for the Last Millennium simulation
* `README`: This README file
* `region_library`: Contains regions that can be called from plotting functions
* `requirements.txt`: Python environment description
* `wrap_areaexposed_script`: Wrap functions to plot area exposed to a given value of climate speed (Supplementary Fig. 1)
* `wrap_ensplots_script`: Wrap plotting functions to make timeseries with ensemble visualizations (Fig. 2)
* `wrap_plot_slice_globe`: Wrap plotting functions to make time-slice maps (Fig. 1, Fig. 2b-d, Supplementary Fig. 2-8)
* `wrap_plot_vector_globe`: Wrap plotting functions to make vector maps (Supplementary Fig. 11)
* `wrap_rangeplot_script`: Wrap plotting functions to plot range of ensemble variability (Fig. 2a, Supplementary Fig. 9)
* `wrap_wrae_script`: Plot warming rate vs. area exposed to threshold value (Fig. 4, Supplementary Fig. 10)

## Questions and answers
### What's with all the variable names referencing "control" and "feedback"?
"Control" refers to the runs following a climate change trajectory WITHOUT SAI ("no-SAI" in Hueholt et al.), while "Feedback" refers to runs where SAI is also deployed in the model ("SAI" in Hueholt et al.). "Control" and "Feedback" were terms that come out of the development of GLENS and ARISE. Choosing better variable names would certainly be part of a hypothetical rewrite of this code base!

### What do the "fun_" and "wrap_" filenames mean?
Generally, code in this repository follows the naming schema:
*    "fun_" contains individual functions
*    "wrap_" wraps these functions to be run from the command line

## Sources and Credit
Unless specified otherwise, all code and documentation was written by [Daniel Hueholt](https://www.hueholt.earth) as a Graduate Research Assistant advised by Prof. [Elizabeth Barnes](https://barnes.atmos.colostate.edu/) and Prof. [James Hurrell](https://sites.google.com/rams.colostate.edu/hurrellgroup/home) at [Colorado State University](https://www.colostate.edu/).  

This work was supported by the National Science Foundation (NSF) Graduate Research Fellowship (Grant 006784) (DMH) and the Defense Advanced Research Projects Agency (DARPA, Grant HR00112290071) (JWH, EAB, ALM). The views, opinions, and/or findings expressed are those of the authors and should not be interpreted as representing the official views or policies of the Department of Defense or the U.S. government.

* Code in this project is licensed under the Open Software License 3.0, included with this repository as LICENSE.txt
* Figures and text associated with this project are licensed under Creative Commons Attribution Share Alike 4.0 International
