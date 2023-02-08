# SAI-ESM for Hueholt et al. 2023 "Assessing Outcomes in Stratospheric Aerosol Injection Scenarios Shortly After Deployment"
Code for analyzing model output from stratospheric aerosol injection (SAI) experiments in Earth system models, particularly the Community Earth System model (CESM). The version of this code available in `main` specifically accompanies:  
**Hueholt, D.M.**, E.A. Barnes, J.W. Hurrell, J.H. Richter, & L. Sun. "Assessing Outcomes in Stratospheric Aerosol Injection
Scenarios Shortly After Deployment" submitted to *Earth's Future*, preprint link to be provided.

For this project, we utilized the following experiments:
1. [Geoengineering Large ENSemble (GLENS)](https://www.cesm.ucar.edu/community-projects/GLENS/), consisting of 21 RCP8.5 Control runs and 21 G1.2(8.5) Feedback runs where SAI is used to target 1.2 Celsius above preindustrial against the RCP8.5 forcing. Citation: [Tilmes et al. 2018](https://doi.org/10.1175/BAMS-D-17-0267.1)
2. [Assessing Responses and Impacts of Solar climate intervention on the Earth system with stratospheric aerosol injection (ARISE-SAI)](https://www.cesm.ucar.edu/community-projects/arise-sai), consisting of 10 G1.5(2-4.5) runs where SAI is used to target 1.5 Celsius above preindustrial with a SSP2-4.5 scenario. Citation: [Richter et al. 2022](https://doi.org/10.5194/gmd-15-8221-2022)
3. CESM2-WACCM Historical, which is used here primarily to extend the ARISE control to 2010 as in GLENS. Citation: [Danabasoglu et al. 2020](https://doi.org/10.1029/2019MS001916)

While this code is targeted at these specific SAI-related model runs, it can (and in the future, will) be modified to work with other experiments and Earth system models.

## Table of Contents
* [Replicating Hueholt et al. 2023](#replicating-hueholt-et-al-2023)  
* [Making other plots](#making-other-plots)  
    * [Difference globe example](#difference-globe-example)  
    * [Timeseries example](#timeseries-example)
* [Get data](#get-data)
    * [Raw from NCAR](#raw-from-ncar)
    * [Pre-processed data](#pre-processed-data)
* [Process data](#process-data)
    * [Climate Data Operators (CDO)](#climate-data-operators-cdo)
    * [Python](#python)
* [Plot data](#plot-data)
* [Requirements](#requirements)
* [Brief description of code within package](#brief-description-of-code-within-package)
* [Questions and answers](#questions-and-answers)
    * [What's with all the variable names referencing "control" and "feedback"?](#whats-with-all-the-variable-names-referencing-control-and-feedback)
    * [What do the "fun_", "run_", and "wrap_" in the filenames mean?](#what-do-the-fun_-run_-and-wrap_-in-the-filenames-mean)
* [Sources and credit](#sources-and-credit)

## Replicating Hueholt et al. 2023
`wrap_paperplots_basicplots_script` generates all difference globes figures, e.g., Figure 1, 3-8, S2. `wrap_paperplots_ensplots_script` yields all timeseries, e.g., Figure 2, S3, S4.

## Making other plots
To make figures that aren't specifically in the paper, use `wrap_basicplots_script` for difference globes and `wrap_ensplots_script` for timeseries.

### Difference globe example
![Difference globe](images/hueholtetal-f03.png)
`wrap_basicplots_script` generates these figures without the title or colorbar; in Hueholt et al. 2023, these are added manually using Keynote. This can easily be modified by the user by editing the functions in `fun_basic_plot`.

### Timeseries example
![Timeseries](images/hueholt-ts-example.png)
The code generates these figures without the title or annotations, which are added manually using Keynote. This behavior can be changed to add titles automatically by modifying the inputs in `wrap_ensplots_script` as stated in the docstring for this file. Similarly, the appearance can also be controlled manually by editing the functions in `fun_ens_plot`.

## Get data
### Raw from NCAR
Scripts in the `get_data` folder are designed to retrieve raw files stored within the glade file system on NCAR's Cheyenne and Casper machines. This must be run from Cheyenne/Casper, so you'll need an account on these machines to obtain data this way!
### Pre-processed data
Pre-processed data from Hueholt et al. 2023 which is ready for plotting is archived at the [Open Science Foundation](https://osf.io/5a2zf/). This repository also contains an [Earth science datasheet](https://github.com/dmhuehol/Datasheets-for-Earth-Science-Datasets) to document the data used in this study.

## Process data
### Climate Data Operators (CDO)
In this repository, [CDO](https://code.mpimet.mpg.de/projects/cdo) is used for many fundamental processing tasks (e.g., annual mean from monthly data, selecting levels or variables) for computational efficiency. Scripts in the `cdo_mproc` folder wrap CDO functions through Python.

### Python
Various more complex processing tasks are carried out through Python scripts.
*   `wrap_ocean_script`: Remaps ocean data from the Parallel Ocean Model (POP) B-grid to lat/lon coordinates, using methods borrowed and adapted from code by [Emily Gordon](https://sites.google.com/view/emilygordon) and [Zachary Labe](https://zacklabe.com/).
*   `wrap_derive_data_script`: Derives and saves data from a base dataset to new files, used for derivations that require multiple datasets or time-consuming calculations (such as Climdex extremes)
*   `wrap_derive_mhw_definitionFile`: Defines MHW baseline for reference period at a given location.
*   `fun_convert_unit`: Contains simple in-line conversions and calculations that can be run directly from within the plotting code (e.g., converting Kelvin to degrees Celsius)

## Plot data
Code to plot data consists of `wrap_basicplots_script` for difference globes and `wrap_ensplots_script` for timeseries. As described above, `wrap_paperplots_basicplots_script` and `wrap_paperplots_ensplots_script` provide one-click replication of all figures from Hueholt et al. 2023, which this repository accompanies.

## Requirements
Required Python packages and versions available on `pip` are listed in the `requirements.txt` file. Additionally, the [marineHeatWaves](https://github.com/ecjoliver/marineHeatWaves) package by [Eric Oliver](https://github.com/ecjoliver) is required to run code related to marine heatwaves.

## Brief description of code within package
All code written in Python unless specified otherwise.
* `CustomExceptions`: Custom exceptions, written partly as a coding exercise
* `fun_basic_plot`: Contains the difference globe plotting functions
* `fun_convert_unit`: Functions for simple in-line unit conversions and calculations
* `fun_derive_data`: Functions to derive and save data from a base dataset to a new netCDF file
* `fun_ens_plot`: Functions to plot data as timeseries with ensemble visualizations (spaghetti, spread)
* `fun_plot_tools`: Plotting functions wrapped by `wrap_basicplots_script` and `wrap_paperplots_basicplots_script`
* `fun_process_data`: Functions for operations used across many scripts like extracting metadata from filenames, opening files in, managing area inputs and realizations, etc.
* `fun_regrid_pop`: Functions to remap ocean from POP B-grid to lat/lon using methods from [Emily Gordon](https://sites.google.com/view/emilygordon) and [Zachary Labe](https://zacklabe.com/)
* `fun_robustness`: Functions to calculate robustness of trends in data (see Hueholt et al. 2023 Supplemental for description of robustness)
* `fun_special_plot`: Special plotting functions for specific uses (e.g., plotting colorbars alone)
* `README`: This README file
* `region_library`: Contains regions that can be called from plotting functions
* `requirements.txt`: Required Python packages available on `pip`
* `run_derive_data_script`: Shell script to run wrap_derive_data_script on NCAR Casper
* `run_ocean_script`: Shell script to run wrap_ocean_script on NCAR Casper
* `wrap_basicplots_script`: Wrap difference globe plotting functions
* `wrap_derive_data_script`: Wrap functions to derive data
* `wrap_derive_mhw_definitionFile`: Wrap function to define MHW baseline for reference period at a given location
* `wrap_ensplots_script`: Wrap plotting functions to make timeseries with ensemble visualizations
* `wrap_ocean_script`: Wrap functions to remap ocean data from POP B-grid to lat/lon coordinates
* `wrap_paperplots_basicplots_script`: Script to instantly replicate all difference globe figures (e.g., Figure 1, 3-8, S2) from Hueholt et al. 2023
* `wrap_paperplots_ensplots_script`: Script to instantly replicate all timeseries figures  (e.g., Figure 2, S3, S4) from Hueholt et al. 2023
* `wrap_plotregions_script`: Make plots of input region(s) on a map
* `wrap_stat_robustness`: Calculate statistical metrics for robustness
* `wrap_test_script`: This is just a file I use as "scratch paper" when coding  

* `cdo_mproc`: Folder of code that wraps CDO functions in Python to process raw ESM data
    * `run_mproc_cdo_prep`: Shell script to run wrap_mproc_cdo_script on NCAR Casper
    * `wrap_mproc_cdo_prep`: Wraps functions which call CDO to carry out various foundational data processing tasks like making annual mean from monthly data
    * `mproc_bees`: This folder contains the individual CDO functions for each task, these are documented within each file  

* `get_data`: Folder of code to obtain data on NCAR Casper  

* `helper_scripts`: Contains a few additional shell scripts for one-off tasks that occasionally need to be run.
    * `do_move_images`: Shell script to move images into folders by region (of limited use to anyone but me)
    * `do_sumtwo_makenew`: Shell script to sum data from two files to make a new variable, such as making total precipitation from convective and large-scale precipitation variables in GLENS
    * `run_selyear`: Shell script to select years from data using CDO
    * `run_sumtwo_makenew`: Shell script to run `do_sumtwo_makenew` on NCAR Casper

* `images`: Folder containing images used in the README file.

## Questions and answers
### What's with all the variable names referencing "control" and "feedback"?
"Control" refers to the runs following a climate change trajectory WITHOUT SAI ("no-SAI" in Hueholt et al. 2023). "Feedback" refers to runs where SAI is also deployed in the model ("SAI" in Hueholt et al. 2023). "Control" and "Feedback" were terms that come out of the development of GLENS and ARISE. I wrote the code in this repository using that standard terminology before settling on the "SAI/no-SAI" terminology we used in the paper.

### What do the "fun_", "run_", and "wrap_" in the filenames mean?
Generally, code in this repository follows the naming schema:
*    "fun_" contains individual functions
*    "wrap_" wraps these functions to be run from the command line
*    "run_" runs scripts from the job scheduler on NCAR's Casper and Cheyenne

## Sources and Credit
Unless specified otherwise, all code and documentation was written by [Daniel Hueholt](https://www.hueholt.earth) as a Graduate Research Assistant advised by Prof. [Elizabeth Barnes](https://barnes.atmos.colostate.edu/) and Prof. [James Hurrell](https://sites.google.com/rams.colostate.edu/hurrellgroup/home) at [Colorado State University](https://www.colostate.edu/).  

This work was supported by the LAD Climate Fund, DARPA (grant no. HR00112290071), and the National Science Foundation Graduate Research Fellowship Program.  

* Code in this project is licensed under the Open Software License 3.0, included with this repository as LICENSE.txt
* Figures and text associated with this project are licensed under Creative Commons Attribution Share Alike 4.0 International
