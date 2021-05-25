# GLENS
Repository for code analyzing stratospheric aerosol injection model output. This code is designed to work with output from the [Geoengineering Large ENSemble (GLENS)](https://www.cesm.ucar.edu/projects/community-projects/GLENS/). However, the medium-term intention is that it will apply to output from the [upcoming NCAR model runs performed under the Safe Climate Research Initiative](https://federallabs.org/news/ncar-noaa-lead-efforts-to-understand-risks-and-benefits-of-solar-geoenginering) with minor tuning.

## Code description

### Data preprocessing
`casper_starter`: Contains the generic header and export command required for all jobs run on Casper or Cheyenne.  
`do_cdo_prep`: Merges GLENS monthly netcdf files, apply time shift, and calculates annual mean fields.  
`run_cdo_prep`: Runs `do_cdo_prep` on a set of control and feedback data files. (Submit using `qsub` on Casper/Cheyenne.)  
`run_untar_prep`: Untars a set of control and feedback files and moves to scratch directory. (Submit using `qsub` on Casper/Cheyenne.)

### Analysis and plotting
#### "Bread and butter"
`plot_difference_globe_script`: Plots differences for GLENS output variable on a 4-panel globe  
`plot_single_difference_globe_script`: Plots differences for GLENS variable on a single globe  
`plot_pdf_script`: Plots pdfs (kde, histogram, or step plots)  
`plot_timeseries_script`: Plots timeseries from GLENS output

### Other

#### Helpers
`difference_over_time`: module for calculating values related to changes over time  
`fun_convert_unit`: module for unit conversions  
`plotting_tools`: module for plotting GLENS output  
`process_glens_fun`: module for processing GLENS output (summing levels, etc.)  
`region_library`: Contains lat/lon bounds and other info for a variety of useful regions in easily-callable format.  
`run_test_script`: This is the coding equivalent of scratch paper.

#### Legacy
`plot_decadalpdf_script_libbyformat_ensmn`: legacy kept for method of ensemble means  
`prep_files_script`: uses Python-based CDO to prepare files, kept in case of offline processing

## Sources and Credit
Unless specified otherwise, all code and documentation was written by Daniel Hueholt as a Graduate Research Assistant advised by Profs. [Elizabeth Barnes](https://sites.google.com/rams.colostate.edu/barnesresearchgroup/home) and [James Hurrell](https://sites.google.com/rams.colostate.edu/hurrellgroup/home) at [Colorado State University](https://www.colostate.edu/).
