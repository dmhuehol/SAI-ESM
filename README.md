# GLENS
Repository for code analyzing stratospheric aerosol injection model output. This code is designed to work with output from the [Geoengineering Large ENSemble (GLENS)](https://www.cesm.ucar.edu/projects/community-projects/GLENS/). However, flexibility is one of my goals for the code in this repository, and the medium-term intention is that it will apply to output from the [upcoming NCAR model runs performed under the Safe Climate Research Initiative](https://federallabs.org/news/ncar-noaa-lead-efforts-to-understand-risks-and-benefits-of-solar-geoenginering) with minimal tuning.

## Code description

### Data preprocessing
`casper_starter`: Contains the generic header and export command required for all jobs run on Casper or Cheyenne.  
`do_cdo_prep`: Merges GLENS monthly netcdf files, apply time shift, and calculates annual mean fields.  
`run_cdo_prep`: Runs `do_cdo_prep` on a set of control and feedback data files. (Submit using `qsub` on Casper/Cheyenne.)

### Analysis and plotting
#### "Bread and butter"

#### Other

#### Helpers
`region_library`: Contains lat/lon bounds and other info for a variety of useful regions in easily-callable format.  
`run_test_script`: This is the coding equivalent of scratch paper.

## Sources and Credit
Unless specified otherwise, all code and documentation was written by Daniel Hueholt as a Graduate Research Assistant under the advisement of Profs. [Elizabeth Barnes](https://sites.google.com/rams.colostate.edu/barnesresearchgroup/home) and [James Hurrell](https://sites.google.com/rams.colostate.edu/hurrellgroup/home) at Colorado State University.
