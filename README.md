# SAI-CESM
Repository for code analyzing stratospheric aerosol injection model (SAI) output utilizing the Community Earth System Model (CESM). This code utilizes output from the following modeling experiments:
1. [Geoengineering Large ENSemble (GLENS)](https://www.cesm.ucar.edu/projects/community-projects/GLENS/), consisting of 21 RCP8.5 Control runs and 21 G1.2(8.5) Feedback runs where SAI is used to target 1.2 Celsius above preindustrial against the RCP8.5 forcing
2. Assessing Responses and Impacts of Solar climate intervention on the Earth system with stratospheric aerosol injection (ARISE-SAI), consisting of 10 G1.5(2-4.5) runs where SAI is used to target 1.5 Celsius above preindustrial with a SSP2-4.5 scenario
3. CMIP6 CESM2-WACCM SSP2-4.5, which functions as the ARISE control with 5 members spanning 2015-2100
4. CESM2-WACCM Historical, which extends the ARISE control back in time to 1850

While this code is targeted at these specific SAI-related model runs, it can be modified to work with other CESM output. The version of this code available in `main` specifically accompanies:
Hueholt et al. 2023 "Assessing Outcomes in Stratospheric Aerosol Injection
Scenarios Shortly After Deployment" submitted to *Earth's Future*.

## Replicating paper results
`wrap_paperplots_basicplots_script` generates all difference globes figures, e.g., Figure 1, 3-8, S2. `wrap_paperplots_ensplots_script` yields all timeseries, e.g., Figure 2, S3, S4.

## Example plots
These figures can be made using `wrap_basicplots_script` and `wrap_ensplots_script`.
### Difference globe example
![Difference globe](images/hueholtetal-f03.png)
The code generates these figures without the title or colorbar; in Hueholt et al. 2023, these are added manually using Keynote.


### Timeseries example
![Timeseries](images/hueholt-ts-example.png)
The code generates these figures without the title or annotations, which are added manually using Keynote.


## Sources and Credit
Unless specified otherwise, all code and documentation was written by [Daniel Hueholt](https://www.hueholt.earth) as a Graduate Research Assistant advised by Prof. [Elizabeth Barnes](https://barnes.atmos.colostate.edu/) and Prof. [James Hurrell](https://sites.google.com/rams.colostate.edu/hurrellgroup/home) at [Colorado State University](https://www.colostate.edu/).
