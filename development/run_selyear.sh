#!/bin/bash -l
### Job Name
#PBS -N selyear_defsPeriod_dailySST
### Project code
#PBS -A P06010014
#PBS -l walltime=10:00
#PBS -q casper
### Merge output and error files
#PBS -j oe
### Select 1 nodes with 1 CPUs each
#PBS -l select=1:ncpus=1:mem=60GB
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M dhueholt@rams.colostate.edu
exec &> logfile_selyear.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

module load cdo

cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.001.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.001.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.002.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.002.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.003.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.003.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.004.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.004.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.005.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.005.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.006.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.006.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.007.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.007.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.008.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.008.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.009.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.009.pop.h.nday1.SST.20150101-20191231_RG.nc'
cdo selyear,2015/2019 '/glade/scratch/dhueholt/daily_SST/selname/regrid/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.010.pop.h.nday1.SST.20150101-20650101_RG.nc' '/glade/scratch/dhueholt/daily_SST/selname/regrid/defsPeriod/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.010.pop.h.nday1.SST.20150101-20191231_RG.nc'
