#!/bin/bash -l
### Job Name
#PBS -N selyear_hist_ocntemp
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

cdo selyear,2010/2014 '/glade/scratch/dhueholt/monthly_OCNTEMP/historical/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.001.pop.h.TEMP.185001-201412.nc' '/glade/scratch/dhueholt/monthly_OCNTEMP/historical/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.001.pop.h.TEMP.201001-201412.nc'
cdo selyear,2010/2014 '/glade/scratch/dhueholt/monthly_OCNTEMP/historical/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.002.pop.h.TEMP.185001-201412.nc' '/glade/scratch/dhueholt/monthly_OCNTEMP/historical/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.002.pop.h.TEMP.201001-201412.nc'
cdo selyear,2010/2014 '/glade/scratch/dhueholt/monthly_OCNTEMP/historical/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.003.pop.h.TEMP.185001-201412.nc' '/glade/scratch/dhueholt/monthly_OCNTEMP/historical/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.003.pop.h.TEMP.201001-201412.nc'
