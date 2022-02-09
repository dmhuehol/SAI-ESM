#!/bin/bash -l
### Job Name
#PBS -N copy_ssp245
### Project code
#PBS -A P06010014
#PBS -l walltime=30:00
#PBS -q casper
### Merge output and error files
#PBS -j oe
### Select 1 nodes with 1 CPUs each
#PBS -l select=1:ncpus=1:mem=10GB
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M dhueholt@rams.colostate.edu
exec &> logfile_ssp245_copy.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

### Directory where output is stored
OUT_PATH="/glade/scratch/dhueholt/ssp245/"

### Raw futures
CMN_PATHRF="/glade/campaign/collections/cmip/CMIP6/timeseries-cmip6/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM."
EMEMRF=(
"001"
"002"
"003"
"004"
"005"
)
CMN_SUB_PATHRF="/atm/proc/tseries/month_1/"
# CMN_PATHRF + EMEMRF + CMN_SUB_PATHRF = directory structure for each ens member
IN_TOKENRF="*.PRECT.*" #* .VARNAME.*

for emrf in ${EMEMRF[@]}; do
    FILE_TO_COPYRF=$CMN_PATHRF$emrf$CMN_SUB_PATHRF$IN_TOKENRF
    cp $FILE_TO_COPYRF $OUT_PATH
done

### Raw historical
CMN_PATHH="/glade/campaign/collections/cmip/CMIP6/timeseries-cmip6/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM."
EMEMH=(
"001"
"002"
"003"
)
CMN_SUB_PATHH="/atm/proc/tseries/month_1/"
# CMN_PATHH + EMEMH + CMN_SUB_PATHH = directory structure for each ens member
IN_TOKENH="*.TREFHT.*" #* .VARNAME.*

for emh in ${EMEMH[@]}; do
    FILE_TO_COPYH=$CMN_PATHH$emh$CMN_SUB_PATHH$IN_TOKENH
    cp $FILE_TO_COPYH $OUT_PATH
done

### CMIP6-processed futures
# CMN_PATHF="/glade/collections/cdg/data/CMIP6/ScenarioMIP/NCAR/CESM2-WACCM/ssp245/"
# EMEMF=(
# "r1"
# "r2"
# "r3"
# "r4"
# "r5"
# )
# CMN_FOLD_STRF="i1p1f1"
# # CMN_PATHF + CMN_FOLD_STRF + EMEMF = directory structure for each ens member
# CMN_SUB_PATHF="/Amon/tas/gn/latest/"
# IN_TOKENF="*.nc"
#
# for emf in ${EMEMF[@]}; do
#     FILE_TO_COPYF=$CMN_PATHF$emf$CMN_FOLD_STRF$CMN_SUB_PATHF$IN_TOKENF
#     cp $FILE_TO_COPYF $OUT_PATH
# done
