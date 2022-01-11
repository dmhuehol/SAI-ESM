#!/bin/bash -l
### Job Name
#PBS -N copy_acntrl_var
### Project code
#PBS -A P06010014
#PBS -l walltime=10:00
#PBS -q casper
### Merge output and error files
#PBS -j oe
### Select 1 nodes with 1 CPUs each
#PBS -l select=1:ncpus=1:mem=10GB
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M dhueholt@rams.colostate.edu
exec &> logfile_get_hist.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

IN_TOKEN=$1
MOD_TOKEN=$2
TIME_TOKEN=$3
OUT_PATH=$4

### Raw historical
CMN_PATHH="/glade/campaign/collections/cmip/CMIP6/timeseries-cmip6/b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM."
EMEMH=(
"001"
"002"
"003"
)
PROC="/proc/tseries/"
S="/"
# CMN_PATHH + EMEMH + /MOD_TOKEN + PROC + TIME_TOKEN = directory structure for each ens member

for emh in ${EMEMH[@]}; do
    FILE_TO_COPYH=$CMN_PATHH$emh$S$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN
    cp $FILE_TO_COPYH $OUT_PATH
done
