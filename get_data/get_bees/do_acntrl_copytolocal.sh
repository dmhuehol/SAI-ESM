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
exec &> logfile_get_arisecntrl.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

IN_TOKEN=$1
MOD_TOKEN=$2
TIME_TOKEN=$3
OUT_PATH=$4

### Raw futures
CMN_PATHRF="/glade/campaign/collections/cmip/CMIP6/timeseries-cmip6/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM."
EMEMRF=(
"001"
"002"
"003"
"004"
"005"
)
PROC="/proc/tseries/"
S="/"
# CMN_PATHRF + EMEMRF + /MOD_TOKEN + PROC + TIME_TOKEN = directory structure for each ens member

for emrf in ${EMEMRF[@]}; do
    FILE_TO_COPYRF=$CMN_PATHRF$emrf$S$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN
    cp $FILE_TO_COPYRF $OUT_PATH
done

### Raw ARISE dedicated futures
# CMN_PATHAF="/glade/campaign/cesm/development/wawg/WACCM6-TSMLT-SSP245/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM."
# EMEMAF=(
# "006"
# "007"
# "008"
# "009"
# "010"
# )
# # CMN_PATHAF + EMEMAF + MOD_TOKEN + PROC + TIME_TOKEN = directory structure for each ens member
#
# for emaf in ${EMEMAF[@]}; do
#     FILE_TO_COPYAF=$CMN_PATHAF$emaf$S$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN
#     echo $FILE_TO_COPYAF
#     # cp $FILE_TO_COPYAF $OUT_PATH
# done
