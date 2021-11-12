#!/bin/bash -l
### Job Name
#PBS -N copy_arise_tsa
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
exec &> logfile_arise_copy.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

OUT_PATH="/glade/scratch/dhueholt/monthly_TSA/"

CMN_PATH="/glade/campaign/cesm/collections/ARISE-SAI-1.5/"
CMN_FOLD_STR="b.e21.BW.f09_g17.SSP245-TSMLT-GAUSS-DEFAULT."
EMEM=(
"001"
"002"
"003"
"004"
"005"
"006"
"007"
"008"
"009"
"010"
)
# CMN_PATH + CMN_FOLD_STR + EMEM = directory structure for each ens member
CMN_SUB_PATH="/lnd/proc/tseries/month_1/"
IN_TOKEN="*.TSA.*"

for em in ${EMEM[@]}; do
    FILE_TO_COPY=$CMN_PATH$CMN_FOLD_STR$em$CMN_SUB_PATH$IN_TOKEN
    cp $FILE_TO_COPY $OUT_PATH
done
