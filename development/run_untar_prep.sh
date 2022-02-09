#!/bin/bash -l
### Job Name
#PBS -N untar_glens_tsa
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

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

TAR_PATH_CNTRL="/glade/campaign/cesm/collections/GLENS/Control/lnd/proc/tseries/monthly/TSA/"
TAR_PATH_FDBCK="/glade/campaign/cesm/collections/GLENS/Feedback/lnd/proc/tseries/monthly/TSA/"
OUT_PATH="/glade/scratch/dhueholt/monthly_TSA/"

for CNTRLF in "$TAR_PATH_CNTRL"*.tar; do
    tar -C $OUT_PATH -xvf "$CNTRLF"
done

for FDBCKF in "$TAR_PATH_FDBCK"*.tar; do
    tar -C $OUT_PATH -xvf "$FDBCKF"
done
