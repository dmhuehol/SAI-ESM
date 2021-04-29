#!/bin/bash -l
### Job Name
#PBS -N untar_2020_O3
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

TAR_PATH_CNTRL="/glade/campaign/cesm/collections/GLENS/Control/atm/proc/tseries/monthly/O3/"
TAR_PATH_FDBCK="/glade/campaign/cesm/collections/GLENS/Feedback/atm/proc/tseries/monthly/O3/"
OUT_PATH="/glade/scratch/dhueholt/"

for CNTRLF in "$TAR_PATH_CNTRL"*021*.tar; do
    tar -C $OUT_PATH -xvf "$CNTRLF"
done

for FDBCKF in "$TAR_PATH_FDBCK"*021*.tar; do
    tar -C $OUT_PATH -xvf "$FDBCKF"
done
