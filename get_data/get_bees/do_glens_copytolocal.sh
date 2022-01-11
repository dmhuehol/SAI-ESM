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
exec &> logfile_get_glens.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

IN_TOKEN=$1
MOD_TOKEN=$2
TIME_TOKEN=$3
OUT_PATH=$4

CNTRL_PATH="/glade/campaign/cesm/collections/GLENS/Control/"
FDBCK_PATH="/glade/campaign/cesm/collections/GLENS/Feedback/"
PROC="/proc/tseries/"
# GLENS data is divided by Feedback/Control, stored with all ensemble members in
# a single directory
NC="/*.nc"

CNTRL_TO_COPY=$CNTRL_PATH$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN$NC
FDBCK_TO_COPY=$FDBCK_PATH$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN$NC
cp $CNTRL_TO_COPY $OUT_PATH
cp $FDBCK_TO_COPY $OUT_PATH
