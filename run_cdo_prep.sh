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
exec &> logfile.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

# Load modules
module load cdo

IN_PATH_CNTRL="/glade/scratch/dhueholt/"
IN_TOKEN_CNTRL="*control.001.*"
OUT_PATH_CNTRL="/glade/scratch/dhueholt/annual_o3/"
sh do_cdo_prep.sh $IN_PATH_CNTRL $IN_TOKEN_CNTRL $OUT_PATH_CNTRL

IN_PATH_FDBCK="/glade/scratch/dhueholt/"
IN_TOKEN_FDBCK="*feedback.001.*"
OUT_PATH_FDBCK="/glade/scratch/dhueholt/annual_o3/"
sh do_cdo_prep.sh $IN_PATH_FDBCK $IN_TOKEN_FDBCK $OUT_PATH_FDBCK
