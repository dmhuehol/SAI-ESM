#!/bin/bash -l
### Job Name
#PBS -N cdo_cesmarisetemp_sellevel_mproc
### Project code
#PBS -A P06010014
#PBS -l walltime=25:00
#PBS -q casper
### Merge output and error files
#PBS -j oe
### Select 1 nodes with 5 CPUs each
#PBS -l select=1:ncpus=5:mem=80GB
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M dhueholt@rams.colostate.edu
exec &> logfile_cdo_mproc.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

### Load modules
module load cdo
module load conda/latest
conda activate dh-env

python wrap_mproc_cdo_prep.py