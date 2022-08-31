#!/bin/bash -l
### Job Name
#PBS -N makeDefsFile
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
exec &> logfile_make_defs.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

### Load modules
module load cdo
module load conda/latest
conda activate dh-env

python wrap_derive_mhw_definitionFile.py
