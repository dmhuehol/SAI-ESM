#!/bin/bash -l
### Job Name
#PBS -N derive_bmhw
### Project code
#PBS -A P06010014
#PBS -l walltime=45:00
#PBS -q casper
### Merge output and error files
#PBS -j oe
### Select 1 nodes with 1 CPUs each
#PBS -l select=1:ncpus=30:mem=180GB
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M dhueholt@rams.colostate.edu
exec &> logfile_derive_bmhw.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

# Load modules
module load python
module load ncarenv
ncar_pylib ncar_pylib_dhueholt

python wrap_derive_binarymhw_script.py
