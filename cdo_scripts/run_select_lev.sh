#!/bin/bash -l
### Job Name
#PBS -N cdo_select_sst
### Project code
#PBS -A P06010014
#PBS -l walltime=30:00
#PBS -q casper
### Merge output and error files
#PBS -j oe
### Select 1 nodes with 1 CPUs each
#PBS -l select=1:ncpus=1:mem=20GB
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M dhueholt@rams.colostate.edu
exec &> logfile.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

### Load modules
module load cdo

IN_PATH="/Users/dhueholt/Documents/GLENS_data/annual_OCNTEMP/"
IN_TOKEN="*.nc"
IN_LEV=500
OUT_PATH="/Users/dhueholt/Documents/GLENS_data/annual_OCNTEMP/cdo_test/"

sh do_select_lev.sh $IN_PATH $IN_TOKEN $IN_LEV $OUT_PATH
