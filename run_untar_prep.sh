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

tar -C /glade/scratch/dhueholt/ -xvf /glade/campaign/cesm/collections/GLENS/Control/atm/proc/tseries/monthly/O3/*2090*.tar
tar -C /glade/scratch/dhueholt/ -xvf /glade/campaign/cesm/collections/GLENS/Feedback/atm/proc/tseries/monthly/O3/*2090*.tar

# python prep_files_script.py # wait on that
