#!/bin/tcsh
### Job Name
PBS -N test_cdo_python
### Project code
PBS -A P06010014
PBS -l walltime=01:00:00
PBS -q regular
### Merge output and error files
PBS -j oe
PBS -k eod
### Select 2 nodes with 36 CPUs each for a total of 72 MPI processes
PBS -l select=2:ncpus=36:mpiprocs=36
### Send email on abort, begin and end
PBS -m abe
### Specify mail recipient
PBS -M dhueholt@rams.colostate.edu

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

tar -C /glade/scratch/dhueholt/ -xvf /glade/campaign/cesm/collections/GLENS/Control/atm/proc/tseries/monthly/O3/*2090*.tar
tar -C /glade/scratch/dhueholt/ -xvf /glade/campaign/cesm/collections/GLENS/Feedback/atm/proc/tseries/monthly/O3/*2090*.tar

# python prep_files_script.py # wait on that
