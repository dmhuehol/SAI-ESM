#!/bin/bash -l
### Job Name
#PBS -N copy_ssp245
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
exec &> logfile_ssp245_copy.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

OUT_PATH="/glade/scratch/dhueholt/ssp245/"

CMN_PATH="/glade/collections/cdg/data/CMIP6/ScenarioMIP/NCAR/CESM2-WACCM/ssp245/"
EMEM=(
"r1"
"r2"
"r3"
"r4"
"r5"
)
CMN_FOLD_STR="i1p1f1"
# CMN_PATH + CMN_FOLD_STR + EMEM = directory structure for each ens member
CMN_SUB_PATH="/Amon/tas/gn/latest/"
IN_TOKEN="*.nc"

for em in ${EMEM[@]}; do
    FILE_TO_COPY=$CMN_PATH$em$CMN_FOLD_STR$CMN_SUB_PATH$IN_TOKEN
    cp $FILE_TO_COPY $OUT_PATH
done
