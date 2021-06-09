#!/bin/bash -l
### Job Name
#PBS -N cdo_T
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

EMEM=(
"001.cam.h0.T.*"
"002.cam.h0.T.*"
"003.cam.h0.T.*"
"004.cam.h0.T.*"
"005.cam.h0.T.*"
"006.cam.h0.T.*"
"007.cam.h0.T.*"
"008.cam.h0.T.*"
"009.cam.h0.T.*"
"010.cam.h0.T.*"
"011.cam.h0.T.*"
"012.cam.h0.T.*"
"013.cam.h0.T.*"
"014.cam.h0.T.*"
"015.cam.h0.T.*"
"016.cam.h0.T.*"
"017.cam.h0.T.*"
"018.cam.h0.T.*"
"019.cam.h0.T.*"
"020.cam.h0.T.*"
"021.cam.h0.T.*"
)

IN_PATH_CNTRL="/glade/scratch/dhueholt/"
IN_TOKEN_CNTRL="*control."
OUT_PATH_CNTRL="/glade/scratch/dhueholt/annual_T/"
for em in ${EMEM[@]}; do
    IN_TOKEN_CNTRLEM=$IN_TOKEN_CNTRL$em
    echo IN_TOKEN_CNTRLEM
    sh do_cdo_prep.sh $IN_PATH_CNTRL $IN_TOKEN_CNTRLEM $OUT_PATH_CNTRL
done

IN_PATH_FDBCK="/glade/scratch/dhueholt/"
IN_TOKEN_FDBCK="*feedback."
OUT_PATH_FDBCK="/glade/scratch/dhueholt/annual_T/"
for em in ${EMEM[@]}; do
    IN_TOKEN_FDBCKEM=$IN_TOKEN_FDBCK$em
    echo IN_TOKEN_FDBCKEM
    sh do_cdo_prep.sh $IN_PATH_FDBCK $IN_TOKEN_FDBCKEM $OUT_PATH_FDBCK
done
