# #!/bin/bash -l
# ### Job Name
# #PBS -N get_hi
# ### Project code
# #PBS -A P06010014
# #PBS -l walltime=30:00
# #PBS -q casper
# ### Merge output and error files
# #PBS -j oe
# ### Select 1 nodes with 1 CPUs each
# #PBS -l select=1:ncpus=1:mem=10GB
# ### Send email on abort, begin and end
# #PBS -m abe
# ### Specify mail recipient
# #PBS -M dhueholt@rams.colostate.edu
# # exec &> logfile_getallfiles.txt
#
# export TMPDIR=/glade/scratch/dhueholt/temp
# mkdir -p $TMPDIR

IN_VAR="hi"
VAR_TYPE="ocn" #atm,ice,lnd,ocn,rof
OUT_PATH="/glade/scratch/dhueholt/monthly_hi/" #Files will be put here

### Copy GLENS
GLENS_CMN_PATH="/glade/campaign/cesm/collections/GLENS/"
GLENS_CNTRL="Control"
GLENS_FDBCK="Feedback"
GLENS_CMN_SUBPATH="/proc/tseries/monthly/"
#GLENS_CMN_PATH + GLENS_CNTRL + VAR_TYPE + GLENS_CMN_SUBPATH + VAR_TYPE + *.nc
#GLENS_CMN_PATH + GLENS_FDBCK + VAR_TYPE + GLENS_CMN_SUBPATH + VAR_TYPE + *.nc
