# !/bin/bash -l
# ## Job Name
# PBS -N cdo_PRECL
# ## Project code
# PBS -A P06010014
# PBS -l walltime=8:30:00
# PBS -q casper
# ## Merge output and error files
# PBS -j oe
# ## Select 1 nodes with 1 CPUs each
# PBS -l select=1:ncpus=1:mem=50GB
# ## Send email on abort, begin and end
# PBS -m abe
# ## Specify mail recipient
# PBS -M dhueholt@rams.colostate.edu
# exec &> logfile.txt

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

### Load modules
module load cdo

EMEM=(
"r1*"
"r2*"
"r3*"
"r4*"
)

IN_PATH_CNTRL="/Users/dhueholt/Documents/GLENS_data/annual_T/ssp245/"
IN_TOKEN_CNTRL="*CESM2-WACCM*"
OUT_PATH_CNTRL="/Users/dhueholt/Documents/GLENS_data/annual_T/ssp245/cdo_test/"
for em in ${EMEM[@]}; do
    IN_TOKEN_CNTRLEM=$IN_TOKEN_CNTRL$em
    sh do_cdo_prep.sh $IN_PATH_CNTRL $IN_TOKEN_CNTRLEM $OUT_PATH_CNTRL
done
#
# IN_PATH_FDBCK="/Users/dhueholt/Documents/GLENS_data/SCIRIS/T/"
# IN_TOKEN_FDBCK="*SSP245*"
# OUT_PATH_FDBCK="/Users/dhueholt/Documents/GLENS_data/SCIRIS/T/cdo_test/"
# for em in ${EMEM[@]}; do
#     IN_TOKEN_FDBCKEM=$IN_TOKEN_FDBCK$em
#     echo $IN_TOKEN_FDBCKEM
#     sh do_cdo_prep.sh $IN_PATH_FDBCK $IN_TOKEN_FDBCKEM $OUT_PATH_FDBCK
# done
