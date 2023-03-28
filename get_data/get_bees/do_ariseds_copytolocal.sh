#######################################
# Copy ARISE-DelayedStart SAI files to local directory.
# Uses "pre-release" directories, will need to be updated
# once moved to final location.
# Globals:
#   None
# Arguments:
#   IN_PATH
#   MOD_TOKEN
#   TIME_TOKEN
#   OUT_PATH
# Written by Daniel Hueholt
# Graduate Research Assistant at Colorado State University
#######################################

export TMPDIR=/glade/scratch/dhueholt/temp
mkdir -p $TMPDIR

IN_TOKEN=$1
MOD_TOKEN=$2
TIME_TOKEN=$3
OUT_PATH=$4

CMN_PATH="/glade/scratch/bdobbins/ARISE/"
CMN_FOLD_STR="b.e21.BW.f09_g17.SSP245-TSMLT-GAUSS-DELAYED-2045."
EMEM=(
"001"
"002"
"003"
"004"
"005"
"006"
"007"
"008"
"009"
"010"
)
# CMN_PATH + CMN_FOLD_STR + EMEM = directory structure for each ens member
PROC="/proc/tseries/"
S="/"

for em in ${EMEM[@]}; do
    FILE_TO_COPY=$CMN_PATH$CMN_FOLD_STR$em$S$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN
    cp $FILE_TO_COPY $OUT_PATH
done
