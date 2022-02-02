#######################################
# Copy ARISE control files to local directory.
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

### Input variables
IN_TOKEN=$1
MOD_TOKEN=$2
TIME_TOKEN=$3
OUT_PATH=$4

### Raw futures
CMN_PATHRF="/glade/campaign/collections/cmip/CMIP6/timeseries-cmip6/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM."
EMEMRF=(
"001"
"002"
"003"
"004"
"005"
)
PROC="/proc/tseries/"
S="/"
# CMN_PATHRF + EMEMRF + /MOD_TOKEN + PROC + TIME_TOKEN = directory structure for each ens member

for emrf in ${EMEMRF[@]}; do
    FILE_TO_COPYRF=$CMN_PATHRF$emrf$S$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN
    cp $FILE_TO_COPYRF $OUT_PATH
done

### Raw ARISE dedicated futures
CMN_PATHAF="/glade/campaign/cesm/collections/CESM2-WACCM-SSP245/b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM."
EMEMAF=(
"006"
"007"
"008"
"009"
"010"
)
# CMN_PATHAF + EMEMAF + MOD_TOKEN + PROC + TIME_TOKEN = directory structure for each ens member

for emaf in ${EMEMAF[@]}; do
    FILE_TO_COPYAF=$CMN_PATHAF$emaf$S$MOD_TOKEN$PROC$TIME_TOKEN$IN_TOKEN
    # echo $FILE_TO_COPYAF
    cp $FILE_TO_COPYAF $OUT_PATH
done
