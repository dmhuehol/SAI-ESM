#######################################
# Merge GLENS monthly netcdf files, shift time, and calculate annual mean.
# Output files are named automatically as:
# type_ensnumber_variable_YYYYMM-YYYYMM[first]_..._YYYYMM-YYYYMM[last].nc
# Globals:
#   None
# Arguments:
#   IN_PATH
#   IN_TOKEN
#   OUT_PATH
# Written by Daniel Hueholt
# Graduate Research Assistant at Colorado State University
#######################################

### Input variables
IN_PATH=$1 #Path to data
IN_TOKEN=$2 #Token to match data files, e.g. "*.001.*.nc"
OUT_PATH=$3 #Path to save files

### Script
IN_CARD="$IN_PATH$IN_TOKEN"
PATH_LENGTH=${#IN_PATH}
echo $PATH_LENGTH

RUN_FNAMES=()
RUN_TIMES=()
RUN_ENSNUMS=()
for f in $IN_CARD; do
  ACTIVE_FNAME=${f:$PATH_LENGTH}
  if [[ "$ACTIVE_FNAME" == *"CESM2-WACCM"* ]]; then #CMIP6 format
    ACTIVE_FNAME=${ACTIVE_FNAME//_/.}
    RUN_FNAMES+=( $ACTIVE_FNAME )
    ACTIVE_TIME=$(echo $ACTIVE_FNAME | cut -d'.' -f7)
    RUN_TIMES+=( $ACTIVE_TIME )
    ACTIVE_ENSNUM=$(echo $ACTIVE_FNAME | cut -d'.' -f5)
    RUN_ENSNUMS+=( $ACTIVE_ENSNUM )
    RUN_TYPE=$(echo $ACTIVE_FNAME | cut -d'.' -f4)
    RUN_VAR=$(echo $ACTIVE_FNAME | cut -d'.' -f1)
  else #GLENS or SCIRIS format
    RUN_FNAMES+=( $ACTIVE_FNAME )
    ACTIVE_TIME=$(echo $ACTIVE_FNAME | cut -d'.' -f10)
    RUN_TIMES+=( $ACTIVE_TIME )
    ACTIVE_ENSNUM=$(echo $ACTIVE_FNAME | cut -d'.' -f6)
    RUN_ENSNUMS+=( $ACTIVE_ENSNUM )
    RUN_TYPE=$(echo $ACTIVE_FNAME | cut -d'.' -f5)
    RUN_VAR=$(echo $ACTIVE_FNAME | cut -d'.' -f9)
  fi
done

### Troubleshooting
# echo ${RUN_FNAMES[@]}
# echo ${RUN_TIMES[@]}
# echo $RUN_TYPE
# echo ${RUN_ENSNUMS[@]}
# echo $RUN_VAR

OUT_FNAME="${RUN_TYPE}_${RUN_ENSNUMS[0]}_${RUN_VAR}"
# echo $OUT_FNAME
for t in ${RUN_TIMES[@]}; do
  OUT_FNAME="${OUT_FNAME}_${t}"
done
echo $OUT_FNAME

OUT_MERGE="${OUT_PATH}${OUT_FNAME}_merge.nc"
OUT_SHIFT="${OUT_PATH}${OUT_FNAME}_shift.nc"
OUT_ANNUAL="${OUT_PATH}${OUT_FNAME}_annual.nc"

cdo -L -yearmonmean -shifttime,'-1days' -mergetime ${IN_CARD} ${OUT_ANNUAL}
# cdo -L -selmon,9 -mergetime ${IN_CARD} ${OUT_ANNUAL}
