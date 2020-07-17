#!/usr/bin/env bash

# Script to merge HD5 tensors

################### VARIABLES ############################################

DESTINATION=
SOURCES=
NUM_JOBS=96
INTERSECT=""
INPLACE=""
SAMPLE_IDS_START=1000000
SAMPLE_IDS_END=6030000

SCRIPT_NAME=$( echo $0 | sed 's#.*/##g' )


################### HELP TEXT ############################################

usage()
{
    cat <<USAGE_MESSAGE

    This script can be used to create tensors from the UKBB data.

    Usage: ${SCRIPT_NAME}    -t <tensor_path>
                          [-i <id_string>] [-n <num_jobs>]
                          [-s <sample_id_start>] [-e <sample_id_end>]
                          [-x <xml_field_id>] [-m <mri_field_id>] [-c <CONTINUOUS_FIELD_IDS>] [-a <CATEGORICAL_FIELD_IDS>]
                          [-h]

    Example: ./${SCRIPT_NAME} -t /mnt/disks/data/generated/tensors/test/2019-02-05/ -i my_run -n 96 -s 1000000 -e 6030000 -x "20205 6025" -m "20208 20209"

        -d      <path>      (Required) Absolute path to directory to write output tensors to.

        -s      <id>        (Required) List of source directories containing tensors in quotes

        -n      <num>       Number of jobs to run in parallel. Default: 96.

        -i      <str>       Intersect mode (intersect or nothing)

        -p      <str>       Inplace mode (inplace or nothing)

        -b      <id>        Smallest sample ID to start with. Default: 1000000.

        -e      <id>        Largest sample ID to end with. Default: 6030000.

        -h                  Print this help text

USAGE_MESSAGE
}

display_time() {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  (( $D > 0 )) && printf '%d days ' $D
  (( $H > 0 )) && printf '%d hours ' $H
  (( $M > 0 )) && printf '%d minutes ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d seconds\n' $S
}

################### OPTION PARSING #######################################

if [[ $# -eq 0 ]]; then
    echo "ERROR: No arguments were specified." 1>&2
    usage
    exit 1
fi

while getopts ":d:s:n:b:e:hip" opt ; do
    case ${opt} in
        h)
            usage
            exit 1
            ;;
        d)
            DESTINATION=$OPTARG
            ;;
        s)
            SOURCES=$OPTARG
            ;;
        n)
            NUM_JOBS=$OPTARG
            ;;
        b)
            SAMPLE_IDS_START=$OPTARG
            ;;
        e)
            SAMPLE_IDS_END=$OPTARG
            ;;
        i)
            INTERSECT="--intersect"
            ;;
        p)
            INPLACE="--inplace"
            ;;
        :)
            echo "ERROR: Option -${OPTARG} requires an argument." 1>&2
            usage
            exit 1
            ;;
        *)
            echo "ERROR: Invalid option: -${OPTARG}" 1>&2
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

################### SCRIPT BODY ##########################################

# Keep track of and display the elapsed time
START_TIME=$(date +%s)

# Variables used to bin sample IDs so we can tensorize them in parallel
INCREMENT=$(( ( $SAMPLE_IDS_END - $SAMPLE_IDS_START ) / $NUM_JOBS ))
COUNTER=1
MIN_SAMPLE_ID=$SAMPLE_IDS_START
MAX_SAMPLE_ID=$(($MIN_SAMPLE_ID + $INCREMENT))

# Run every parallel job within its own container -- 'tf.sh' handles the Docker launching
while [[ $COUNTER -lt $(( $NUM_JOBS + 1 )) ]]; do
    echo -e "\nLaunching job for sample IDs starting with $MIN_SAMPLE_ID and ending with $MAX_SAMPLE_ID via:"

        cat <<LAUNCH_CMDLINE_MESSAGE
                $HOME/ml/scripts/tf.sh -r -c $HOME/ml/ml4cvd/tensorize/merge_hd5s.py \
		--destination $DESTINATION \
		--sources $SOURCES \
		--min_sample_id $MIN_SAMPLE_ID \
		--max_sample_id $MAX_SAMPLE_ID \
		$INTERSECT $INPLACE &
LAUNCH_CMDLINE_MESSAGE

    $HOME/ml/scripts/tf.sh -r -c $HOME/ml/ml4cvd/tensorize/merge_hd5s.py \
		--destination $DESTINATION \
		--sources $SOURCES \
		--min_sample_id $MIN_SAMPLE_ID \
		--max_sample_id $MAX_SAMPLE_ID \
		$INTERSECT $INPLACE &

    let COUNTER=COUNTER+1
    let MIN_SAMPLE_ID=MIN_SAMPLE_ID+INCREMENT
    let MAX_SAMPLE_ID=MAX_SAMPLE_ID+INCREMENT

    sleep 1s
done

################### DISPLAY TIME #########################################

END_TIME=$(date +%s)
ELAPSED_TIME=$(($END_TIME - $START_TIME))
printf "\nDispatched $((COUNTER - 1)) merge jobs in "
display_time $ELAPSED_TIME
