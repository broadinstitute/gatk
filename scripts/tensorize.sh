#!/usr/bin/env bash

# Script to generate tensors from UKBB data

################### VARIABLES ############################################

TENSOR_PATH=
ID="tensorization_run"
NUM_JOBS=96
SAMPLE_IDS_START=1000000
SAMPLE_IDS_END=6030000
XML_FIELD=  # exclude ecg data
MRI_FIELD=  # exclude mri data

SCRIPT_NAME=$( echo $0 | sed 's#.*/##g' )

# TODO: Remove the hardcoded field IDs below
CONTINUOUS_FIELD_IDS=" 12679 12680 12681 12682 12683 12684 12685 12702 22605 22631 22641 22651 90012 90013 4207 21021 23098 23099 23100 23101 23102 23104 23105 23106 23107 23108 23109 23110 23111 23112 23113 23114 23115 23116 23117 23118 23119 23120 23121 23122 23123 23124 23125 23126 23127 23128 23129 23130 3160 20015 21001 21002 5084 5085 5086 5087 5096 5097 5098 5099 5116 5117 5118 5119 5132 5133 5134 5135 5156 5157 5158 5159 5160 5161 5162 5163 5215 5274 5254 5255 5256 5257 5262 5263 5264 5265 5078 5079 5199 5201 5204 5206 5208 5211 3083 3084 3085 3086 3143 3144 3146 3147 3148 4100 4101 4103 4104 4105 4106 4119 4120 4122 4123 4124 4125 4138 4139 4140 4141 4142 4143 4144 4145 4146 4147 3062 3063 20019 20021 20022 20007 20009 20011 30000 30010 30020 30030 30040 30050 30060 30070 30080 30090 30100 30110 30120 30130 30140 30150 30160 30170 30180 30190 30200 30210 30220 30230 30240 30250 30260 30270 30280 30290 30300 30500 30510 30520 30530 40008 40007 189 22003 22004 22005 22009 12651 20191 20159 12673 12674 12675 12676 12677 12678 12686 12687 12697 12698 12699 22200 22599 22602 22603 22642 22643 22644 22645 22652 22653 22654 22655 22661 22663 22664 41083 41084 41085 41086 41087 41088 41089 41090 41091 41092 41093 41094 41095 41096 41097 41098 41099 41100 41101 41102 41103 41109 41110 41111 41112 41132 41076 41078 41079 41104 41105 41142 41143 41144 41080 41082 41106 41107 41108 41146 41148 4194 4195 4196 4198 4199 4200 93 94 95 102 4079 4080 6032 6033 6038 6039 5088 5089 5100 5101 5102 5103 5104 5105 5106 5107 5108 5109 5110 5111 5112 5113 5114 5115 5190 5193 5221 5237 5251 5276 5292 5306 6071 6073 5074 5075 5076 5077 5186 5188 5200 5202 5207 5209 3064 21003 20016 4282 4283 4285 396 397 398 399 400 4288 4290 4291 404 20023 1807 1845 1873 1883 2946 3526 3972 3982 5057 2355 3809 2217 4689 4700 5430 5901 5923 5945 2966 2976 3627 3761 3786 3894 3992 4012 4022 4056 4269 4272 4276 4279 1568 1578 1588 1598 1608 4407 4418 4429 4440 4451 4462 5364 1289 1299 1309 1319 1438 1458 1488 1498 1528 3680 864 874 884 894 904 914 1070 1080 1090 1050 1060 1737 2277 1269 1279 2867 2887 2897 2926 3436 3456 3659 2684 2704 2714 2734 2744 2754 2764 2794 2804 2824 3536 3546 3581 3700 3710 3829 3839 3849 3872 3882 2405 129 130 84 87 134 135 137 92 136 30314 30324 30334 30344 30354 30364 30374 30384 30404 30414 40009 34 "
CATEGORICAL_FIELD_IDS=" 42007 42009 42011 42013 42001 42003 42005 12139 12187 12188 12141 12253 12254 12140 12223 12224 20165 20167 20169 20171 20173 20175 20177 20179 20181 20183 20185 20187 20189 20193 12700  22604 22606 22607 22608 22609 22610 22611 22612 22613 22614 22615 22616 22617 22618 22619 22620 22630 22640 22650 22660 4186 4204 6014 6015 6016 6017 6019 6020 6024 6034 5189 5191 5194 5196 6070 6072 5185 5187 5205 5212 6074 6075 4092 4093 4095 4096 3088 3089 3090 3159 4924 5556 5779 4259 4281 4287 4292 4293 4294 20018 1647 1677 1687 1697 1707 1767 1777 1787 1797 1835 3912 3942 4501 20107 20110 20111 2316 4717 4728 5452 5463 5474 5485 5496 5507 5518 5529 5540 2335 3606 3616 3751 2345 2365 2207 2227 5408 5419 5441 5610 5832 5843 5855 5877 5890 5912 5934 6119 2178 2188 2296 2306 2247 2257 3393 4792 4803 4814 4825 4836 2443 2453 2463 2473 2986 3005 4041 2492 2415 2844 2956 3404 3414 3571 3741 3773 3799 4067 4232 4243 4268 4270 4275 4277 4849 1558 1618 1628 2664 3731 3859 20117 1329 1339 1349 1359 1369 1379 1389 1408 1418 1428 1448 1468 1478 1508 1518 1538 1548 2654 924 943 971 981 991 1001 1011 1021 1100 2624 2634 3637 3647 1717 1727 1747 1757 2267 1239 1249 1259 2644 2877 2907 2936 3446 3466 3476 3486 3496 3506 5959 20116 21000 2674 2694 2724 2774 2784 2814 2834 3591 3720 2375 2385 2395 120 20115 3140 3079 31 5181 5182 5183 5324 5325 5326 5327 5328 22000 22001 22006 6147 6148 6150 6151 6152 6153 6154 6155 6177 6179 6149 6159 6144 6162 6164 6157 6158 20001 20002 20003 20004 "

################### HELP TEXT ############################################

usage()
{
    cat <<USAGE_MESSAGE

    This script can be used to create tensors from the UKBB data.

    Usage: ${SCRIPT_NAME}    -t <tensor_path>
                          [-i <id_string>] [-n <num_jobs>]
                          [-s <sample_id_start>] [-e <sample_id_end>]
                          [-x <xml_field_id>] [-m <mri_field_id>]
                          [-h]

    Example: ./${SCRIPT_NAME} -t /mnt/disks/data/generated/tensors/test/2019-02-05/ -i my_run -n 96 -s 1000000 -e 6030000 -x "20205 6025" -m "20208 20209"

        -t      <path>      (Required) Absolute path to directory to write output tensors to.

        -i      <id>        ID string used to distinguish the run. Default: 'tensorization_run'.

        -n      <num>       Number of jobs to run in parallel. Default: 96.

        -s      <id>        Smallest sample ID to start with. Default: 1000000.

        -e      <id>        Largest sample ID to start with. Default: 6030000.

        -x      <ids>       ECG field IDs. Default: None - no ECG data is used.

        -m      <ids>       MRI field IDs. Default: None - no MRI data is used.

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

while getopts ":t:i:n:s:e:x:m:h" opt ; do
    case ${opt} in
        h)
            usage
            exit 1
            ;;
        t)
            TENSOR_PATH=$OPTARG
            ;;
        i)
            ID=$OPTARG
            ;;
        n)
            NUM_JOBS=$OPTARG
            ;;
        s)
            SAMPLE_IDS_START=$OPTARG
            ;;
        e)
            SAMPLE_IDS_END=$OPTARG
            ;;
        x)
            XML_FIELD=$OPTARG
            ;;
        m)
            MRI_FIELD=$OPTARG
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
MAX_SAMPLE_ID=$(( $MIN_SAMPLE_ID + $INCREMENT - 1 ))

# Run every parallel job within its own container -- 'tf.sh' handles the Docker launching
while [[ $COUNTER -lt $(( $NUM_JOBS + 1 )) ]]; do
    echo -e "\nLaunching job for sample IDs starting with $MIN_SAMPLE_ID and ending with $MAX_SAMPLE_ID via:"

        cat <<LAUNCH_CMDLINE_MESSAGE
            $HOME/ml/scripts/tf.sh -nt $HOME/ml/ml4cvd/recipes.py
                --mode tensorize
                --tensors $TENSOR_PATH
                --output_folder $TENSOR_PATH
                --id $ID
                --include_heart_zoom
                --mri_field_id $MRI_FIELD
                --xml_field_id $XML_FIELD
                --categorical_field_ids <too long to list>
                --continuous_field_ids <too long to list>
                --dicoms ./dicoms_$MIN_SAMPLE_ID/
                --min_sample_id $MIN_SAMPLE_ID
                --max_sample_id $MAX_SAMPLE_ID &
LAUNCH_CMDLINE_MESSAGE

    $HOME/ml/scripts/tf.sh -nt $HOME/ml/ml4cvd/recipes.py \
		--mode tensorize \
		--tensors $TENSOR_PATH \
		--output_folder $TENSOR_PATH \
		--id $ID \
		--include_heart_zoom \
		--mri_field_id $MRI_FIELD \
		--xml_field_id $XML_FIELD \
		--categorical_field_ids $CATEGORICAL_FIELD_IDS \
		--continuous_field_ids $CONTINUOUS_FIELD_IDS \
		--dicoms ./dicoms_$MIN_SAMPLE_ID/ \
		--min_sample_id $MIN_SAMPLE_ID \
		--max_sample_id $MAX_SAMPLE_ID &

    let COUNTER=COUNTER+1
    let MIN_SAMPLE_ID=MIN_SAMPLE_ID+INCREMENT
    let MAX_SAMPLE_ID=MAX_SAMPLE_ID+INCREMENT

    sleep 2s
done

################### DISPLAY TIME #########################################

END_TIME=$(date +%s)
ELAPSED_TIME=$(($END_TIME - $START_TIME))
printf "\nDispatched $((COUNTER - 1)) tensorization jobs in "
display_time $ELAPSED_TIME
