#!/bin/bash
# NOTE: this script is intended to be placed in google cloud storage
# and invoked by adding the following line to your cromwell workflow
# options:
#    "monitoring_script": "gs://bucket/path/to/cromwell_monitoring_script.sh"
# Upon task completion "monitoring.log" will be added to the appropriate
# cloud storage folder.
set -Eeuo pipefail

MONITOR_MOUNT_POINT=${MONITOR_MOUNT_POINT:-"/cromwell_root"}
SLEEP_TIME=${SLEEP_TIME:-"10"}

function getCpuUsage() {
    # get the summary cpu statistics (i.e. for all cpus) since boot
    # get the numeric values in an array, dropping the first field (the
    # string, "cpu")
    CPU_TIMES=(`sed -n 's/^cpu\s//p' /proc/stat`)
    # idle time (in system units) is the 3rd numeric field
    IDLE_TIME=${CPU_TIMES[3]}
    # total cpu time is sum of all fields
    TOTAL_TIME=0
    for T in ${CPU_TIMES[@]}; do
        ((TOTAL_TIME += T))
    done
    
    # get the previous times from temp file
    read PREVIOUS_IDLE PREVIOUS_TOTAL < $TEMP_CPU
    
    # write current times to temp file
    echo "$IDLE_TIME $TOTAL_TIME" > $TEMP_CPU
    
    # get the difference in idle and total times since the previous
    # update, and report the usage as: non-idle time as a percentage
    # of total time
    awk -v IDLE=$((IDLE_TIME-PREVIOUS_IDLE)) \
        -v TOTAL=$((TOTAL_TIME-PREVIOUS_TOTAL)) \
        'BEGIN { printf "%.1f%%", 100 * (1 - IDLE / TOTAL) }'
}

function getMem() {
    # get desired memory value from /proc/meminfo, in GiB, and also
    # as a percentage of total
    # argument is the label of the desired memory value
    cat /proc/meminfo \
        | awk -v MEM_FIELD="$1" '{
            f[substr($1, 1, length($1)-1)] = $2
        } END {
            printf "%.2f GiB", f[MEM_FIELD] / 1048576
        }'     
}

function getMemUnavailable() {
    # get unavailable memory from /proc/meminfo, in GiB
    cat /proc/meminfo \
        | awk '{
            f[substr($1, 1, length($1)-1)] = $2
        } END {
            
            if("MemAvailable" in f) {
                mem_available = f["MemAvailable"]
            } else {
                mem_available = f["MemFree"] + f["Buffers"] + f["Cached"]
            }
            mem_in_use = f["MemTotal"] - mem_available
            printf "%.2f GiB %.1f%%", mem_in_use / 1048576, 100 * mem_in_use / f["MemTotal"] 
        }'    
}

# old version using "free -m" are kept in case a container somehow has
# weird values in /proc/meminfo
function getMem_with_free() {
    # get memory info from "free" command. Convert to float in GB.
    # First argument is desired row of output table.
    # Second argument is desired column.
    MEM_ROW=$(echo "$1" | awk '{print tolower($1)}')
    MEM_COLUMN=$(echo "$2" | awk '{print tolower($1)}')
    free -m | awk -v MEM_ROW=$MEM_ROW -v MEM_COLUMN=$MEM_COLUMN \
        'NR=1 {
            for(i=1; i<=NF; i++) { f[tolower($i)]=NF+1-i }
        }
        {
            regex="^"MEM_ROW
            if(tolower($1) ~ regex) {
                print $(NF+1-f[MEM_COLUMN])/1024 " GiB"
            }
        }'
}

# old version using "free -m" are kept in case a container somehow has
# weird values in /proc/meminfo
function getMemUnavailable_using_free() {
    # get memory that is in active use (not just cached) from "free"
    # command. Convert to float in GiB, followed by percent of total.
    # NOTE: weird computation with awk due to variety of output from
    # free on different systems. Rows and columns differ, and on some
    # systems the desired quantity is used "used" memory, on most it's
    # "used" - "buffers" - "cached". If "buffers" and "cached" don't
    # exist, then awk will subtract 0 so the correct result is returned.
    free -m \
        | awk '\
            NR=1 {
                for(i=1; i<=NF; i++) { f[tolower($i)]=NF+1-i }
            }
            {
                if(tolower($1) ~ "^mem") {
                    IN_USE=($(NF+1-f["used"]) - $(NF+1-f["buffers"]) - $(NF+1-f["cached"]))
                    printf "%.3f GiB %.1f%%", IN_USE/1024, 100*IN_USE/$(NF+1-f["total"])
                }
            }'
}


function getDisk() {
    # get information about disk usage from "df" command.
    DISK_COLUMN=$(echo "$1" | awk '{print tolower($1)}')
    MOUNT_POINT=$2
    # extract desired value
    VALUE=$(\
        df -h "$MOUNT_POINT" \
        | sed 's/Mounted on/Mounted-on/' \
        | awk -v DISK_COLUMN=$DISK_COLUMN '
            FNR==1 {
                NF_HEADER=NF
                for(i=1; i<=NF; i++) { f[tolower($i)]=NF-i }
            }
            FNR>1 {
                FIELD_NUM=NF-f[DISK_COLUMN]
                if(FIELD_NUM > 0) {
                    VALUE=$(FIELD_NUM)
                    print VALUE
                } else if(f[DISK_COLUMN] == NF_HEADER-1 && NF == 1) {
                    VALUE=$(1)
                    print VALUE
                }
            }' \
    )
    # If value is a number follwed by letters, it is a value with units
    # and needs to be converted. Otherwise just print value
    if [[ "$VALUE" =~ [0-9.]+[A-z]+ ]]; then
        echo "$VALUE"\
        | sed -E 's/([0-9.]*)([^0-9.]*)/\1 \2/' \
        | awk '{
            UNIT=substr($2, 1, 1)
            if(UNIT == "T") {
                SCALE=2^10
            } else if(UNIT == "G") {
                SCALE=1
            } else if(UNIT == "M") {
                SCALE=2^-10
            } else if(UNIT == "K") {
                SCALE=2^-20
            } else if(UNIT == "B") {
                SCALE=2^-30
            } else {
                SCALE=1
            }
            printf "%.3f GiB", $1 * SCALE
        }'
    else
        echo "$VALUE"
    fi
}

function findBlockDevice() {
    MOUNT_POINT=$1
    FILESYSTEM=$(grep -E "$MOUNT_POINT\s" /proc/self/mounts \
                | awk '{print $1}')
    DEVICE_NAME=$(basename "$FILESYSTEM")
    FS_IN_BLOCK=$(find -L /sys/block/ -mindepth 2 -maxdepth 2 -type d \
                       -name "$DEVICE_NAME")
    if [ -n "$FS_IN_BLOCK" ]; then
        # found path to the filesystem in the block devices. get the
        # block device as the parent dir
        dirname "$FS_IN_BLOCK"
    elif [ -d "/sys/block/$DEVICE_NAME" ]; then
        # the device is itself a block device
        echo "/sys/block/$DEVICE_NAME"
    else
        # couldn't find, possibly mounted by mapper.
        # look for block device that is just the name of the symlinked
        # original file. if not found, echo empty string (no device found)
        BLOCK_DEVICE=$(ls -l "$FILESYSTEM" 2>/dev/null \
                        | cut -d'>' -f2 \
                        | xargs basename 2>/dev/null \
                        || echo)
        if [[ -z "$BLOCK_DEVICE" ]]; then
            1>&2 echo "Unable to find block device for filesystem $FILESYSTEM."
            if [[ -d /sys/block/sdb ]] && ! grep -qE "^/dev/sdb" /etc/mtab; then
                1>&2 echo "Guessing present but unused sdb is the correct block device."
                echo "/sys/block/sdb"
            else           
                1>&2 echo "Disk IO will not be monitored."
            fi
        fi
    fi
}

function handle_integer_wrap() {
    if [ $1 -ge 0 ]; then
        echo $1
    else
        WRAPPED=$1
        echo "$((WRAPPED + 2**30))"
    fi
}



function getBlockDeviceIO() {
    # get read and write IO rate by looking at appropriate block device
    STAT_FILE="$1"
    if [[ -f "$STAT_FILE" ]]; then
        # get IO stats as comma-separated list to extract 3rd and 7th fields
        STATS=$(sed -E 's/[[:space:]]+/,/g' $STAT_FILE | sed -E 's/^,//'\
                | cut -d, -f3,7 | sed -E 's/,/ /g')
        # get results of previous poll
        read OLD_READ OLD_WRITE < $TEMP_IO
        # save new poll results
        read READ_SECTORS WRITE_SECTORS <<<$STATS
        echo "$READ_SECTORS $WRITE_SECTORS" > $TEMP_IO
        # update read and write sectors as difference since previous poll
        READ_SECTORS=$(handle_integer_wrap $((READ_SECTORS - OLD_READ)))
        WRITE_SECTORS=$(handle_integer_wrap $((WRITE_SECTORS - OLD_WRITE)))

        # output change in read/write sectors in kiB/s
        echo "$READ_SECTORS $WRITE_SECTORS" \
            | awk -v T=$SLEEP_TIME -v B=$SECTOR_BYTES \
                '{ printf "%.3f MiB/s %.3f MiB/s",  $1*B/T/1048576, $2*B/T/1048576 }'
    else
        echo "N/A MiB/s N/A MiB/s"
    fi
}


function runtimeInfo() {
    echo [$(date)]
    echo \* CPU usage: $(getCpuUsage)
    echo \* Memory usage: $(getMemUnavailable)
    echo \* Disk usage: $(getDisk Used $MONITOR_MOUNT_POINT) $(getDisk Use% $MONITOR_MOUNT_POINT)
    echo \* Read/Write IO: $(getBlockDeviceIO "$BLOCK_DEVICE_STAT_FILE")
}

# print out header info
echo ==================================
echo =========== MONITORING ===========
echo ==================================
echo --- General Information ---
echo \#CPU: $(nproc)
echo Total Memory: $(getMem MemTotal)
echo Total Disk space: $(getDisk Size "$MONITOR_MOUNT_POINT")
echo 
echo --- Runtime Information ---


# make a temp file to store io information, remove it on exit
TEMP_IO=$(mktemp "${TMPDIR:-/tmp/}$(basename $0).XXXXXXXXXXXX")
# make a temp file to store cpu information, remove it on exit
# remove temp files on exit
TEMP_CPU=$(mktemp "${TMPDIR:-/tmp/}$(basename $0).XXXXXXXXXXXX")
trap "rm -f $TEMP_IO $TEMP_CPU" EXIT


# find the block device
BLOCK_DEVICE=$(findBlockDevice "$MONITOR_MOUNT_POINT")
if [[ -z "$BLOCK_DEVICE" ]] \
        || [[ ! -f "$BLOCK_DEVICE/queue/hw_sector_size" ]]; then
    # no block device found, can't get IO info
    SECTOR_BYTES=0
    BLOCK_DEVICE_STAT_FILE=""
else
    SECTOR_BYTES=$(cat "$BLOCK_DEVICE/queue/hw_sector_size")
    BLOCK_DEVICE_STAT_FILE="$BLOCK_DEVICE/stat"
fi


# since getBlockDeviceIO looks at differences in stat file, run the
# update so the first reported update has a sensible previous result to
# compare to
echo "0 0" > $TEMP_IO
getBlockDeviceIO "$BLOCK_DEVICE_STAT_FILE" > /dev/null

# same thing for getCpuUsage
echo "0 0" > $TEMP_CPU
getCpuUsage > /dev/null


while true; do
    runtimeInfo
    sleep $SLEEP_TIME
done
