#!/bin/bash

set -eu

#################################################
function ins_against_ins {

	echo "insertion against PacBio insertion"
	echo " $1"
	echo " $2"

	POSTFIX=$5

	WINDOW_SIZE=50

	bedtools window \
	    -a "$1" -b "$2" \
	    -w "$WINDOW_SIZE" -c | \
		awk '{print $NF}' | sort -n | uniq | xargs echo

	bedtools window \
	    -a "$1" -b "$2" \
	    -w "$WINDOW_SIZE" |
        grep -F "$3" > \
        "$OUTPUT_DIR""insIns.""$4_$POSTFIX"".window_w50.txt"
    wc -l "$OUTPUT_DIR""insIns.""$4_$POSTFIX"".window_w50.txt" | awk '{print $1}'

    if [[ $4 == "GATKvsPacbio" ]]; then
    	awk 'BEGIN {OFS="	";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""insIns.""$4_$POSTFIX"".window_w50.txt" > \
            "$OUTPUT_DIR""insIns.""$4_$POSTFIX"".window_w50.IDs.txt"
    else
		awk 'BEGIN {OFS="	";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""insIns.""$4_$POSTFIX"".window_w50.txt" > \
            "$OUTPUT_DIR""insIns.""$4_$POSTFIX"".window_w50.IDs.txt"
    fi
    echo
}

function dup_against_ins {

	echo "tandup against PacBio insertion"
	echo "$1"
	echo "$2"

	POSTFIX=$5

	WINDOW_SIZE=50

	bedtools window \
    	-a "$1" -b "$2" \
    	-l 0 -r "$WINDOW_SIZE" -c | \
    	awk '{print $NF}' | sort -n | uniq | xargs echo

    bedtools window \
    	-a "$1" -b "$2" \
    	-l 0 -r "$WINDOW_SIZE" | \
    	grep -F "$3" > \
    	"$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.txt"
    wc -l "$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.txt" | awk '{print $1}'

    awk 'BEGIN {OFS="	";} {print $3, $11}' "$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.txt" > \
        "$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.IDs.txt"

    if [[ $4 == "GATKvsPacbio" ]]; then
    	awk 'BEGIN {OFS="	";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.txt" > \
            "$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.IDs.txt"
    else
		awk 'BEGIN {OFS="	";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.txt" > \
            "$OUTPUT_DIR""dupIns.""$4_$POSTFIX"".window_l0r50.IDs.txt"
    fi
    echo
}

function rpl_against_ins {

	echo "replacement against PacBio insertion"
	echo "$1"
	echo "$2"

	POSTFIX=$5

	WINDOW_SIZE=50

	bedtools window \
	    -a "$1" -b "$2" \
	    -w "$WINDOW_SIZE" -c | \
		awk '{print $NF}' | sort -n | uniq | xargs echo

	bedtools window \
	    -a "$1" -b "$2" \
	    -w "$WINDOW_SIZE" |
        grep -F "$3" > \
        "$OUTPUT_DIR""rplIns.""$4_$POSTFIX"".window_w50.txt"
    wc -l "$OUTPUT_DIR""rplIns.""$4_$POSTFIX"".window_w50.txt" | awk '{print $1}'

    if [[ $4 == "GATKvsPacbio" ]]; then
    	awk 'BEGIN {OFS="	";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""rplIns.""$4_$POSTFIX"".window_w50.txt" > \
            "$OUTPUT_DIR""rplIns.""$4_$POSTFIX"".window_w50.IDs.txt"
    else
		awk 'BEGIN {OFS="	";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""rplIns.""$4_$POSTFIX"".window_w50.txt" > \
            "$OUTPUT_DIR""rplIns.""$4_$POSTFIX"".window_w50.IDs.txt"
    fi
    echo
}


function rpl_against_del {

    echo "replacement against PacBio deletion"
    echo "$1"
    echo "$2"

    POSTFIX=$5

    WINDOW_SIZE=50

    bedtools window \
        -a "$1" -b "$2" \
        -w "$WINDOW_SIZE" -c | \
        awk '{print $NF}' | sort -n | uniq | xargs echo

    bedtools window \
        -a "$1" -b "$2" \
        -w "$WINDOW_SIZE" |
        grep -F "$3" > \
        "$OUTPUT_DIR""rplDel.""$4_$POSTFIX"".window_w50.txt"
    wc -l "$OUTPUT_DIR""rplDel.""$4_$POSTFIX"".window_w50.txt" | awk '{print $1}'

    if [[ $4 == "GATKvsPacbio" ]]; then
        awk 'BEGIN {OFS="   ";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""rplDel.""$4_$POSTFIX"".window_w50.txt" > \
            "$OUTPUT_DIR""rplDel.""$4_$POSTFIX"".window_w50.IDs.txt"
    else
        awk 'BEGIN {OFS="   ";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""rplDel.""$4_$POSTFIX"".window_w50.txt" > \
            "$OUTPUT_DIR""rplDel.""$4_$POSTFIX"".window_w50.IDs.txt"
    fi
    echo
}