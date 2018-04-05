#!/bin/bash

set -eu

#################################################
function ins_against_ins {

	echo "insertion against PacBio insertion"
	echo " $1" # call set to be evaluated
	echo " $2" # PacBio callset;

	POSTFIX="$3_$4"

	WINDOW_SIZE=50

	# intersection strategy: up/downstream 50 bases from call set to be evaluated
	# first check what the intersection scenario is, because a variable may "intersect" with several from PacBio
	bedtools window \
		-a "$1" -b "$2" \
		-w "$WINDOW_SIZE" -c | \
		awk '{print $NF}' | sort -n | uniq | xargs echo

	bedtools window \
		-a "$1" -b "$2" \
		-w "$WINDOW_SIZE" > \
		"$OUTPUT_DIR""insIns.""$POSTFIX"".window_w50.txt"
	wc -l "$OUTPUT_DIR""insIns.""$POSTFIX"".window_w50.txt" | awk '{print $1}'

	# output 
	if [[ $5 == "gatk" ]]; then
		awk 'BEGIN {OFS="	";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""insIns.""$POSTFIX"".window_w50.txt" > \
			"$OUTPUT_DIR""insIns.""$POSTFIX"".window_w50.IDs.txt"
	else
		awk 'BEGIN {OFS="	";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""insIns.""$POSTFIX"".window_w50.txt" > \
			"$OUTPUT_DIR""insIns.""$POSTFIX"".window_w50.IDs.txt"
	fi
	echo
}

function dup_against_ins {

	echo "tandup against PacBio insertion"
	echo "$1"
	echo "$2"

	POSTFIX="$3_$4"

	WINDOW_SIZE=50

	# intersection strategy: up/downstream 0/50 bases from call set to be evaluated
	# first check what the intersection scenario is, because a variable may "intersect" with several from PacBio
	bedtools window \
		-a "$1" -b "$2" \
		-l 0 -r "$WINDOW_SIZE" -c | \
		awk '{print $NF}' | sort -n | uniq | xargs echo

	bedtools window \
		-a "$1" -b "$2" \
		-l 0 -r "$WINDOW_SIZE" > \
		"$OUTPUT_DIR""dupIns.""$POSTFIX"".window_l0r50.txt"
	wc -l "$OUTPUT_DIR""dupIns.""$POSTFIX"".window_l0r50.txt" | awk '{print $1}'

	if [[ $5 == "gatk" ]]; then
		awk 'BEGIN {OFS="	";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""dupIns.""$POSTFIX"".window_l0r50.txt" > \
			"$OUTPUT_DIR""dupIns.""$POSTFIX"".window_l0r50.IDs.txt"
	else
		awk 'BEGIN {OFS="	";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""dupIns.""$POSTFIX"".window_l0r50.txt" > \
			"$OUTPUT_DIR""dupIns.""$POSTFIX"".window_l0r50.IDs.txt"
	fi
	echo
}

function rpl_against_ins {

	echo "replacement against PacBio insertion"
	echo "$1"
	echo "$2"

	POSTFIX="$3_$4"

	WINDOW_SIZE=50

	# intersection strategy: up/downstream 50 bases from call set to be evaluated
	# first check what the intersection scenario is, because a variable may "intersect" with several from PacBio
	bedtools window \
		-a "$1" -b "$2" \
		-w "$WINDOW_SIZE" -c | \
		awk '{print $NF}' | sort -n | uniq | xargs echo

	bedtools window \
		-a "$1" -b "$2" \
		-w "$WINDOW_SIZE" > \
		"$OUTPUT_DIR""rplIns.""$POSTFIX"".window_w50.txt"
	wc -l "$OUTPUT_DIR""rplIns.""$POSTFIX"".window_w50.txt" | awk '{print $1}'

	if [[ $5 == "gatk" ]]; then
		awk 'BEGIN {OFS="	";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""rplIns.""$POSTFIX"".window_w50.txt" > \
			"$OUTPUT_DIR""rplIns.""$POSTFIX"".window_w50.IDs.txt"
	else
		awk 'BEGIN {OFS="	";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""rplIns.""$POSTFIX"".window_w50.txt" > \
			"$OUTPUT_DIR""rplIns.""$POSTFIX"".window_w50.IDs.txt"
	fi
	echo
}


function rpl_against_del {

	echo "replacement against PacBio deletion"
	echo "$1"
	echo "$2"

	POSTFIX="$3_$4"

	WINDOW_SIZE=50

	# intersection strategy: up/downstream 50 bases from call set to be evaluated
	# first check what the intersection scenario is, because a variable may "intersect" with several from PacBio
	bedtools window \
		-a "$1" -b "$2" \
		-w "$WINDOW_SIZE" -c | \
		awk '{print $NF}' | sort -n | uniq | xargs echo

	bedtools window \
		-a "$1" -b "$2" \
		-w "$WINDOW_SIZE" > \
		"$OUTPUT_DIR""rplDel.""$POSTFIX"".window_w50.txt"
	wc -l "$OUTPUT_DIR""rplDel.""$POSTFIX"".window_w50.txt" | awk '{print $1}'

	if [[ $5 == "gatk" ]]; then
		awk 'BEGIN {OFS="	";} {id=$9":"$10; print $3, id}' "$OUTPUT_DIR""rplDel.""$POSTFIX"".window_w50.txt" > \
			"$OUTPUT_DIR""rplDel.""$POSTFIX"".window_w50.IDs.txt"
	else
		awk 'BEGIN {OFS="	";} {id=$11":"$12; print $3, id}' "$OUTPUT_DIR""rplDel.""$POSTFIX"".window_w50.txt" > \
			"$OUTPUT_DIR""rplDel.""$POSTFIX"".window_w50.IDs.txt"
	fi
	echo
}