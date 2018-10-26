#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script takes a FASTA file and reformats it such that the sequences 
# themselves are all on a single line below the sequence names.
# This is useful when you need to batch process FASTA files using shell scripts.
# 
# There is a bit of awk wizardry of which I am particularly proud and have no
# qualms revealing that I was inspired by some bits on Stack Overflow, though 
# these bits were amended to work for this purpose.  
#
# EXAMPLE:
#     ./reformatFastaSequencesToOneLine.sh FASTA_FILE.fa
#
# AUTHOR: Jonn Smith
#
###############################################################################

INFILE=$1

[[ ! -f $INFILE ]] && echo "ERROR: File does not exist: $INFILE" 1>&2 && exit 1

## Inspired by https://stackoverflow.com/questions/33120215/remove-all-occurrence-of-new-line-between-two-patterns-sed-or-awk
awk 'BEGIN { first=1; }
      /^>/ { 
				if (!first) { 
					print "";
				} 
				print; 
				first = 0; 
				next; 
			}
      { 
				printf sep""$0; 
			}' ${INFILE}
