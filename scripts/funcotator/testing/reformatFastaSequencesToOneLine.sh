#!/usr/bin/env bash

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
