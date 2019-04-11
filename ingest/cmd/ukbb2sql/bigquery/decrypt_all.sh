#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

#this script takes a folder of .enc files from UKBB, keys, and produces .csv.gz
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
enc_directory=${__dir}/uk_biobank_4_1_2019

for file in 28112 23300 23301 23302 
do
	${__dir}/ukbunpack ${enc_directory}/ukb${file}.enc ${__dir}/k17488_${file}.key
	${__dir}/ukbconv ${enc_directory}/ukb${file}.enc_ukb csv
	gzip ${enc_directory}/ukb${file}.csv
done
