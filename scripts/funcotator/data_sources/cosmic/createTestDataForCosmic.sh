#!/usr/bin/env bash

# Creates a set of test data for the Funcotator test suite based on
# a COSMIC sqlite3 database file.

LIMIT=500

COSMIC_DB="Cosmic.db"
OUT_CSV_FILE="CosmicTest.csv"
OUT_DB_FILE="CosmicTest.db"

[ -f ${OUT_CSV_FILE} ] && rm ${OUT_CSV_FILE}
[ -f ${OUT_DB_FILE} ] && rm ${OUT_DB_FILE}

sqlite3 Cosmic.db <<EOF
.echo off 
.headers on
.mode csv
.output ${OUT_CSV_FILE} 
SELECT * FROM Cosmic LIMIT ${LIMIT};
EOF

sqlite3 ${OUT_DB_FILE} <<EOF
.echo off
.mode csv
.import ${OUT_CSV_FILE} Cosmic
CREATE INDEX GeneIndex ON Cosmic("Gene name");
VACUUM;
EOF

