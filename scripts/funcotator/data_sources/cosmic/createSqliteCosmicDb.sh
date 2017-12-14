#!/usr/bin/env bash

COSMIC_FILE=CosmicCompleteTargetedScreensMutantExport.tsv
OUT_DB_FILE="Cosmic.db"

[ -f ${OUT_DB_FILE} ] && rm ${OUT_DB_FILE}

sqlite3 ${OUT_DB_FILE} <<EOF
.echo on
.mode tabs
.import ${COSMIC_FILE} RawCosmic
CREATE TABLE Cosmic AS SELECT * FROM RawCosmic WHERE ("Mutation AA" != "" OR "Mutation genome position" != "");
DROP TABLE RawCosmic;
CREATE INDEX GeneIndex ON Cosmic("Gene name");
VACUUM;
EOF

