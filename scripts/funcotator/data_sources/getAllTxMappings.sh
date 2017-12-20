#!/usr/bin/env bash

# Creates a set of files that map records between many known databases of genes.
# Pulled directly from the ensemble database.
# Currently not used by Funcotator.

outFileBaseName="all_TX_mappings"
outExt=".tsv"

hg19db="homo_sapiens_core_75_37"
hg38db="homo_sapiens_core_90_38"

hg19FileName=${outFileBaseName}_v75_37.hg19${outExt}
hg38FileName=${outFileBaseName}_v90_38.hg38${outExt}

################################################################################

query="SELECT CONCAT(transcript.stable_id, '.', transcript.version) AS transcript_id, xref.display_label, external_db.db_name FROM translation, transcript, object_xref, xref,external_db WHERE transcript.transcript_id = translation.transcript_id AND translation.translation_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = 'Translation' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id;"

echo "Getting HG19 TX Mappings..." 
[ -f ${hg19FileName} ] && rm ${hg19FileName}
time mysql -u anonymous -h ensembldb.ensembl.org -e "use ${hg19db};${query}" | awk 'NR<2{print $0;next}{print $0| "sort -r"}' >> ${hg19FileName}

echo "Getting HG38 TX Mappings..." 
[ -f ${hg38FileName} ] && rm ${hg38FileName}
time mysql -u anonymous -h ensembldb.ensembl.org -e "use ${hg38db};${query}" | awk 'NR<2{print $0;next}{print $0| "sort -r"}' >> ${hg38FileName}

echo 'Done!'
