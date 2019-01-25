#!/usr/bin/env bash

# Add in the data files first:

# Create sql insert statement:
insertionQuery="INSERT INTO DSPRegressionTesting.InputFiles (path, type, md5sum) VALUES " 

INPUTFILEURL=gs://broad-dsp-methods-regression-testing/inputData

echo "Creating input file list from ${INPUTFILEURL}"
for f in $( gsutil ls ${INPUTFILEURL} ) ; do
  if [[ ${f} =~ ^.*.txt$ ]] ; then
    continue
  fi

  FT=""
  bn=$( basename $f )
  if [[ $bn =~ ^.*.bam$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="GENOME_READS"
  elif [[ $bn =~ ^.*.bam$ ]] ; then
    FT="EXOME_READS"
  elif [[ $bn =~ ^.*.bai$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="GENOME_READS_INDEX"
  elif [[ $bn =~ ^.*.bai$ ]] ; then
    FT="EXOME_READS_INDEX"
  elif [[ $bn =~ ^.*.vcf.gz$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.gz$ ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF_INDEX"
  elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] ; then
    FT="SOMATIC_VCF_INDEX"
  fi

  MD5=$( grep "${bn}[ \t]" FileProvenance.txt | awk '{print $2}' )
  if [[ ${#MD5} -eq 0 ]] ; then
    echo "Could not get MD5 for $bn - IGNORING" 1>&2
    continue
  fi

  insertionQuery="${insertionQuery}('${f}', '${FT}', '${MD5}'), "
done

# Get rid of trailing ", " and add a semicolon:
insertionQuery=$( echo "${insertionQuery}" | sed 's#, $#;#g' )

echo
echo "About to execute query.  You will be asked for the DB password for DSPRegressionTesting" 
echo "Query:"
echo ${insertionQuery}

mysql -h mysql-prd2.broadinstitute.org -u DSPRegressionTesting -p -e "${insertionQuery}"

echo "Done."
################################################################################
################################################################################
################################################################################

# Create the Tools:


