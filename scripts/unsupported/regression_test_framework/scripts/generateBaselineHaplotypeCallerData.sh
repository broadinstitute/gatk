#!/usr/bin/env bash

GATK4="broadinstitute/gatk:4.0.0.0"



function generateJson() 
{
  local gatk_docker=${1}

  echo "{"
  echo "  \"ToolComparisonWdl.gatk_docker\": \"${1}\","
  echo "  \"ToolComparisonWdl.analysis_docker\": \"broadinstitute/gatk-nightly:2019-02-26-4.1.0.0-31-g23bd0a2f8-SNAPSHOT\","
  echo ""
  echo "  \"ToolComparisonWdl.input_bucket_location\": \"gs://broad-dsp-methods-regression-testing/inputData/\","
  echo "  \"ToolComparisonWdl.truth_bucket_location\": \"gs://broad-dsp-methods-regression-testing/inputData/\","
  echo ""
  echo "  \"ToolComparisonWdl.input_bams\": [ \"NexPond-359781.bam\", \"NexPond-359877.bam\", \"NexPond-360361.bam\", \"NexPond-360457.bam\", \"NexPond-361337.bam\", \"NexPond-361433.bam\", \"NexPond-362428.bam\", \"NexPond-363907.bam\", \"NexPond-445394.bam\" ],"
  echo ""
  echo "  \"ToolComparisonWdl.contamination\": 0.0,"
  echo "  \"ToolComparisonWdl.interval_padding\": 0,"
  echo "  \"ToolComparisonWdl.gvcf_mode\": \"false\","
  echo ""
  echo "  \"ToolComparisonWdl.disk_space_gb\": 512,"
  echo "  \"ToolComparisonWdl.boot_disk_size_gb\": 64,"
  echo "  \"ToolComparisonWdl.mem_gb\": 32"
  echo "}  "
}



################################################################################

BADLIST=""

for gatkVersion in \
  2019-03-07-4.1.0.0-50-g342569ac0-SNAPSHOT  \
  2019-03-08-4.0.4.0-SNAPSHOT \
  2019-03-08-4.0.6.0-SNAPSHOT \
  2019-03-08-4.0.7.0-SNAPSHOT \
  2019-03-08-4.0.8.0-SNAPSHOT \
  2019-03-08-4.0.8.1-SNAPSHOT \
  2019-03-08-4.0.9.0-SNAPSHOT \
  2019-03-08-4.0.10.1-SNAPSHOT \
  2019-03-08-4.0.10.0-SNAPSHOT \
  2019-03-08-4.0.12.0-SNAPSHOT \
  2019-03-08-4.1.0.0-SNAPSHOT  \
  ; do

  gatkDockerImage="jonnsmith/gatk_test_builds:${gatkVersion}"

  echo "Creating json for ${gatkVersion}"
  generateJson ${gatkDockerImage} > ${gatkVersion}.json
  
  # Submit the job:
  echo "Submitting job for ${gatkVersion}"
  ~/Development/cromshell/cromshell submit -w \
    ~/Development/gatk/scripts/unsupported/regression_test_framework/wdl/compareHaplotypeCallerRuns.wdl \
    ${gatkVersion}.json \
    ~/Development/gatk/scripts/unsupported/regression_test_framework/json/options.json \
    ~/Development/gatk/scripts/unsupported/regression_test_framework/wdl/compareHaplotypeCallerRuns_Subworkflows.zip

  if [[ $? -ne 0 ]] ; then
    BADLIST="${BADLIST} ${gatkVersion}"
  fi

  rm ${gatkVersion}.json
  sleep 2

  echo ""

done

if [[ ${#BADLIST} -ne 0 ]] ; then
  echo "YOU MUST RESUBMIT THE FOLLOWING JOBS:"
  echo "${BADLIST}"
fi

echo "Done"

