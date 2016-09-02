#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

java -jar ~/cromwell-0.19.3.jar run /home/travis/build/broadinstitute/gatk-protected/scripts/wdl/pon_gatk_workflow.wdl pon_gatk_workflow_exome.json
java -jar ~/cromwell-0.19.3.jar run /home/travis/build/broadinstitute/gatk-protected/scripts/wdl/pon_gatk_workflow.wdl pon_gatk_workflow_wgs.json
java -jar ~/cromwell-0.19.3.jar run /home/travis/build/broadinstitute/gatk-protected/scripts/wdl/case_gatk_acnv_workflow.wdl case_gatk_acnv_workflow_wgs.json
java -jar ~/cromwell-0.19.3.jar run /home/travis/build/broadinstitute/gatk-protected/scripts/wdl/case_gatk_acnv_workflow.wdl case_gatk_acnv_workflow_exome.json
