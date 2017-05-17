#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

ln -fs /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/cnv_common_tasks.wdl
ln -fs /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_panel_creation_workflow.wdl
ln -fs /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_single_sample_calling_workflow.wdl
ln -fs /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_cohort_calling_workflow.wdl

CROMWELL_JAR="cromwell-0.26.jar"

# Panel WES
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_panel_creation_workflow.wdl gCNV_panel_creation_workflow_wes.json
# Panel WGS
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_panel_creation_workflow.wdl gCNV_panel_creation_workflow_wgs.json

# Single sample WES calling
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_single_sample_calling_workflow.wdl gCNV_single_sample_calling_workflow_wes.json
# Single sample WGS calling
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_single_sample_calling_workflow.wdl gCNV_single_sample_calling_workflow_wgs.json

# Cohort WES calling
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_cohort_calling_workflow.wdl gCNV_cohort_calling_workflow_wes.json
# Cohort WGS calling
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk-protected/scripts/cnv_wdl/germline/gCNV_cohort_calling_workflow.wdl gCNV_cohort_calling_workflow_wgs.json