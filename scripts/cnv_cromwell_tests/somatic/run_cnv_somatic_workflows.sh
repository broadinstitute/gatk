#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/cnv_somatic_tasks.wdl
ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_copy_ratio_bam_workflow.wdl
ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_allele_fraction_pair_workflow.wdl

CROMWELL_JAR="cromwell-0.26.jar"

# Panel WES
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wes_workflow.json
# Panel WGS
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wgs_workflow.json

# Pair WES
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wes_workflow.json
# Pair WGS
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wgs_workflow.json
# Pair WES tumor-only
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wes_tumor-only_workflow.json
# Pair WGS tumor-only
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wgs_tumor-only_workflow.json