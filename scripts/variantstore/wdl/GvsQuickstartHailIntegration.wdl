version 1.0

import "GvsUtils.wdl" as Utils
import "GvsExtractAvroFilesForHail.wdl" as ExtractAvroFilesForHail
import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration

workflow GvsQuickstartHailIntegration {
    input {
        String branch_name
        String hail_wheel = "gs://gvs-internal-scratch/hail-wheels/2022-10-18/0.2.102-964bee061eb0/hail-0.2.102-py3-none-any.whl"
    }

    String project_id = "gvs-internal"

    call QuickstartVcfIntegration.GvsQuickstartVcfIntegration {
        input:
            branch_name = branch_name,
            drop_state = "NONE",
            extract_do_not_filter_override = false,
            dataset_suffix = "hail",
    }

    call ExtractAvroFilesForHail.GvsExtractAvroFilesForHail {
        input:
            go = GvsQuickstartVcfIntegration.done,
            project_id = project_id,
            dataset_name = GvsQuickstartVcfIntegration.dataset_name,
            filter_set_name = GvsQuickstartVcfIntegration.filter_set_name,
            scatter_width = 10,
    }

    call CreateAndTieOutVds {
        input:
            branch_name = branch_name,
            hail_wheel = hail_wheel,
            avro_prefix = GvsExtractAvroFilesForHail.avro_prefix,
            vds_destination_path = GvsExtractAvroFilesForHail.vds_output_path,
            tieout_vcfs = GvsQuickstartVcfIntegration.output_vcfs,
            tieout_vcf_indexes = GvsQuickstartVcfIntegration.output_vcf_indexes,
    }

    output {
        Array[File] output_vcfs = GvsQuickstartVcfIntegration.output_vcfs
        Array[File] output_vcf_indexes = GvsQuickstartVcfIntegration.output_vcf_indexes
        Float total_vcfs_size_mb = GvsQuickstartVcfIntegration.total_vcfs_size_mb
        File manifest = GvsQuickstartVcfIntegration.manifest
        String vds_output_path = GvsExtractAvroFilesForHail.vds_output_path
        Boolean done = true
    }
}


task CreateAndTieOutVds {
    input {
        File hail_wheel
        String branch_name
        String avro_prefix
        String vds_destination_path
        Array[File] tieout_vcfs
        Array[File] tieout_vcf_indexes
    }
    parameter_meta {
        tieout_vcfs: {
             localization_optional: true
         }
        tieout_vcf_indexes: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Copy the versions of the Hail import and tieout scripts for this branch from GitHub.
        script_url_prefix="https://raw.githubusercontent.com/broadinstitute/gatk/~{branch_name}/scripts/variantstore/wdl/extract"
        for script in hail_gvs_import.py hail_join_vds_vcfs.py gvs_vds_tie_out.py
        do
            curl --silent --location --remote-name "${script_url_prefix}/${script}"
        done

        # Create a manifest of VCFs and indexes to bulk download with `gcloud storage cp`.
        touch vcf_manifest.txt
        # This is extremely noisy and not interesting, turn off xtrace.
        set +o xtrace
        for file in ~{sep=' ' tieout_vcfs} ~{sep=' ' tieout_vcf_indexes}
        do
            echo $file >> vcf_manifest.txt
        done
        # xtrace back on
        set -o xtrace

        # Copy VCFs and indexes to the current directory.
        cat vcf_manifest.txt | gcloud storage cp -I .

        # `avro_prefix` includes a trailing `avro` so don't add another `avro` here.
        gcloud storage cp --recursive ~{avro_prefix} $PWD

        export REFERENCES_PATH=$PWD/references
        mkdir -p ${REFERENCES_PATH}

        gcloud storage cp 'gs://hail-common/references/Homo_sapiens_assembly38.fasta*' ${REFERENCES_PATH}

        # Temurin Java 8
        apt-get -qq install wget apt-transport-https gnupg
        wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
        echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
        apt-get -qq update
        apt -qq install -y temurin-8-jdk

        pip install ~{hail_wheel}
        export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'

        export WORK=$PWD/work
        mkdir ${WORK}

        export TEMP_PATH=$WORK/temp
        mkdir ${TEMP_PATH}

        export VDS_PATH=$WORK/gvs_import.vds
        export AVRO_PATH=$PWD/avro

        python3 ./hail_gvs_import.py --avro-path ${AVRO_PATH} --vds-path ${VDS_PATH} --temp-path ${TEMP_PATH} --references-path ${REFERENCES_PATH}

        export JOINED_MATRIX_TABLE_PATH=${WORK}/joined.mt

        python3 ./hail_join_vds_vcfs.py --vds-path ${VDS_PATH} --joined-matrix-table-path ${JOINED_MATRIX_TABLE_PATH} *.vcf.gz

        # Copy up the VDS
        gcloud storage cp --recursive ${VDS_PATH} ~{vds_destination_path}

        pip install pytest
        ln -s ${WORK}/joined.mt .
        pytest ./gvs_vds_tie_out.py
    >>>
    runtime {
        # `slim` here to be able to use Java
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:409.0.0-slim"
        disks: "local-disk 2000 HDD"
        memory: "30 GiB"
    }
    output {
        Boolean done = true
    }
}
