version 1.0

import "GvsUtils.wdl" as Utils
import "GvsExtractAvroFilesForHail.wdl" as ExtractAvroFilesForHail
import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration

workflow GvsQuickstartHailIntegration {
    input {
        String branch_name
        Boolean use_VQSR_lite = true
        File interval_list
        Boolean use_classic_VQSR = true
        Boolean extract_do_not_filter_override
        String dataset_suffix = "hail"
        String? gatk_override
        String expected_output_prefix
    }

    String project_id = "gvs-internal"

    call QuickstartVcfIntegration.GvsQuickstartVcfIntegration {
        input:
            branch_name = branch_name,
            drop_state = "NONE",
            use_VQSR_lite = use_VQSR_lite,
            extract_do_not_filter_override = extract_do_not_filter_override,
            dataset_suffix = dataset_suffix,
            gatk_override = gatk_override,
            interval_list = interval_list,
            expected_output_prefix = expected_output_prefix,
    }

    call ExtractAvroFilesForHail.GvsExtractAvroFilesForHail {
        input:
            go = GvsQuickstartVcfIntegration.done,
            project_id = project_id,
            use_VQSR_lite = use_VQSR_lite,
            dataset_name = GvsQuickstartVcfIntegration.dataset_name,
            filter_set_name = GvsQuickstartVcfIntegration.filter_set_name,
            scatter_width = 10,
            call_set_identifier = branch_name,
    }

    call CreateAndTieOutVds {
        input:
            branch_name = branch_name,
            use_VQSR_lite = use_VQSR_lite,
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
        String branch_name
        Boolean use_VQSR_lite
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
        for script in hail_gvs_import.py hail_join_vds_vcfs.py gvs_vds_tie_out.py import_gvs.py
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

        export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
        pip install --upgrade pip
        pip install hail

        export WORK=$PWD/work
        mkdir ${WORK}

        export TEMP_PATH=$WORK/temp
        mkdir ${TEMP_PATH}

        export VDS_PATH=$WORK/gvs_import.vds
        export AVRO_PATH=$PWD/avro

        python3 ./hail_gvs_import.py \
            --avro-path ${AVRO_PATH} \
            --vds-path ${VDS_PATH} \
            --temp-path ${TEMP_PATH} \
            --references-path ${REFERENCES_PATH} \
            ~{true='--use-vqsr-lite' false='' use_VQSR_lite}

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
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-slim"
        disks: "local-disk 2000 HDD"
        memory: "30 GiB"
    }
    output {
        Boolean done = true
    }
}
