version 1.0

import "../GvsUtils.wdl" as Utils

workflow GvsTieOutVds {
    input {
        String vds_path
        Array[File] tieout_vcfs
        Array[File] tieout_vcf_indexes
        String tieout_vcf_suffix = ".gz"
        String? git_branch_or_tag
        String? cloud_sdk_slim_docker
        String? hail_version
    }

    if (!defined(git_branch_or_tag) || !defined(cloud_sdk_slim_docker) || !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])
    String effective_git_branch_or_tag = select_first([git_branch_or_tag, GetToolVersions.git_hash])

    call TieOutVDS {
        input:
            go = true,
            tieout_vcfs = tieout_vcfs,
            tieout_vcf_indexes = tieout_vcf_indexes,
            vds_path = vds_path,
            tieout_vcf_suffix = tieout_vcf_suffix,
            git_branch_or_tag = effective_git_branch_or_tag,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            hail_version = effective_hail_version,
    }
}

task TieOutVDS {
    input {
        Boolean go
        String git_branch_or_tag
        String vds_path
        Array[File] tieout_vcfs
        Array[File] tieout_vcf_indexes
        String tieout_vcf_suffix
        String cloud_sdk_slim_docker
        String hail_version
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
        script_url_prefix="https://raw.githubusercontent.com/broadinstitute/gatk/~{git_branch_or_tag}/scripts/variantstore/scripts"
        for script in hail_gvs_import.py hail_gvs_util.py hail_join_vds_vcfs.py gvs_vds_tie_out.py import_gvs.py import_gvs_ploidy.py
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

        export REFERENCES_PATH=$PWD/references
        mkdir -p ${REFERENCES_PATH}

        gcloud storage cp 'gs://hail-common/references/Homo_sapiens_assembly38.fasta*' ${REFERENCES_PATH}

        # Versions of Hail near 0.2.117 demand Java 8 or Java 11, and refuse to run on Java 17. (This is because Google Dataproc is still on Java 11)
        # Temurin Java 8
        apt-get -qq install wget apt-transport-https gnupg
        wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
        echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
        apt-get -qq update
        apt -qq install -y temurin-8-jdk

        export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
        pip install --upgrade pip
        pip install hail==~{hail_version}

        export WORK=$PWD/work
        mkdir ${WORK}

        export TEMP_PATH=$WORK/temp
        mkdir ${TEMP_PATH}

        export VDS_PATH=$WORK/gvs_export.vds
        mkdir ${VDS_PATH}

        gcloud storage cp -r ~{vds_path}  ${WORK}


        export JOINED_MATRIX_TABLE_PATH=${WORK}/joined.mt

        python3 ./hail_join_vds_vcfs.py --vds-path ${VDS_PATH} --joined-matrix-table-path ${JOINED_MATRIX_TABLE_PATH} *.vcf~{tieout_vcf_suffix}

        pip install pytest
        ln -s ${WORK}/joined.mt .
        pytest ./gvs_vds_tie_out.py
    >>>
    runtime {
        # `slim` here to be able to use Java
        docker: cloud_sdk_slim_docker
        maxRetries: 2
        disks: "local-disk 2000 HDD"
        memory: "30 GiB"
    }
    output {
        Boolean done = true
    }
}
