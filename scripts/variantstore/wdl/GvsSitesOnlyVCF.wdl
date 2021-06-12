version 1.0
workflow GvsSitesOnlyVCF {
   input {
        Array[File] gvs_extract_cohort_filtered_vcfs
        String output_sites_only_file_name
        String output_merged_file_name
        String output_annotated_file_name
        File nirvana_data_directory
        File? gatk_override
    }

    scatter(i in range(length(gvs_extract_cohort_filtered_vcfs)) ) {
        call SitesOnlyVcf {
                    input:
                        vcf_bgz_gts                 = gvs_extract_cohort_filtered_vcfs[i],
                        output_filename             = "${output_sites_only_file_name}_${i}.sites_only.vcf.gz",
                }
    }

    call MergeVCFs { # why are we merging before annotating? couldn't we keep them separate?
      input:
          input_vcfs = SitesOnlyVcf.output_vcf,
          input_vcf_indices = SitesOnlyVcf.output_vcf_idx,
          output_merged_file_name = output_merged_file_name,
    }

    call AnnotateVCF {
      input:
          input_vcf = MergeVCFs.merged_vcf,
          output_annotated_file_name = output_annotated_file_name,
          nirvana_data_tar = nirvana_data_directory
    }
}

################################################################################
task SitesOnlyVcf {
    input {
        File vcf_bgz_gts
        String output_filename
    }
    String output_vcf_idx = basename(output_filename) + ".tbi" # or will this be .idx if from .vcf.gz?
    command <<<
        set -e
        gatk --java-options "-Xmx2048m" \
            SelectVariants \
                -V ~{vcf_bgz_gts} \
                --exclude-filtered \
                --sites-only-vcf-output \
                -O ~{output_filename}
     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf="~{output_filename}"
        File output_vcf_idx="~{output_vcf_idx}"
    }
}

task MergeVCFs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_merged_file_name
    }
    String output_vcf = basename(output_merged_file_name) + ".vcf.gz"
    String output_vcf_idx = basename(output_vcf) + ".tbi"
    command <<<
        set -e
        gatk --java-options "-Xmx2048m" \
            MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        File nirvana_data_tar
    }
    String annotation_json_name = basename(output_annotated_file_name) + ".json.gz"
    String annotation_json_name_jsi = annotation_json_name + ".jsi"

    String nirvana_location = "/opt/nirvana/Nirvana.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        set -e

        # NOTE: Validate a lil so that we don't waste time copying down the data sources if there's an error.
        # =======================================
        # Handle our data sources:

        # Extract the tar.gz:
        echo "Extracting annotation data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ~{nirvana_data_tar} -C datasources_dir  --strip-components 2
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"

        dotnet ~{nirvana_location} \
             -c $DATA_SOURCES_FOLDER~{path} \
             --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
             -r $DATA_SOURCES_FOLDER~{path_reference} \
             -i ~{input_vcf} \
             -o ~{output_annotated_file_name}
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "annotation/nirvana:3.14"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotation_json = "~{annotation_json_name}"
        File annotation_json_jsi = "~{annotation_json_name_jsi}"
    }
}


