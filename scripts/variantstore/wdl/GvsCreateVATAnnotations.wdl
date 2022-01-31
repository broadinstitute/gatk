version 1.0
workflow GvsCreateVATAnnotations {
   input {
        File nirvana_data_directory
        String output_path

        String? service_account_json_path
        File? gatk_override
    }


      ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
        call AnnotateVCF {
          input:
            input_vcf = ExtractAnAcAfFromVCF.output_vcf, ## TODO these need to be grabbed from the bucket
            input_vcf_index = ExtractAnAcAfFromVCF.output_vcf_index, ## TODO these need to be grabbed from the bucket
            output_annotated_file_name = "${input_vcf_name}_annotated",
            nirvana_data_tar = nirvana_data_directory,
            custom_annotations_file = ExtractAnAcAfFromVCF.annotations_file,
        }

        call PrepAnnotationJson {
          input:
            annotation_json = AnnotateVCF.annotation_json,
            output_file_suffix = "${input_vcf_name}.json.gz",
            output_path = output_path,
            service_account_json_path = service_account_json_path
        }


################################################################################

task AnnotateVCF {
    input {
        File input_vcf
        File input_vcf_index
        String output_annotated_file_name
        File nirvana_data_tar
        File custom_annotations_file
    }
    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String annotation_json_name_jsi = annotation_json_name + ".jsi"
    String nirvana_location = "/opt/nirvana/Nirvana.dll"
    String custom_creation_location = "/opt/nirvana/SAUtils.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        set -e

        # =======================================
        # Handle our data sources:

        # Extract the tar.gz:
        echo "Extracting annotation data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ~{nirvana_data_tar} -C datasources_dir  --strip-components 2
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"


        # =======================================
        # Create custom annotations:
        echo "Creating custom annotations"
        mkdir customannotations_dir
        CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

        # Add AC/AN/AF as custom annotations
        ## use --skip-ref once you are on a later version of nirvana
        dotnet ~{custom_creation_location} customvar \
             -r $DATA_SOURCES_FOLDER~{path_reference} \
             -i ~{custom_annotations_file} \
             -o $CUSTOM_ANNOTATIONS_FOLDER

        # =======================================
        # Create Nirvana annotations:


        dotnet ~{nirvana_location} \
             -c $DATA_SOURCES_FOLDER~{path} \
             --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
             --sd $CUSTOM_ANNOTATIONS_FOLDER \
             -r $DATA_SOURCES_FOLDER~{path_reference} \
             -i ~{input_vcf} \
             -o ~{output_annotated_file_name}

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "annotation/nirvana:3.14"
        memory: "32 GB"
        cpu: "2"
        preemptible: 5
        disks: "local-disk 500 SSD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotation_json = "~{annotation_json_name}"
        File annotation_json_jsi = "~{annotation_json_name_jsi}"
    }
}

task PrepAnnotationJson {
    input {
        File annotation_json
        String output_file_suffix
        String output_path
        String? service_account_json_path
    }

    String output_vt_json = "vat_vt_bq_load" + output_file_suffix
    String output_genes_json = "vat_genes_bq_load" + output_file_suffix
    String output_vt_gcp_path = output_path + 'vt/'
    String output_genes_gcp_path = output_path + 'genes/'
    String output_annotations_gcp_path = output_path + 'annotations/'

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    ## TODO these temp files do not currently get cleaned up. Some of them may be helpful for recovery.

    command <<<
        set -e

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_variant_annotation_table.py \
          --annotated_json ~{annotation_json} \
          --output_vt_json ~{output_vt_json} \
          --output_genes_json ~{output_genes_json}

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
        fi

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'
        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

        # for debugging purposes only
        gsutil cp ~{annotation_json} '~{output_annotations_gcp_path}'

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211101"
        memory: "8 GB"
        preemptible: 5
        cpu: "1"
        disks: "local-disk 250 SSD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File vat_vt_json="~{output_vt_json}"
        File vat_genes_json="~{output_genes_json}"
        Boolean done = true
    }
}