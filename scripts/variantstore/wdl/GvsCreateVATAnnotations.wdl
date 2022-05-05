version 1.0
workflow GvsCreateVATAnnotations {
   input {
        File input_vcf
        File input_vcf_index
        String input_vcf_name
        File ancestry_mapping_list
        File nirvana_data_directory
        String output_path

        String? service_account_json_path
        File custom_annotations_template
        File ref
    }

      ## Create a sites-only VCF from the original GVS jointVCF
      ## Calculate AC/AN/AF for subpopulations and extract them for custom annotations
        call ExtractAnAcAfFromVCF {
          input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            service_account_json_path = service_account_json_path,
            subpopulation_sample_list = ancestry_mapping_list,
            custom_annotations_template = custom_annotations_template,
            ref = ref,
            output_path = output_path
        }


      ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
        call AnnotateVCF {
          input:
            input_vcf = ExtractAnAcAfFromVCF.output_vcf,
            input_vcf_index = ExtractAnAcAfFromVCF.output_vcf_index,
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

        # ------------------------------------------------
        # Outputs:
        output {
            Int count_variants = ExtractAnAcAfFromVCF.count_variants
            File track_dropped = ExtractAnAcAfFromVCF.track_dropped
            Boolean done = true
        }
}

################################################################################

task ExtractAnAcAfFromVCF {
    input {
        File input_vcf
        File input_vcf_index
        String? service_account_json_path
        File subpopulation_sample_list
        File custom_annotations_template
        File ref
        String output_path
    }
    parameter_meta {
        input_vcf: {
          localization_optional: true
        }
        input_vcf_index: {
          localization_optional: true
        }
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'
    String custom_annotations_file_name = "ac_an_af.tsv"
    String local_input_vcf = basename(input_vcf)
    String local_input_vcf_index = basename(input_vcf_index)
    String normalized_vcf_compressed = "normalized.vcf.gz"
    String normalized_vcf_indexed = "normalized.vcf.gz.tbi"

    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting and calculating the an/ac/af & sc by subpopulation into a tsv
    command <<<
        set -e

        # custom function to prepend the current datetime to an echo statement
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json

        fi

        echo_date "VAT: Custom localization of inputs"

        cp ~{custom_annotations_template} ~{custom_annotations_file_name}

        gsutil cp ~{input_vcf} ~{local_input_vcf}
        gsutil cp ~{input_vcf_index} ~{local_input_vcf_index}
        gsutil cp ~{ref} Homo_sapiens_assembly38.fasta

        # expected_subpopulations = [
        # "afr",
        # "amr",
        # "eas",
        # "eur",
        # "mid",
        # "oth",
        # "sas"
        #]

        echo_date "VAT: Convert input to BCF format"
        bcftools convert --threads 4 -O b -o original.bcf ~{local_input_vcf}
        rm ~{local_input_vcf}

        echo_date "VAT: Calculating number of +50 alt alleles on N sites"

        ## track the dropped variants with +50 alt alleles or N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        bcftools view --threads 4 -i 'N_ALT>50 || REF~"N"' -O u original.bcf | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        echo_date "VAT: filter out sites with too many alt alleles"
        bcftools view --threads 4 -e 'N_ALT>50 || REF~"N"' --no-update original.bcf -O u | \
        ## filter out the non-passing sites
        bcftools view --threads 4 -f 'PASS,.' --no-update -O u | \
        ## remove extraneous INFO and FORMAT fields

        bcftools annotate -x ^FORMAT/FT,FORMAT/GT -O b -o filtered.bcf ##  ^INFO/AC,INFO/AF,INFO/AN,

        echo_date "VAT: normalize, left align and split multi allelic sites to new lines, remove duplicate lines"
        bcftools norm --threads 4 -m- --check-ref w -f Homo_sapiens_assembly38.fasta filtered.bcf -O b -o normalized.bcf
        rm filtered.bcf

        echo_date "VAT: filter out spanning deletions and variants with an AC of 0, respect the FT flag"
        bcftools view --threads 4 -e 'ALT[0]="*" || AC=0' --no-update normalized.bcf -O u | \
        ## ensure that we respect the FT tag
        bcftools filter --threads 4 -i "FORMAT/FT='PASS,.'" --set-GTs . -O b -o normalized.filtered.bcf

        ## clean up unneeded file
        rm normalized.bcf

        echo_date "VAT: detecting and removing duplicate rows"

        ## During normalization, sometimes duplicate variants appear but with different calculations. This seems to be a bug in bcftools. For now we are dropping all duplicate variants
        ## to locate the duplicates, we first make a file of just the first 5 columns
        grep -v '^#' normalized.filtered.bcf | cut -f1,2,4,5 | sort | uniq -d > duplicates.tsv

        # If there ARE dupes to remove
        if [ -s duplicates.tsv ]; then
            ## remove those rows (that match up to the first 5 cols)
            bcftools view --threads 4 normalized.filtered.bcf | grep -v -wFf duplicates.tsv | bcftools view --threads 4 -O b -o deduplicated.bcf
        else
            # There are no duplicates to remove
            cp normalized.filtered.bcf deduplicated.bcf
        fi
        rm normalized.filtered.bcf

        ## add duplicates to the file that's tracking dropped variants
        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv ## clean up unneeded file

        echo_date "VAT: calculate annotations for all subpopulations"
        ## AC_het,AC_hom and AC_Hemi are used to calculate the participant count
        bcftools plugin fill-tags --threads 4 -- deduplicated.bcf -S ~{subpopulation_sample_list} -t AC,AF,AN,AC_het,AC_hom,AC_Hemi | bcftools query -f \
        '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\t%AC_Hom\t%AC_Het\t%AC_Hemi\t%AC_afr\t%AN_afr\t%AF_afr\t%AC_Hom_afr\t%AC_Het_afr\t%AC_Hemi_afr\t%AC_amr\t%AN_amr\t%AF_amr\t%AC_Hom_amr\t%AC_Het_amr\t%AC_Hemi_amr\t%AC_eas\t%AN_eas\t%AF_eas\t%AC_Hom_eas\t%AC_Het_eas\t%AC_Hemi_eas\t%AC_eur\t%AN_eur\t%AF_eur\t%AC_Hom_eur\t%AC_Het_eur\t%AC_Hemi_eur\t%AC_mid\t%AN_mid\t%AF_mid\t%AC_Hom_mid\t%AC_Het_mid\t%AC_Hemi_mid\t%AC_oth\t%AN_oth\t%AF_oth\t%AC_Hom_oth\t%AC_Het_oth\t%AC_Hemi_oth\t%AC_sas\t%AN_sas\t%AF_sas\t%AC_Hom_sas\t%AC_Het_sas\t%AC_Hemi_sas\n' \
        >> ~{custom_annotations_file_name}

        ## for validation of the pipeline
        wc -l ~{custom_annotations_file_name} | awk '{print $1 -7}'  > count.txt

        echo_date "VAT: produce sites-only VCf"

        ## compress the vcf and index it, make it sites-only for the next step
        bcftools view --threads 4 --no-update --drop-genotypes deduplicated.bcf -O z -o ~{normalized_vcf_compressed}
        ## if we can spare the IO and want to pass a smaller file we can also drop the info field w bcftools annotate -x INFO
        bcftools index --tbi ~{normalized_vcf_compressed}

        echo_date "VAT: finished"
    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_05_02"
        maxRetries: 3
        memory: "16 GB"
        preemptible: 3
        cpu: "4"
        disks: "local-disk 500 HDD"
    }
    # ------------------------------------------------
    # Outputs:
    output {
        File annotations_file = "~{custom_annotations_file_name}"
        Int count_variants = read_int("count.txt")
        File track_dropped = "track_dropped.tsv"
        File output_vcf = "~{normalized_vcf_compressed}"
        File output_vcf_index = "~{normalized_vcf_indexed}"
    }
}

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

        echo "Extracting annotation data sources tar/gzip file..."
        mkdir datasources_dir
        tar zxvf ~{nirvana_data_tar} -C datasources_dir  --strip-components 2
        DATA_SOURCES_FOLDER="$PWD/datasources_dir"


        # =======================================
        echo "Creating custom annotations"
        mkdir customannotations_dir
        CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

        # Add AC/AN/AF as custom annotations
        ## use --skip-ref once you are on a version of nirvana later than 3.14 (once they have created a docker image for it)
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
        disks: "local-disk 250 HDD"
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

    ## note: these temp files do not currently get cleaned up as some of them may be helpful for recovery.

    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
        fi

        # for debugging purposes only
        gsutil cp ~{annotation_json} '~{output_annotations_gcp_path}'

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_variant_annotation_table.py \
          --annotated_json ~{annotation_json} \
          --output_vt_json ~{output_vt_json} \
          --output_genes_json ~{output_genes_json}

        gsutil cp ~{output_vt_json} '~{output_vt_gcp_path}'
        gsutil cp ~{output_genes_json} '~{output_genes_gcp_path}'

     >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_05_02"
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