version 1.0

workflow GvsCreateVATSubpopulationCalculations {
   input {
        File inputFileofFileNames
        File inputFileofIndexFileNames
        String output_path

        String? service_account_json_path
        File? gatk_override

        File AnAcAf_annotations_template
        File ancestry_file
        File reference
   }

    ## Scatter across the shards from the GVS jointVCF
    scatter(i in range(length(MakeSubpopulationFiles.input_vcfs)) ) {
        ## Create a sites-only VCF from the original GVS jointVCF
        ## Calculate AC/AN/AF for subpopulations and extract them for custom annotations
        ## To prevent premature failures from this brittle step, the output will be saved to GCP
        String input_vcf_name = basename(MakeSubpopulationFiles.input_vcfs[i], ".vcf.gz")
        call ExtractAnAcAfFromVCF {
            input:
              input_vcf = MakeSubpopulationFiles.input_vcfs[i],
              input_vcf_index = MakeSubpopulationFiles.input_vcf_indices[i],
              service_account_json_path = service_account_json_path,
              subpopulation_sample_list = MakeSubpopulationFiles.ancestry_mapping_list,
              custom_annotations_template = AnAcAf_annotations_template,
              ref = reference,
              output_path = output_path
        }
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
    String normalized_vcf = "normalized.vcf"
    String normalized_vcf_compressed = "normalized.vcf.gz"
    String normalized_vcf_indexed = "normalized.vcf.gz.tbi"

    String output_directory_name = output_path + 'med/'


    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting and calculating the an/ac/af & sc by subpopulation into a tsv
    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
          gsutil cp ~{service_account_json_path} local.service_account.json
          export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
          gcloud auth activate-service-account --key-file=local.service_account.json

        fi

        cp ~{custom_annotations_template} ~{custom_annotations_file_name}

        gsutil cp ~{input_vcf} ~{local_input_vcf}
        gsutil cp ~{input_vcf_index} ~{local_input_vcf_index}
        gsutil cp ~{ref} Homo_sapiens_assembly38.fasta

        ## compare the ancestry sample list with the vcf sample list
        # TODO throw an error (but dont fail the job) if there are samples that are in one, but not the other. Throw two different errors.
        # Currently commented out to save on io with the AoU beta VAT creation
        # awk '{print $1}' ~{subpopulation_sample_list} | tail -n +2 | sort -u > collected_subpopulation_samples.txt
        # bcftools query --list-samples ~{local_input_vcf} | sort -u > collected_samples.txt
        # diff collected_subpopulation_samples.txt collected_samples.txt

        # expected_subpopulations = [
        # "afr",
        # "amr",
        # "eas",
        # "eur",
        # "mid",
        # "oth",
        # "sas"
        #]

        ## track the dropped variants with +50 alt alleles or N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        bcftools view -i 'N_ALT>50 || REF~"N"' ~{local_input_vcf} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        wc -l track_dropped.tsv

        ## filter out sites with too many alt alleles
        bcftools view -e 'N_ALT>50 || REF~"N"' --no-update ~{local_input_vcf} -Ou | \
        ## filter out the non-passing sites
        bcftools view  -f 'PASS,.' --no-update -Oz -o filtered.vcf.gz
        ## normalize, left align and split multi allelic sites to new lines, remove duplicate lines
        bcftools norm -m- --check-ref w -f Homo_sapiens_assembly38.fasta filtered.vcf.gz -Oz -o normalized.vcf.gz
        rm ~{local_input_vcf}
        ## filter out spanning deletions and variants with an AC of 0
        bcftools view -e 'ALT[0]="*" || AC=0' --no-update normalized.vcf.gz -Ou | \
        ## ensure that we respect the FT tag
        bcftools filter -i "FORMAT/FT='PASS,.'" --set-GTs . -Oz -o ~{normalized_vcf}

        ## clean up unneeded file
        rm normalized.vcf.gz

        ## During normalization, sometimes duplicate variamts appear but with different calculations. This seems to be a bug in bcftools. For now we are dropping all duplicate variants
        ## The way in which this is done is a bit hamfisted and should be optimized in the future.
        ## to locate the duplicates, we first make a file of just the first 5 columns
        bcftools query ~{normalized_vcf} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > check_duplicates.tsv
        ## check it for duplicates and put them in a new file
        sort check_duplicates.tsv | uniq -d | cut -f1,2,3,4,5  > duplicates.tsv
        rm check_duplicates.tsv ## clean up
        ## remove those rows (that match up to the first 5 cols)
        zgrep -v -wFf duplicates.tsv ~{normalized_vcf} > deduplicated.vcf.gz
        rm ~{normalized_vcf} ## clean up

        ## add duplicates to the file that's tracking dropped variants
        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv ## clean up unneeded file

        ## calculate annotations for all subpopulations
        ## AC_het,AC_hom and AC_Hemi are used to calculate the participant count
        bcftools plugin fill-tags  -- deduplicated.vcf.gz -S ~{subpopulation_sample_list} -t AC,AF,AN,AC_het,AC_hom,AC_Hemi | bcftools query -f \
        '%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\t%AC_Hom\t%AC_Het\t%AC_Hemi\t%AC_afr\t%AN_afr\t%AF_afr\t%AC_Hom_afr\t%AC_Het_afr\t%AC_Hemi_afr\t%AC_amr\t%AN_amr\t%AF_amr\t%AC_Hom_amr\t%AC_Het_amr\t%AC_Hemi_amr\t%AC_eas\t%AN_eas\t%AF_eas\t%AC_Hom_eas\t%AC_Het_eas\t%AC_Hemi_eas\t%AC_eur\t%AN_eur\t%AF_eur\t%AC_Hom_eur\t%AC_Het_eur\t%AC_Hemi_eur\t%AC_mid\t%AN_mid\t%AF_mid\t%AC_Hom_mid\t%AC_Het_mid\t%AC_Hemi_mid\t%AC_oth\t%AN_oth\t%AF_oth\t%AC_Hom_oth\t%AC_Het_oth\t%AC_Hemi_oth\t%AC_sas\t%AN_sas\t%AF_sas\t%AC_Hom_sas\t%AC_Het_sas\t%AC_Hemi_sas\n' \
        >> ~{custom_annotations_file_name}

        ## for validation of the pipeline
        wc -l ~{custom_annotations_file_name} | awk '{print $1 -7}'  > count.txt

        ## compress the vcf and index it, make it sites-only for the next step
        bcftools view --no-update --drop-genotypes deduplicated.vcf.gz -Oz -o ~{normalized_vcf_compressed}
        ## if we can spare the IO and want to pass a smaller file we can also drop the info field w bcftools annotate -x INFO
        bcftools index --tbi  ~{normalized_vcf_compressed}


        ## To prevent premature failures from this brittle step, the output will be saved to GCP
        ## There are 5 (!!!) files that will need to get saved---so let's make a directory for them
        # String output_directory_name
        ## and the five files will be:
        # custom_annotations_file_name, count.txt, track_dropped.tsv, normalized_vcf_compressed, normalized_vcf_indexed
        ## and they will be loaded with gsutil commands
        # gsutil cp ~{custom_annotations_file_name} '~{output_directory_name}/~{custom_annotations_file_name}'
        # gsutil cp count.txt '~{output_directory_name}/count.txt'
        # gsutil cp track_dropped.tsv '~{output_directory_name}/track_dropped.tsv}'
        # gsutil cp ~{normalized_vcf_compressed} '~{output_directory_name}/~{normalized_vcf_compressed}'
        # gsutil cp ~{normalized_vcf_indexed} '~{output_directory_name}/~{normalized_vcf_indexed}'

    >>>
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_20211101"
        maxRetries: 3
        memory: "32 GB"
        preemptible: 3
        cpu: "4"
        disks: "local-disk 500 SSD"
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
