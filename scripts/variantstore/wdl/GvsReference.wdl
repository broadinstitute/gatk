version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsReference {
    input {
        String reference_name
        File? custom_reference

        # If `git_branch_or_tag` is not specified by a caller (i.e. integration tests), default to the current branch or
        # tag as specified in `GvsUtils.GetToolVersions`.
        String? git_branch_or_tag
        # Potentially specified by a calling integration WDL.
        String? git_hash

        String? basic_docker
        String? variants_docker
        String? gatk_docker

        String? workspace_bucket
        String? submission_id
    }

    if ((reference_name != "hg38") || defined(custom_reference)) {
        if (!defined(git_hash) ||
            !defined(basic_docker) || !defined(variants_docker) || !defined(gatk_docker) ||
            !defined(workspace_bucket) || !defined(submission_id)) {
            call Utils.GetToolVersions {
                input:
                    git_branch_or_tag = git_branch_or_tag,
            }
        }

        String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
        String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
        String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
        String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
        String effective_submission_id = select_first([submission_id, GetToolVersions.submission_id])

        if ((reference_name != "hg38") && (reference_name != "CUSTOM")) {
            call Utils.TerminateWorkflow as UnknownReferenceName {
                input:
                    basic_docker = effective_basic_docker,
                    message = "Unrecognized value 'reference_name' defined. Must be either 'hg38' or 'CUSTOM`. Exiting."
            }
        }

        if ((reference_name == "hg38") && defined(custom_reference)) {
            call Utils.TerminateWorkflow as Hg38CustomReferenceDefined {
                input:
                    basic_docker = effective_basic_docker,
                    message = "If reference_name is set to 'hg38' then 'custom_reference' should not be defined. Exiting."
            }
        }

        if ((reference_name == "CUSTOM") && !defined(custom_reference)) {
            call Utils.TerminateWorkflow as CustomMustHaveCustReferenceDefined {
                input:
                    basic_docker = effective_basic_docker,
                    message = "If 'reference_name' is set to 'CUSTOM', then 'custom_reference' must be defined. Exiting."
            }
        }

        # TODO - I'm not wild about storing these things on thge cloud. (e.g. index, dictionary) - maybe make it optional?

        if (defined(custom_reference)) {
            call GenerateBgzSequenceDictionaryAndIndex {
                input:
                    reference_fasta = select_first([custom_reference]),
                    output_gcs_dir = effective_workspace_bucket + "/submissions/" + effective_submission_id,
                    gatk_docker = effective_gatk_docker,
            }

            call GenerateContigMapping {
                input:
                    sequence_dictionary = GenerateBgzSequenceDictionaryAndIndex.sequence_dictionary,
                    variants_docker = effective_variants_docker,
            }
        }
    }

    output {
        Boolean is_custom_reference = if (defined(custom_reference)) then true else false
        String reference_version = if (defined(custom_reference)) then "CUSTOM" else "38"
        File reference_fasta = select_first([custom_reference, "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"])
        File reference_fasta_index = select_first([GenerateBgzSequenceDictionaryAndIndex.fasta_index, "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"])
        File reference_dict = select_first([GenerateBgzSequenceDictionaryAndIndex.sequence_dictionary, "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"])
        File? custom_contig_mapping_file = GenerateContigMapping.custom_contig_mapping_file
        File wgs_calling_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
        File exome_calling_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list"
    }
}

task GenerateBgzSequenceDictionaryAndIndex {
    input {
        File reference_fasta
        String output_gcs_dir
        String gatk_docker
    }
    parameter_meta {
        reference_fasta: {
             help: "Reference FASTA file, can be compressed with bgzip, gzip or uncompressed."
         }
    }

    # The base name of the fasta file without the (possible) .bgz or .gz extension
    String base_filename = sub(basename(reference_fasta), "\\.bgz$|\\.gz$", "")

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "Hello!"
        echo "Reference FASTA file: ~{reference_fasta}"
        echo "Base filename from that: ~{base_filename}"

        mkdir output

        if [[ ~{reference_fasta} == *.bgz ]]
        then
            echo "Input reference FASTA file is already compressed with bgzip."
            mv ~{reference_fasta} output/
        elif [[ ~{reference_fasta} == *.gz ]]
        then
            # A test to see if the input reference FASTA file is compressed with bgzip.
            if bgzip --reindex ~{reference_fasta} >/dev/null 2>&1
            then
                echo "Input reference FASTA file is compressed with bgzip."
                base=$(basename ~{reference_fasta} .gz)
                mv ~{reference_fasta} output/${base}.bgz
            else
                echo "Input reference FASTA file is not compressed with bgzip, uncompressing with gzip and recompressing with bgzip."
                base=$(basename ~{reference_fasta} .gz)
                gzip --decompress --stdout ~{reference_fasta} > ${base}
                bgzip ${base} --stdout > output/${base}.bgz
            fi
        else
            echo "Input reference FASTA file is uncompressed, compressing with bgzip."
            base=$(basename ~{reference_fasta})
            bgzip ~{reference_fasta} --stdout > output/${base}.bgz
        fi

        # Ensure no trailing slash on the GCS output directory
        output_gcs_dir="~{output_gcs_dir}"
        output_gcs_dir="${output_gcs_dir%/}"

        base=$(basename output/*.bgz)

        # Generate sequence dictionary using samtools. Make sure to specify the URI corresponding to where the bgzipped
        # FASTA will be stored in the GCS output directory.
        samtools dict --uri ${output_gcs_dir}/${base} output/*.bgz > output/$(basename output/*.bgz).dict

        # Generate FASTA index using samtools
        samtools faidx output/*.bgz

        echo "hello"
        ls -l output/
        echo "A"
        ls -l output/~{base_filename}.bgz
        echo "B"
        ls -l output/~{base_filename}.bgz.fai
        echo "C"
        ls -l output/~{base_filename}.bgz.dict

        gcloud storage cp output/*.bgz output/*.fai output/*.dict ${output_gcs_dir}/

        echo "{
            \"fasta_bgz\": \"$output_gcs_dir/$(basename output/*.bgz)\",
            \"fasta_index\": \"$output_gcs_dir/$(basename output/*.fai)\",
            \"sequence_dictionary\": \"$output_gcs_dir/$(basename output/*.dict)\"
            }" > reference_files.json
    >>>

    runtime {
        docker: gatk_docker
        memory: "4 GB"
        disks: "local-disk " + ceil(5 * size(reference_fasta, "GB")) + " HDD"
        cpu: 1
    }

    output {
        File fasta_bgz = "output/" + base_filename + ".bgz"
        File fasta_index = "output/" + base_filename + ".bgz.fai"
        File sequence_dictionary = "output/" + base_filename + ".bgz.dict"
        File reference_files_json = "reference_files.json"
    }
}

task GenerateContigMapping {
    input {
        File sequence_dictionary
        String variants_docker
    }

    String contig_mapping_filename = basename(sequence_dictionary, ".dict") + ".contig_mapping.tsv"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "Hello!"
        echo ~{contig_mapping_filename}
        python3 /app/generate_custom_reference_mappings.py \
            ~{sequence_dictionary} > ~{contig_mapping_filename}

    >>>

    runtime {
        docker: variants_docker
        memory: "4 GB"
        disks: "local-disk 500 HDD"
        cpu: 1
    }

    output {
        File custom_contig_mapping_file = contig_mapping_filename
    }
}

task CreateWeightedBedFile {
    input {
        Boolean go = true
        String project_id
        String dataset_name
        String variants_docker
        File reference_dictionary
        File contig_mapping
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Check if the vet_001 table exists
        if ! bq --apilog=false --project_id=~{project_id} show ~{dataset_name}.vet_001 &> /dev/null; then
            echo "ERROR: Table ~{dataset_name}.vet_001 does not exist. Please ensure the table exists before running this task."
            exit 1
        fi

        echo "Generating vet weight bins from BigQuery dataset ~{dataset_name} in project ~{project_id}."

        bq --apilog=false --project_id=~{project_id} query --max_rows 1000000000 --format=csv --use_legacy_sql=false '
            SELECT CAST(TRUNC(location / 1000) * 1000 AS INT64) bin, count(*) entries
            FROM ~{dataset_name}.vet_001
            GROUP BY bin ORDER BY bin' > vet_weight_bins.tsv

        echo "Converting vet weight bins to BED format."
        python /app/convert_weight_file_CSV_to_bed.py \
            --bin-size 1 \
            --mapping-file ~{contig_mapping} \
            --infile vet_weight_bins.tsv \
            --outfile vet_weight_bins.bed

        echo "Closing gaps in BED file."
        python /app/close_bed_file_gaps.py \
            --bin-size 1 \
            --infile vet_weight_bins.bed \
            --outfile vet_weight_bins_closed.bed \
            --interval-list ~{reference_dictionary}

        echo "Successfully created weighted.bed file 'vet_weight_bins_closed.bed'."

    >>>

    runtime {
        docker: variants_docker
        memory: "4 GB"
        disks: "local-disk 500 HDD"
        cpu: 1
    }

    output {
        File weighted_bed_file = "vet_weight_bins_closed.bed"
        File vet_weight_bins_file = "vet_weight_bins.tsv"
    }
}