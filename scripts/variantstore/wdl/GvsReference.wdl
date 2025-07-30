version 1.0

import "GvsUtils.wdl" as Utils

import "../structs/Reference.wdl" as Reference

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
                message = "Unrecognized `reference_name`: " + reference_name + ". Exiting."
        }
    }

    if ((reference_name == "hg38") && defined(custom_reference)) {
        call Utils.TerminateWorkflow as Hg38CustomReferenceDefined {
            input:
                basic_docker = effective_basic_docker,
                message = "If `reference_name` is set to `hg38`, ` then `custom_reference` should not be defined. Exiting."
        }
    }

    if ((reference_name == "CUSTOM") && !defined(custom_reference)) {
        call Utils.TerminateWorkflow as CustomMustHaveCustReferenceDefined {
            input:
                basic_docker = effective_basic_docker,
                message = "If `reference_name` is set to `CUSTMO`, ` then `custom_reference` must be defined. Exiting."
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

        File? custom_reference_index = read_json(GenerateBgzSequenceDictionaryAndIndex.reference_files_json).reference_index
        File? custom_sequence_dictionary = read_json(GenerateBgzSequenceDictionaryAndIndex.reference_files_json).sequence_dictionary

        call GenerateContigMapping {
            input:
                sequence_dictionary = select_first([custom_sequence_dictionary]),
                in_reference_json = GenerateBgzSequenceDictionaryAndIndex.reference_files_json,
                output_gcs_dir = effective_workspace_bucket + "/submissions/" + effective_submission_id,
                variants_docker = effective_variants_docker,
        }

        File? custom_contig_mapping = read_json(GenerateContigMapping.reference_files_json).contig_mapping
    }

    Reference hg38_reference = {
        "reference_version": "38",
        "reference_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
        "reference_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "reference_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
        "wgs_calling_interval_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list",
        "exome_calling_interval_list": "gs://gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list",
    }

    output {
        Reference reference = hg38_reference
        Boolean is_custom_reference = if (defined(custom_reference)) then true else false
        String reference_version = if (defined(custom_reference)) then "CUSTOM" else "38"
        File reference_fasta = select_first([custom_reference, "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"])
        File reference_fasta_index = select_first([custom_reference_index, "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"])
        File reference_dict = select_first([custom_sequence_dictionary, "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"])
        File? custom_contig_mapping_file = custom_contig_mapping
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

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

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
        File reference_files_json = "reference_files.json"
    }
}

task GenerateContigMapping {
    input {
        File sequence_dictionary
        File in_reference_json
        String output_gcs_dir
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        base=$(basename ~{sequence_dictionary} .dict)

        python3 /app/generate_custom_reference_mappings.py \
        ~{sequence_dictionary} > ${base}.contig_mapping.tsv

        # Ensure no trailing slash on the GCS output directory
        output_gcs_dir="~{output_gcs_dir}"
        output_gcs_dir="${output_gcs_dir%/}"

        gcloud storage cp ${base}.contig_mapping.tsv ${output_gcs_dir}/

        # Update the reference files JSON with the contig mapping file
        cloud_contig_mapping_file="${output_gcs_dir}/$(basename ${base}.contig_mapping.tsv)"

        jq " . += {
        \"contig_mapping\": \"$cloud_contig_mapping_file\"
        } " ~{in_reference_json} > updated_reference_files.json

    >>>

    runtime {
        docker: variants_docker
        memory: "4 GB"
        disks: "local-disk 500 HDD"
        cpu: 1
    }

    output {
        File reference_files_json = "updated_reference_files.json"
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