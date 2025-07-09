version 1.0

import "GvsUtils.wdl" as Utils

struct ReferenceFiles {
  String fasta_bgz
  String fasta_index
  String sequence_dictionary
  String contig_mapping
  String vet_weight_bins
  String weighted_bed
}

workflow PrepareReferenceFiles {
  input {
    File reference_fasta
    String output_gcs_dir
    String project_id
    String dataset_name
  }

  call Utils.GetToolVersions {}

  call GenerateBgzSequenceDictionaryAndIndex {
    input:
      reference_fasta = reference_fasta,
      output_gcs_dir = output_gcs_dir,
      gatk_docker = GetToolVersions.gatk_docker,
  }

  call GenerateContigMapping {
    input:
      sequence_dictionary = read_json(GenerateBgzSequenceDictionaryAndIndex.reference_files_json).sequence_dictionary,
      in_reference_json = GenerateBgzSequenceDictionaryAndIndex.reference_files_json,
      output_gcs_dir = output_gcs_dir,
      variants_docker = GetToolVersions.variants_docker,
  }

  call CreateWeightedBedFile {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      reference_dictionary = read_json(GenerateBgzSequenceDictionaryAndIndex.reference_files_json).sequence_dictionary,
      contig_mapping = read_json(GenerateContigMapping.reference_files_json).contig_mapping,
      in_reference_json = GenerateContigMapping.reference_files_json,
      output_gcs_dir = output_gcs_dir,
      variants_docker = GetToolVersions.variants_docker,
  }

  output {
    ReferenceFiles reference_files = read_json(CreateWeightedBedFile.updated_reference_files_json)
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
    String output_gcs_dir
    String variants_docker
    File reference_dictionary
    File contig_mapping
    File in_reference_json
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

    echo "Vet weight bins generated."
    ls -l vet_weight_bins.tsv
    echo "That's the file"

    python /app/weighted_bed_from_vet_bins.py \
      --input-bin-data vet_weight_bins.tsv \
      --reference-dict ~{reference_dictionary} \
      --contig-mapping ~{contig_mapping} \
      > weighted.bed

    # Ensure no trailing slash on the GCS output directory
    output_gcs_dir="~{output_gcs_dir}"
    output_gcs_dir="${output_gcs_dir%/}"

    base=$(basename ~{reference_dictionary} .dict)
    cloud_weight_bins_file="${output_gcs_dir}/${base}.vet_weight_bins.tsv"
    gcloud storage cp vet_weight_bins.tsv ${cloud_weight_bins_file}

    cloud_weighted_bed_file="${output_gcs_dir}/${base}.weighted.bed"
    gcloud storage cp weighted.bed ${cloud_weighted_bed_file}

    jq " . += {
      \"vet_weight_bins\": \"$cloud_weight_bins_file\",
      \"weighted_bed\": \"$cloud_weighted_bed_file\"
    } " ~{in_reference_json} > updated_reference_files.json

  >>>

  runtime {
    docker: variants_docker
    memory: "4 GB"
    disks: "local-disk 500 HDD"
    cpu: 1
  }

  output {
    File updated_reference_files_json = "updated_reference_files.json"
  }
}
