version 1.0

workflow GvsExtractAvroFilesForHail {
    input {
        String project_id
        String dataset
        String filter_set_name
        String gcs_temporary_path
    }

    call OutputPath { input: go = true }

    call ExtractFromNonSuperpartitionedTables {
        input:
            project_id = project_id,
            dataset = dataset,
            filter_set_name = filter_set_name,
            avro_sibling = OutputPath.out
    }

    call CountSamples {
        input:
            project_id = project_id,
            dataset = dataset
    }

    # Superpartitions have max size 4000. The inner '- 1' is so the 4000th (and multiples of 4000) sample lands in the
    # appropriate partition, the outer '+ 1' is to iterate over the correct number of partitions.
    scatter (i in range(((CountSamples.num_samples - 1) / 4000) + 1)) {
        call ExtractFromSuperpartitionedTables {
            input:
                project_id = project_id,
                dataset = dataset,
                filter_set_name = filter_set_name,
                avro_sibling = OutputPath.out,
                table_index = i + 1 #  'i' is 0-based so add 1.
        }
    }

    call GenerateHailScript {
        input:
            go_non_superpartitioned = ExtractFromNonSuperpartitionedTables.done,
            go_superpartitioned = ExtractFromSuperpartitionedTables.done,
            avro_prefix = ExtractFromNonSuperpartitionedTables.output_prefix,
            gcs_temporary_path = gcs_temporary_path
    }
    output {
        String output_prefix = ExtractFromNonSuperpartitionedTables.output_prefix
    }
}


task OutputPath {
    meta {
        description: "Does nothing but produce the cloud path to its stdout."
        # Always make a new output path, otherwise every invocation would clobber the original.
        volatile: true
    }
    input {
        Boolean go = true
    }
    command <<<
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


task CountSamples {
    meta {
        description: "Counts the number of samples in the sample_info table efficiently."
        # There's no last modified check so always count the samples.
        volatile: true
    }
    input {
        String project_id
        String dataset
    }
    command <<<
        python3 <<FIN

        from google.cloud import bigquery

        client = bigquery.Client(project="~{project_id}")
        sample_info_table_id = f'~{project_id}.~{dataset}.sample_info'
        sample_info_table = client.get_table(sample_info_table_id)

        print(str(sample_info_table.num_rows))

        FIN
    >>>

    output {
        Int num_samples = read_int(stdout())
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2022_08_22"
    }
}


task ExtractFromNonSuperpartitionedTables {
    meta {
        description: "Extracts from the non-superpartitioned tables: sample_info, filter_set_info, filter_set_sites"
        # Always run, no last modified check.
        volatile: true
    }
    input {
        String project_id
        String dataset
        String filter_set_name
        String avro_sibling
    }
    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        echo $avro_prefix > "avro_prefix.out"

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/sample_mapping/sample_mapping_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT sample_id, sample_name, '40',
            'gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list' as intervals_file
            FROM \`~{project_id}.~{dataset}.sample_info\`
            ORDER BY sample_id
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vqsr_filtering_data/vqsr_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, type as model, ref, alt, vqslod, yng_status
            FROM \`~{project_id}.~{dataset}.filter_set_info\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/site_filtering_data/site_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, filters
            FROM \`~{project_id}.~{dataset}.filter_set_sites\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vqsr_tranche_data/vqsr_tranche_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT model, truth_sensitivity, min_vqslod, filter_name
            FROM \`~{project_id}.~{dataset}.filter_set_tranches\`
            WHERE filter_set_name = '~{filter_set_name}'
        "
    >>>

    output {
        Boolean done = true
        String output_prefix = read_string("avro_prefix.out")
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
    }
}


task ExtractFromSuperpartitionedTables {
    meta {
        description: "Extracts from the superpartitioned tables: vet_<table index>, ref_ranges_<table index>"
        # Always run, no last modified check.
        volatile: true
    }
    input {
        String project_id
        String dataset
        String filter_set_name
        String avro_sibling
        Int table_index
    }
    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
        table_index: "1-based index for the superpartitioned ref_ranges and vet tables to be extracted by this call."
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        str_table_index=$(printf "%03d" ~{table_index})

        # These bq exports error out if there are any objects at the sibling level to where output files would be written
        # so an extra layer of `vet_${str_table_index}` is inserted here.
        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vets/vet_${str_table_index}/vet_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, sample_id, ref, REPLACE(alt,',<NON_REF>','') alt, call_GT as GT, call_AD as AD, call_GQ as GQ, cast(SPLIT(call_pl,',')[OFFSET(0)] as int64) as RGQ
            FROM \`~{project_id}.~{dataset}.vet_${str_table_index}\`
            ORDER BY location
        "

        bq query --nouse_legacy_sql --project_id=~{project_id} "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/refs/ref_ranges_${str_table_index}/ref_ranges_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, sample_id, length, state
            FROM \`~{project_id}.~{dataset}.ref_ranges_${str_table_index}\`
            ORDER BY location
        "
    >>>

    output {
        Boolean done = true
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
    }
}

task GenerateHailScript {
    input {
        String avro_prefix
        String gcs_temporary_path
        Boolean go_non_superpartitioned
        Array[Boolean] go_superpartitioned
    }
    meta {
        # Do not cache, this doesn't know if the "tree" under `avro_prefix` has changed.
        volatile: true
    }
    parameter_meta {
        gcs_temporary_path: "Path to network-visible temporary directory/bucket for intermediate file storage"
        go_non_superpartitioned: "Sync on completion of non-superpartitioned extract"
        go_superpartitioned: "Sync on completion of all superpartitioned extract shards"
    }

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        vds_output_path="$(dirname ~{avro_prefix})/gvs_export.vds"
        echo $vds_output_path > vds_output_path.txt

        vcf_output_path="$(dirname ~{avro_prefix})/gvs_export.vcf"
        echo $vcf_output_path > vcf_output_path.txt

        gsutil ls -r '~{avro_prefix}' > avro_listing.txt

        python3 /app/generate_hail_gvs_import.py \
            --avro_prefix '~{avro_prefix}' \
            --avro_listing_file avro_listing.txt \
            --vds_output_path "${vds_output_path}" \
            --vcf_output_path "${vcf_output_path}" \
            --gcs_temporary_path ~{gcs_temporary_path} > hail_script.py
    >>>

    output {
        Boolean done = true
        String vds_output_path = read_string('vds_output_path.txt')
        String vcf_output_path = read_string('vcf_output_path.txt')
        String hail_script = read_string('hail_script.py')
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:vs_605_hail_codegen"
    }
}
