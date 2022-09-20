version 1.0

workflow GvsExtractAvroFilesForHail {
    input {
        String project_id
        String dataset
        String filter_set_name
        Int scatter_width = 10
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

    Int num_samples = CountSamples.num_samples
    # First superpartition contains samples 1 to 4000, second 4001 to 8000 etc; add one to quotient unless exactly 4000.
    Int num_superpartitions = if (num_samples % 4000 == 0) then num_samples / 4000 else (num_samples / 4000 + 1)

    scatter (i in range(scatter_width)) {
        call ExtractFromSuperpartitionedTables {
            input:
                project_id = project_id,
                dataset = dataset,
                filter_set_name = filter_set_name,
                avro_sibling = OutputPath.out,
                num_superpartitions = num_superpartitions,
                shard_index = i,
                num_shards = scatter_width
        }
    }

    call GenerateHailScripts {
        input:
            go_non_superpartitioned = ExtractFromNonSuperpartitionedTables.done,
            go_superpartitioned = ExtractFromSuperpartitionedTables.done,
            avro_prefix = ExtractFromNonSuperpartitionedTables.output_prefix,
    }
    output {
        File hail_gvs_import_script = GenerateHailScripts.hail_gvs_import_script
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
        # Not dealing with caching for now as that would introduce a lot of complexity.
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
        # Not dealing with caching for now as that would introduce a lot of complexity.
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
            WHERE withdrawn IS NULL
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
        # Not dealing with caching for now as that would introduce a lot of complexity.
        volatile: true
    }
    input {
        String project_id
        String dataset
        String filter_set_name
        String avro_sibling
        Int num_superpartitions
        Int shard_index
        Int num_shards
    }
    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
        num_superpartitions: "Total number of superpartitions requiring extraact"
        shard_index: "0-based index of this superpartition extract shard"
        num_shards: "Count of all superpartition extract shards"
    }

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        avro_prefix="$(dirname ~{avro_sibling})/avro"

        for superpartition in $(seq ~{shard_index + 1} ~{num_shards} ~{num_superpartitions})
        do
            str_table_index=$(printf "%03d" $superpartition)

            # These bq exports error out if there are any objects at the sibling level to where output files would be written
            # so an extra layer of `vet_${str_table_index}` is inserted here.
            bq query --nouse_legacy_sql --project_id=~{project_id} "
                EXPORT DATA OPTIONS(
                uri='${avro_prefix}/vets/vet_${str_table_index}/vet_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
                SELECT location, v.sample_id, ref, REPLACE(alt,',<NON_REF>','') alt, call_GT as GT, call_AD as AD, call_GQ as GQ, cast(SPLIT(call_pl,',')[OFFSET(0)] as int64) as RGQ
                FROM \`~{project_id}.~{dataset}.vet_${str_table_index}\` v
                INNER JOIN \`~{project_id}.~{dataset}.sample_info\` s ON s.sample_id = v.sample_id
                WHERE withdrawn IS NULL
                ORDER BY location
            "

            bq query --nouse_legacy_sql --project_id=~{project_id} "
                EXPORT DATA OPTIONS(
                uri='${avro_prefix}/refs/ref_ranges_${str_table_index}/ref_ranges_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
                SELECT location, r.sample_id, length, state
                FROM \`~{project_id}.~{dataset}.ref_ranges_${str_table_index}\` r
                INNER JOIN \`~{project_id}.~{dataset}.sample_info\` s ON s.sample_id = r.sample_id
                WHERE withdrawn IS NULL
                ORDER BY location
            "
        done
    >>>

    output {
        Boolean done = true
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:398.0.0"
    }
}

task GenerateHailScripts {
    input {
        String avro_prefix
        Boolean go_non_superpartitioned
        Array[Boolean] go_superpartitioned
    }
    meta {
        # Do not cache, this doesn't know if the "tree" under `avro_prefix` has changed.
        volatile: true
    }
    parameter_meta {
        go_non_superpartitioned: "Sync on completion of non-superpartitioned extract"
        go_superpartitioned: "Sync on completion of all superpartitioned extract shards"
    }

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        # A bit of randomness to not clobber outputs if run multiple times.
        rand=$(openssl rand -hex 4)

        # The write prefix will be a sibling to the Avro "directory" that embeds the current date and some randomness.
        write_prefix="$(dirname ~{avro_prefix})/$(date -Idate)-${rand}"

        vds_output_path="${dir}/gvs_export.vds"
        echo $vds_output_path > vds_output_path.txt

        # vcf_output_path="${dir}/gvs_export.vcf"
        # echo $vcf_output_path > vcf_output_path.txt

        # sites_only_vcf_output_path="${dir}/gvs_sites_only.vcf"
        # echo $sites_only_vcf_output_path > sites_only_vcf_output_path.txt

        cat /app/hail_gvs_import.tmpl.py |
            sed "s/@AVRO_PREFIX@/~{avro_prefix}" |
            sed "s/@WRITE_PREFIX@/$write_prefix/" >
            hail_gvs_import.py

    >>>

    output {
        Boolean done = true
        String vds_output_path = read_string('vds_output_path.txt')
        # String vcf_output_path = read_string('vcf_output_path.txt')
        # String sites_only_vcf_output_path = read_string('sites_only_vcf_output_path.txt')
        File hail_gvs_import_script = 'hail_gvs_import.py'
        # File hail_export_tieout_vcf_script = 'hail_export_tieout_vcf.py'
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:rc_616_var_store_2022_09_06"
    }
}
