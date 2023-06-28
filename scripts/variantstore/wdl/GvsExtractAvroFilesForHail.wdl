version 1.0

import "GvsUtils.wdl" as Utils


workflow GvsExtractAvroFilesForHail {
    input {
        Boolean go = true
        String project_id
        String dataset_name
        String filter_set_name
        String call_set_identifier
        Boolean use_VQSR_lite = true
        Int scatter_width = 10
    }

    String fq_gvs_dataset = "~{project_id}.~{dataset_name}"
    String filter_set_info_tablename = "filter_set_info"
    String fq_filter_set_info_table = "~{fq_gvs_dataset}.~{filter_set_info_tablename}"

    call Utils.ValidateFilterSetName {
        input:
            project_id = project_id,
            fq_filter_set_info_table = "~{fq_filter_set_info_table}",
            filter_set_name = filter_set_name
    }

    call Utils.IsVQSRLite {
        input:
            project_id = project_id,
            fq_filter_set_info_table = "~{project_id}.~{dataset_name}.filter_set_info",
            filter_set_name = filter_set_name
    }

    call OutputPath {
        input: go = ValidateFilterSetName.done
    }

    call ExtractFromNonSuperpartitionedTables {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_info_tablename = filter_set_info_tablename,
            filter_set_name = filter_set_name,
            avro_sibling = OutputPath.out,
            call_set_identifier = call_set_identifier,
            is_vqsr_lite = IsVQSRLite.is_vqsr_lite
    }

    call Utils.CountSuperpartitions {
        input:
            project_id = project_id,
            dataset_name = dataset_name
    }

    scatter (i in range(scatter_width)) {
        call ExtractFromSuperpartitionedTables {
            input:
                project_id = project_id,
                dataset_name = dataset_name,
                call_set_identifier = call_set_identifier,
                avro_sibling = OutputPath.out,
                num_superpartitions = CountSuperpartitions.num_superpartitions,
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
        File hail_create_vat_inputs_script = GenerateHailScripts.hail_create_vat_inputs_script
        String vds_output_path = GenerateHailScripts.vds_output_path
        String sites_only_vcf_output_path = GenerateHailScripts.sites_only_vcf_output_path
        String vat_inputs_output_path = GenerateHailScripts.vat_inputs_output_path
        String avro_prefix = ExtractFromNonSuperpartitionedTables.output_prefix
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
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
        disks: "local-disk 500 HDD"
    }
}


task ExtractFromNonSuperpartitionedTables {
    meta {
        description: "Extracts from the non-superpartitioned tables: sample_info, filter_set_sites, filter_set_info/filter_set_info_vqsr_lite, and filter_set_tranches (if using VQSR Classic)"
        # Not dealing with caching for now as that would introduce a lot of complexity.
        volatile: true
    }
    input {
        String project_id
        String dataset_name
        String filter_set_info_tablename
        String filter_set_name
        String avro_sibling
        String call_set_identifier
        Boolean is_vqsr_lite = true
    }

    String vqs_score_field = if (is_vqsr_lite == true) then 'calibration_sensitivity' else 'vqslod'

    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
    }
    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        echo $avro_prefix > "avro_prefix.out"

        python3 /app/run_avro_query.py --sql "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/sample_mapping/sample_mapping_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT sample_id, sample_name, '40',
            'gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list' as intervals_file
            FROM \`~{project_id}.~{dataset_name}.sample_info\`
            WHERE withdrawn IS NULL AND
            is_control = false
            ORDER BY sample_id
        " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name sample_info --project_id ~{project_id}

        python3 /app/run_avro_query.py --sql "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vqsr_filtering_data/vqsr_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, type as model, ref, alt, ~{vqs_score_field}, yng_status
            FROM \`~{project_id}.~{dataset_name}.filter_set_info\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name ~{filter_set_info_tablename} --project_id ~{project_id}

        python3 /app/run_avro_query.py --sql "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/site_filtering_data/site_filtering_data_*.avro', format='AVRO', compression='SNAPPY') AS
            SELECT location, filters
            FROM \`~{project_id}.~{dataset_name}.filter_set_sites\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name filter_set_sites --project_id ~{project_id}

        if [ ~{is_vqsr_lite} = false ]; then
            python3 /app/run_avro_query.py --sql "
                EXPORT DATA OPTIONS(
                uri='${avro_prefix}/vqsr_tranche_data/vqsr_tranche_data_*.avro', format='AVRO', compression='SNAPPY') AS
                SELECT model, truth_sensitivity, min_vqslod, filter_name
                FROM \`~{project_id}.~{dataset_name}.filter_set_tranches\`
                WHERE filter_set_name = '~{filter_set_name}'
            " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name filter_set_tranches --project_id ~{project_id}
        fi
    >>>

    output {
        Boolean done = true
        String output_prefix = read_string("avro_prefix.out")
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        disks: "local-disk 500 HDD"
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
        String dataset_name
        String avro_sibling
        String call_set_identifier
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
        avro_prefix="$(dirname ~{avro_sibling})/avro"

        for superpartition in $(seq ~{shard_index + 1} ~{num_shards} ~{num_superpartitions})
        do
            str_table_index=$(printf "%03d" $superpartition)

            # These bq exports error out if there are any objects at the sibling level to where output files would be written
            # so an extra layer of `vet_${str_table_index}` is inserted here.
            python3 /app/run_avro_query.py --sql "
                EXPORT DATA OPTIONS(
                uri='${avro_prefix}/vets/vet_${str_table_index}/vet_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
                SELECT location, v.sample_id, ref, REPLACE(alt,',<NON_REF>','') alt, call_GT as GT, call_AD as AD, call_GQ as GQ, cast(SPLIT(call_pl,',')[OFFSET(0)] as int64) as RGQ
                FROM \`~{project_id}.~{dataset_name}.vet_${str_table_index}\` v
                INNER JOIN \`~{project_id}.~{dataset_name}.sample_info\` s ON s.sample_id = v.sample_id
                WHERE withdrawn IS NULL AND
                is_control = false
                ORDER BY location
            " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name vet_${str_table_index} --project_id ~{project_id}

            python3 /app/run_avro_query.py --sql "
                EXPORT DATA OPTIONS(
                uri='${avro_prefix}/refs/ref_ranges_${str_table_index}/ref_ranges_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY') AS
                SELECT location, r.sample_id, length, state
                FROM \`~{project_id}.~{dataset_name}.ref_ranges_${str_table_index}\` r
                INNER JOIN \`~{project_id}.~{dataset_name}.sample_info\` s ON s.sample_id = r.sample_id
                WHERE withdrawn IS NULL AND
                is_control = false
                ORDER BY location
            " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name ref_ranges_${str_table_index} --project_id ~{project_id}
        done
    >>>

    output {
        Boolean done = true
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        disks: "local-disk 500 HDD"
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

        # 4 random hex bytes to not clobber outputs if this is run multiple times for the same avro_prefix.
        # Unlike many implementations, at the time of this writing this works on both Debian and Alpine based images
        # so it should continue to work even if the `variantstore` image switches to Alpine:
        # https://stackoverflow.com/a/34329799
        rand=$(od -vN 4 -An -tx1 /dev/urandom | tr -d " ")

        # The write prefix will be a sibling to the Avro "directory" that embeds the current date and some randomness.
        write_prefix="$(dirname ~{avro_prefix})/$(date -Idate)-${rand}"

        vds_output_path="${write_prefix}/gvs_export.vds"
        echo $vds_output_path > vds_output_path.txt

        tmpfile=$(mktemp)
        # `sed` can use delimiters other than `/`. This is required here since the replacement GCS paths will
        # contain `/` characters.
        cat /app/hail_gvs_import.py |
            sed "s;@AVRO_PREFIX@;~{avro_prefix};" |
            sed "s;@WRITE_PREFIX@;${write_prefix};" > ${tmpfile}
        mv ${tmpfile} hail_gvs_import.py

        vcf_output_path="${write_prefix}/gvs_export.vcf"
        echo $vcf_output_path > vcf_output_path.txt
        sites_only_vcf_output_path="${write_prefix}/gvs_sites_only.vcf"
        echo $sites_only_vcf_output_path > sites_only_vcf_output_path.txt
        vat_tsv_output_path="${write_prefix}/vat_inputs.tsv"
        echo $vat_tsv_output_path > vat_inputs_output_path.txt

        tmpfile=$(mktemp)
        cat /app/hail_create_vat_inputs.py |
            sed "s;@VDS_INPUT_PATH@;${vds_output_path};" |
            sed "s;@SITES_ONLY_VCF_OUTPUT_PATH@;${sites_only_vcf_output_path};" |
            sed "s;@VAT_CUSTOM_ANNOTATIONS_OUTPUT_PATH@;${vat_tsv_output_path};" > ${tmpfile}
        mv ${tmpfile} hail_create_vat_inputs.py

    >>>

    output {
        Boolean done = true
        String vds_output_path = read_string('vds_output_path.txt')
        String sites_only_vcf_output_path = read_string('sites_only_vcf_output_path.txt')
        String vat_inputs_output_path = read_string('vat_inputs_output_path.txt')
        File hail_gvs_import_script = 'hail_gvs_import.py'
        File hail_create_vat_inputs_script = 'hail_create_vat_inputs.py'
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        disks: "local-disk 500 HDD"
    }
}
