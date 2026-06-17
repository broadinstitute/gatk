version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsExtractAvroFilesForHail {
    input {

        # A sample table or view is expected to exist in the dataset that has at least the following columns:
        # sample_id (int), sample_name (string), is_control (boolean), and withdrawn (type unspecified but nullable).
        # `sample_info` fits this bill, but so would a view that filters `sample_info` or joins to other tables, so long
        # as the result contains those four columns and does *not* contain duplicate rows. e.g. the
        # following creates a view that contains only samples that have an allele in the alt_allele table at the same
        # location corresponding to an "accursed" sample with erroneous reference data at that location:
        #
        #  CREATE OR REPLACE VIEW `project.dataset.sample_info_vs_1862_accursed_sample`
        #  AS
        #  SELECT si.*
        #  FROM `project.dataset.sample_info` si WHERE si.sample_id IN (
        #      -- Note the DISTINCT is important here as otherwise split hetvars in alt_allele will cause duplicate rows
        #      -- in the view which would then cause downstream errors.
        #      SELECT distinct aa.sample_id from `project.dataset.alt_allele` aa
        #      WHERE aa.location = 4 * 1000 * 1000 * 1000 * 1000 + 190181397
        #  )

        # `sample_table_or_view_name` must be the *unqualified* name of the sample table or view, i.e. without the `project.dataset` prefix.
        String sample_table_or_view_name = "sample_info"
        String? git_branch_or_tag
        String? git_hash
        # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go = true
        String project_id
        String dataset_name
        String filter_set_name
        String call_set_identifier
        Int scatter_width = 10
        String? basic_docker
        String? cloud_sdk_docker
        String? variants_docker
        String ploidy_table_name = "sample_chromosome_ploidy"
    }

    if (!defined(git_hash) || !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(variants_docker)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

    String fq_gvs_dataset = "~{project_id}.~{dataset_name}"
    String filter_set_info_tablename = "filter_set_info"
    String fq_filter_set_info_table = "~{fq_gvs_dataset}.~{filter_set_info_tablename}"

    call Utils.ValidateFilterSetName {
        input:
            project_id = project_id,
            fq_filter_set_info_table = "~{fq_filter_set_info_table}",
            filter_set_name = filter_set_name,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call OutputPath {
        input:
            go = ValidateFilterSetName.done,
            basic_docker = effective_basic_docker,
    }

    call ExtractFromSampleTable {
        input:
            sample_table_or_view_name = sample_table_or_view_name,
            project_id = project_id,
            dataset_name = dataset_name,
            avro_sibling = OutputPath.out,
            call_set_identifier = call_set_identifier,
            variants_docker = effective_variants_docker,
    }

    call ExtractFromFilterTables {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            filter_set_info_tablename = filter_set_info_tablename,
            filter_set_name = filter_set_name,
            avro_sibling = OutputPath.out,
            call_set_identifier = call_set_identifier,
            variants_docker = effective_variants_docker,
    }

    call ExtractFromPloidyTable {
        input:
            sample_table_or_view_name = sample_table_or_view_name,
            project_id = project_id,
            dataset_name = dataset_name,
            ploidy_table_name = ploidy_table_name,
            avro_sibling = OutputPath.out,
            call_set_identifier = call_set_identifier,
            variants_docker = effective_variants_docker,
    }

    call Utils.CountSuperpartitions {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    call Utils.GetBQTableLastModifiedDatetime as RefTableDatetimeCheck {
        input:
            project_id = project_id,
            fq_table = "~{project_id}.~{dataset_name}.ref_ranges_001",
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }


    call Utils.IsUsingCompressedReferences {
        input:
            query_project_id = project_id,
            dest_project_id = project_id,
            dataset_name = dataset_name,
            ref_table_timestamp = RefTableDatetimeCheck.last_modified_timestamp,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    scatter (i in range(scatter_width)) {
        call ExtractFromSuperpartitionedTables {
            input:
                sample_table_or_view_name = sample_table_or_view_name,
                project_id = project_id,
                dataset_name = dataset_name,
                call_set_identifier = call_set_identifier,
                avro_sibling = OutputPath.out,
                num_superpartitions = CountSuperpartitions.num_superpartitions,
                shard_index = i,
                num_shards = scatter_width,
                variants_docker = effective_variants_docker,
                use_compressed_references = IsUsingCompressedReferences.is_using_compressed_references,
        }
    }

    output {
        String avro_path = ExtractFromFilterTables.output_prefix
        String recorded_git_hash = effective_git_hash
    }
}


task OutputPath {
    meta {
        description: "Does nothing but produce the cloud path to its stdout."
        # Always make a new output path, otherwise every invocation would clobber the original.
        volatile: true
    }
    input {
        # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go = true
        String basic_docker
    }
    command <<<
    >>>
    output {
        File out = stdout()
    }
    runtime {
        docker: basic_docker
        disks: "local-disk 500 HDD"
    }
}

# splitting out the extract sample table Avros into its own task to create partition-based files
# this might not be the most efficient for callsets under a certain size
task ExtractFromSampleTable {
    meta {
        description: "Extracts from the sample table, split up by partition"
        # Not dealing with caching for now as that would introduce a lot of complexity.
        volatile: true
    }
    input {
        String sample_table_or_view_name
        String project_id
        String dataset_name
        String avro_sibling
        String call_set_identifier
        String variants_docker
    }

    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        echo $avro_prefix > "avro_prefix.out"

        python3 /app/run_avro_query_for_sample_table.py \
            --avro_prefix ${avro_prefix} \
            --call_set_identifier ~{call_set_identifier} \
            --dataset_name ~{dataset_name} \
            --project_id=~{project_id} \
            --sample_table_name ~{sample_table_or_view_name}
    >>>

    output {
        Boolean done = true
        String output_prefix = read_string("avro_prefix.out")
    }

    runtime {
        docker: variants_docker
        disks: "local-disk 500 HDD"
    }
}


task ExtractFromFilterTables {
    meta {
        description: "Extracts from the tables: filter_set_sites and filter_set_info"
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
        String variants_docker
    }

    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        echo $avro_prefix > "avro_prefix.out"

        python3 /app/run_avro_query.py --sql "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/vets_filtering_data/vets_filtering_data_*.avro', format='AVRO', compression='SNAPPY', overwrite=true) AS
            SELECT location, type as model, ref, alt, calibration_sensitivity, yng_status
            FROM \`~{project_id}.~{dataset_name}.filter_set_info\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name ~{filter_set_info_tablename} --project_id=~{project_id}

        python3 /app/run_avro_query.py --sql "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/site_filtering_data/site_filtering_data_*.avro', format='AVRO', compression='SNAPPY', overwrite=true) AS
            SELECT location, filters
            FROM \`~{project_id}.~{dataset_name}.filter_set_sites\`
            WHERE filter_set_name = '~{filter_set_name}'
            ORDER BY location
        " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name filter_set_sites --project_id=~{project_id}
    >>>

    output {
        Boolean done = true
        String output_prefix = read_string("avro_prefix.out")
    }

    runtime {
        docker: variants_docker
        disks: "local-disk 500 HDD"
    }
}


task ExtractFromPloidyTable {
    meta {
        description: "Extracts from the sample chromosome ploidy table"
        # Not dealing with caching for now as that would introduce a lot of complexity.
        volatile: true
    }
    input {
        String sample_table_or_view_name
        String project_id
        String dataset_name
        String ploidy_table_name
        String avro_sibling
        String call_set_identifier
        String variants_docker
    }

    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        avro_prefix="$(dirname ~{avro_sibling})/avro"
        echo $avro_prefix > "avro_prefix.out"

        # Note the query below extracts ploidy data for chrX and chrY only as those are the only chromosomes the VDS
        # ploidy logic looks at.

        python3 /app/run_avro_query.py --sql "
            EXPORT DATA OPTIONS(
            uri='${avro_prefix}/ploidy_data/ploidy_data_*.avro', format='AVRO', compression='SNAPPY', overwrite=true) AS
            SELECT (
                CASE (p.chromosome / 1000000000000)
                    WHEN 23 THEN 'chrX'
                    WHEN 24 THEN 'chrY'
                    END) AS location, s.sample_name, p.ploidy
            FROM \`~{project_id}.~{dataset_name}.~{ploidy_table_name}\` p
            JOIN \`~{project_id}.~{dataset_name}.~{sample_table_or_view_name}\` s ON p.sample_id = s.sample_id
            WHERE (p.chromosome / 1000000000000 = 23 or p.chromosome / 1000000000000 = 24)
            AND s.withdrawn IS NULL
            AND s.is_control = false
        " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name ~{ploidy_table_name} --project_id=~{project_id}
    >>>
    output {
        Boolean done = true
        String output_prefix = read_string("avro_prefix.out")
    }

    runtime {
        docker: variants_docker
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
        String sample_table_or_view_name
        String project_id
        String dataset_name
        String avro_sibling
        String call_set_identifier
        Int num_superpartitions
        Int shard_index
        Int num_shards
        String variants_docker
        Boolean use_compressed_references = false
    }

    parameter_meta {
        avro_sibling: "Cloud path to a file that will be the sibling to the 'avro' 'directory' under which output Avro files will be written."
        num_superpartitions: "Total number of superpartitions requiring extraact"
        shard_index: "0-based index of this superpartition extract shard"
        num_shards: "Count of all superpartition extract shards"
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        avro_prefix="$(dirname ~{avro_sibling})/avro"

        # Query the sample table upfront to find only the superpartitions that (a) have at least one
        # non-withdrawn, non-control sample and (b) are assigned to this shard. A superpartition p
        # belongs to shard shard_index when (p - 1) % num_shards == shard_index, which is the same
        # assignment produced by `seq shard_index+1 num_shards num_superpartitions`. Skipping empty
        # superpartitions avoids running EXPORT DATA queries that would produce no output.
        python3 /app/run_avro_query.py \
            --sql '
                SELECT superpartition
                FROM (
                    SELECT DISTINCT CAST(CEIL(sample_id / 4000.0) AS INT64) AS superpartition
                    FROM `~{project_id}.~{dataset_name}.~{sample_table_or_view_name}`
                    WHERE withdrawn IS NULL AND is_control = false
                )
                WHERE MOD(superpartition - 1, ~{num_shards}) = ~{shard_index}
                  AND superpartition <= ~{num_superpartitions}
                ORDER BY superpartition
            ' \
            --call_set_identifier ~{call_set_identifier} \
            --dataset_name ~{dataset_name} \
            --table_name ~{sample_table_or_view_name} \
            --project_id ~{project_id} \
            --output_file superpartitions_to_process.txt

        readarray -t superpartitions_to_process < superpartitions_to_process.txt

        for superpartition in "${superpartitions_to_process[@]}"
        do
            str_table_index=$(printf "%03d" $superpartition)

            # These bq exports error out if there are any objects at the sibling level to where output files would be written
            # so an extra layer of `vet_${str_table_index}` is inserted here.
            python3 /app/run_avro_query.py --sql "
                EXPORT DATA OPTIONS(
                uri='${avro_prefix}/vets/vet_${str_table_index}/vet_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY', overwrite=true) AS
                SELECT location, v.sample_id, ref, REPLACE(alt,',<NON_REF>','') alt, call_GT as GT, call_AD as AD, call_GQ as GQ, cast(SPLIT(call_pl,',')[OFFSET(0)] as int64) as RGQ,
                call_PS as PS
                FROM \`~{project_id}.~{dataset_name}.vet_${str_table_index}\` v
                INNER JOIN \`~{project_id}.~{dataset_name}.~{sample_table_or_view_name}\` s ON s.sample_id = v.sample_id
                WHERE withdrawn IS NULL AND
                is_control = false
                ORDER BY location
            " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name vet_${str_table_index} --project_id=~{project_id}

            if [ ~{use_compressed_references} = false ]; then
                python3 /app/run_avro_query.py --sql "
                    EXPORT DATA OPTIONS(
                    uri='${avro_prefix}/refs/ref_ranges_${str_table_index}/ref_ranges_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY', overwrite=true) AS
                    SELECT location, r.sample_id, length, state
                    FROM \`~{project_id}.~{dataset_name}.ref_ranges_${str_table_index}\` r
                    INNER JOIN \`~{project_id}.~{dataset_name}.~{sample_table_or_view_name}\` s ON s.sample_id = r.sample_id
                    WHERE withdrawn IS NULL AND
                    is_control = false
                    ORDER BY location
                " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name ref_ranges_${str_table_index} --project_id ~{project_id}
            else
                python3 /app/run_avro_query.py --sql "
                    CREATE TEMP FUNCTION intToState(state INT64)
                        RETURNS string
                        AS (
                          CASE state
                            WHEN 7 THEN 'v'
                            WHEN 8 THEN '*'
                            WHEN 9 THEN 'm'
                            WHEN 10 THEN 'u'
                            ELSE CAST(state as string)
                          END
                    );
                    CREATE TEMP FUNCTION UnpackRefRangeInfo(superpackEntry int64)
                        RETURNS STRUCT<location INT64, len INT64, state string>
                        AS (
                          STRUCT(
                          1000000000000 * ((superpackEntry >> 48) & 0xFFFF) + ((superpackEntry >> 16) & 0xFFFFFFFF),
                          (superpackEntry >> 4) & 0xFFF,
                          intToState(superpackEntry & 0xF))
                    );
                    EXPORT DATA OPTIONS(
                    uri='${avro_prefix}/refs/ref_ranges_${str_table_index}/ref_ranges_${str_table_index}_*.avro', format='AVRO', compression='SNAPPY', overwrite=true) AS
                    SELECT UnpackRefRangeInfo(packed_ref_data).location as location, r.sample_id, UnpackRefRangeInfo(packed_ref_data).len as length, UnpackRefRangeInfo(packed_ref_data).state as state
                    FROM \`~{project_id}.~{dataset_name}.ref_ranges_${str_table_index}\` r
                    INNER JOIN \`~{project_id}.~{dataset_name}.~{sample_table_or_view_name}\` s ON s.sample_id = r.sample_id
                    WHERE withdrawn IS NULL AND
                    is_control = false
                    ORDER BY location
                    " --call_set_identifier ~{call_set_identifier} --dataset_name ~{dataset_name} --table_name ref_ranges_${str_table_index} --project_id ~{project_id}
            fi
        done
    >>>

    output {
        Boolean done = true
    }

    runtime {
        docker: variants_docker
        disks: "local-disk 500 HDD"
        noAddress: true
    }
}
