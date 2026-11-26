version 1.0

# Per-locus no-call sets built from the VDS, for computing AN over an arbitrary cohort (VS-2054).
# hail_create_nocall_sets.py, on a Dataproc cluster, writes one no-call set per locus of the densified VDS, with the
# VAT's filters applied, and the VDS keys that are not left-aligned with their carriers. This workflow then gives
# those keys their VAT VIDs and builds, in the given dataset:
#
#   <table_prefix>_sets                  one row per VDS locus: its no-call set as sample indices or as a bitmap
#   <table_prefix>_samples               sample index, name and default ploidy on non-PAR chrX
#   <table_prefix>_vid_overrides         for each VID that is not at its VDS locus, where to look it up and which
#                                        carriers to count as called there
#   <table_prefix>_vid_override_detail   how each override was derived
#   <table_prefix>_cohort_an             table function (sample_names, vids) -> (vid, cohort_size, an)
import "../wdl/GvsUtils.wdl" as Utils
import "GvsCreateParticipantMappingFromVDS.wdl" as ParticipantMapping

workflow GvsCreateNoCallSetsFromVDS {
    input {
        String vds_path
        String project_id
        String dataset_name
        String table_prefix
        String? contigs
        Int? num_test_cohorts
        String? nocall_sets_path
        Boolean use_existing_hail_outputs = false
        String reference_name = "hg38"

        String region = "us-central1"
        # Compute Engine stockouts are per-zone, and Dataproc's default placement picks one zone and fails rather than
        # moving on. "auto" tries every zone in the region.
        String cluster_zones = "auto"
        Int num_primary_workers = 2
        # Wider clusters have been unstable.
        Int max_secondary_workers = 200
        Int num_local_ssds = 1
        Int? num_master_local_ssds
        String worker_machine_type = "n2d-highmem-8"
        String master_machine_type = "n2d-highmem-32"
        Boolean use_tiny_dataproc_cluster = false
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Boolean leave_cluster_running_at_end = false
        Float? master_memory_fraction

        String? git_branch_or_tag
        String? hail_version
        File? hail_wheel
        String? basic_docker
        String? variants_docker
        String? cloud_sdk_slim_docker
        String? workspace_project
    }

    parameter_meta {
        vds_path: {
            help: "The VDS the VAT was built from. The no-call sets describe this VDS, so AN computed from them matches only a VAT built from the same one."
        }
        table_prefix: {
            help: "Prefix of every table this workflow builds, and of the cohort_an table function. Existing tables with the prefix are replaced. The staging tables (_sets_raw, _shifted_keys, _key_to_vid) expire after seven days."
        }
        contigs: {
            help: "Optional comma-separated contigs, e.g. chr21,chrX, to restrict the build to them. Default: the whole VDS."
        }
        num_test_cohorts: {
            help: "Optional. Also record hl.agg.call_stats AN for this many random cohorts at every locus (expected_an), to check AN from the sets against. For validation runs only: it widens every row."
        }
        nocall_sets_path: {
            help: "GCS prefix for the Hail job's outputs: nocall_sets.parquet, shifted_keys.parquet, keys_to_normalize.vcf.bgz, samples.tsv and summary.json. Defaults to a path in the workspace bucket named for the dataset and the table prefix."
        }
        use_existing_hail_outputs: {
            help: "Skip the Hail job and build the tables from outputs already at nocall_sets_path. The Hail job has no resume state of its own, so use this to rerun the table build after it has finished."
        }
        max_secondary_workers: {
            help: "Ceiling on the preemptible secondary workers the autoscaling policy may add."
        }
        hail_version: {
            help: "Optional Hail version. Cannot define both this parameter and `hail_wheel`."
        }
        hail_wheel: {
            help: "Optional Hail wheel file. Cannot define both this parameter and `hail_version`."
        }
        num_master_local_ssds: {
            help: "Local SSDs on the master, overriding num_local_ssds. Valid counts depend on machine family and vCPU count, and larger machines need more: N2D machines with 32 or 48 vCPUs reject 1 and accept 2, 4, 8, 16 or 24, but an N2 machine with 32 vCPUs needs at least 4. Left unset, the master gets the workers' count, except that 1 becomes 2, which suits the default n2d-highmem-32 master and n1-highmem-32. For any other master type, check the allowed counts in the Compute Engine local SSD documentation and set this explicitly."
        }
    }

    # Always called: the workspace bucket comes only from here.
    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_google_project = select_first([workspace_project, GetToolVersions.google_project])

    if (defined(hail_version) && defined(hail_wheel)) {
        call Utils.TerminateWorkflow as BothHailVersionAndHailWheelDefined {
            input:
                message = "Cannot define both `hail_version` and `hail_wheel`, exiting.",
                basic_docker = effective_basic_docker,
        }
    }

    if (use_existing_hail_outputs && !defined(nocall_sets_path)) {
        call Utils.TerminateWorkflow as ExistingOutputsNeedAPath {
            input:
                message = "`use_existing_hail_outputs` needs `nocall_sets_path`, exiting.",
                basic_docker = effective_basic_docker,
        }
    }

    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    String requested_nocall_sets_path = sub(select_first([nocall_sets_path,
        "~{GetToolVersions.workspace_bucket}/nocall_sets/~{dataset_name}/~{table_prefix}"]), "/$", "")

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = effective_basic_docker,
    }

    if (!use_existing_hail_outputs) {
        call Utils.GetHailScripts {
            input:
                variants_docker = effective_variants_docker,
        }

        call CreateNoCallSets {
            input:
                vds_path = vds_path,
                contigs = contigs,
                num_test_cohorts = num_test_cohorts,
                output_path = requested_nocall_sets_path,
                temp_path = "~{GetToolVersions.workspace_bucket}/hail-temp/nocall-sets",
                region = region,
                cluster_zones = cluster_zones,
                num_primary_workers = num_primary_workers,
                max_secondary_workers = max_secondary_workers,
                num_local_ssds = num_local_ssds,
                num_master_local_ssds = num_master_local_ssds,
                worker_machine_type = worker_machine_type,
                master_machine_type = master_machine_type,
                use_tiny_dataproc_cluster = use_tiny_dataproc_cluster,
                cluster_max_idle_minutes = cluster_max_idle_minutes,
                cluster_max_age_minutes = cluster_max_age_minutes,
                leave_cluster_running_at_end = leave_cluster_running_at_end,
                master_memory_fraction = master_memory_fraction,
                hail_version = effective_hail_version,
                hail_wheel = hail_wheel,
                run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
                hail_create_nocall_sets_script = GetHailScripts.hail_create_nocall_sets_script,
                hail_create_participant_mapping_inputs_script = GetHailScripts.hail_create_participant_mapping_inputs_script,
                hail_create_vat_inputs_script = GetHailScripts.hail_create_vat_inputs_script,
                create_vat_inputs_script = GetHailScripts.create_vat_inputs_script,
                workspace_project = effective_google_project,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
        }
    }

    # Through the task's output when the job ran, so that everything below waits for it.
    String hail_outputs = select_first([CreateNoCallSets.finished_output_path, requested_nocall_sets_path])

    call ParticipantMapping.NormalizeKeys {
        input:
            keys_to_normalize_vcf = "~{hail_outputs}/keys_to_normalize.vcf.bgz",
            ref = GetReference.reference.reference_fasta,
            variants_docker = effective_variants_docker,
    }

    call LoadNoCallSets {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            table_prefix = table_prefix,
            hail_outputs = hail_outputs,
            key_to_vid = NormalizeKeys.key_to_vid,
            variants_docker = effective_variants_docker,
    }

    call CreateNoCallTables {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            table_prefix = table_prefix,
            sets_raw_table = LoadNoCallSets.sets_raw_table,
            shifted_keys_table = LoadNoCallSets.shifted_keys_table,
            key_to_vid_table = LoadNoCallSets.key_to_vid_table,
            samples_table = LoadNoCallSets.samples_table,
            variants_docker = effective_variants_docker,
    }

    output {
        File summary = "~{hail_outputs}/summary.json"
        File key_to_vid = NormalizeKeys.key_to_vid
        File override_summary = CreateNoCallTables.override_summary
        String sets_table = CreateNoCallTables.sets_table
        String samples_table = CreateNoCallTables.fq_samples_table
        String vid_overrides_table = CreateNoCallTables.vid_overrides_table
        String vid_override_detail_table = CreateNoCallTables.vid_override_detail_table
        String cohort_an_function = CreateNoCallTables.cohort_an_function
    }
}

task CreateNoCallSets {
    input {
        String vds_path
        String? contigs
        Int? num_test_cohorts
        String output_path
        String temp_path

        String region
        String cluster_zones
        Int num_primary_workers
        Int max_secondary_workers
        Int num_local_ssds
        Int? num_master_local_ssds
        String worker_machine_type
        String master_machine_type
        Boolean use_tiny_dataproc_cluster
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Boolean leave_cluster_running_at_end
        Float? master_memory_fraction

        String hail_version
        File? hail_wheel
        File run_in_hail_cluster_script
        File hail_create_nocall_sets_script
        File hail_create_participant_mapping_inputs_script
        File hail_create_vat_inputs_script
        File create_vat_inputs_script
        String workspace_project
        String cloud_sdk_slim_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        apt-get update
        apt-get install --assume-yes python3.11-venv
        python3 -m venv ./localvenv
        . ./localvenv/bin/activate

        pip3 install --upgrade pip

        if [[ ! -z "~{hail_wheel}" ]]
        then
            pip3 install ~{hail_wheel}
        else
            pip3 install hail~{'==' + hail_version}
        fi

        pip3 install --upgrade google-cloud-dataproc ijson

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | od -h -An | tr -d '[:space:]')-$(head -c2 < /dev/urandom | od -h -An | tr -d '[:space:]')"
        cluster_name="nocall-sets-${hex}"

        # run_in_hail_cluster.py renders each key as `--key value`.
        cat > script-arguments.json <<FIN
        {
            "vds_input_path": "~{vds_path}",
            "output_path": "~{output_path}",
            "temp_path": "~{temp_path}"~{', "contigs": "' + contigs + '"'}~{', "num_test_cohorts": "' + num_test_cohorts + '"'}
        }
        FIN
        cat script-arguments.json

        # The script imports the mapping script for its key classes, and the VAT scripts for the VAT's filters.
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{hail_create_nocall_sets_script} \
            --secondary-script-path-list ~{hail_create_participant_mapping_inputs_script} \
            --secondary-script-path-list ~{hail_create_vat_inputs_script} \
            --secondary-script-path-list ~{create_vat_inputs_script} \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            ~{true='--use-tiny-dataproc-cluster' false='' use_tiny_dataproc_cluster} \
            --num-primary-workers ~{num_primary_workers} \
            --max-secondary-workers ~{max_secondary_workers} \
            --region ~{region} \
            --zones '~{cluster_zones}' \
            --num-local-ssds ~{num_local_ssds} \
            ~{'--num-master-local-ssds ' + num_master_local_ssds} \
            --worker-machine-type ~{worker_machine_type} \
            --master-machine-type ~{master_machine_type} \
            --workspace-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end}
    >>>

    runtime {
        memory: "6.5 GB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: 0
        docker: cloud_sdk_slim_docker
        bootDiskSizeGb: 10
    }

    output {
        String finished_output_path = output_path
    }
}

task LoadNoCallSets {
    input {
        String project_id
        String dataset_name
        String table_prefix
        String hail_outputs
        File key_to_vid
        String variants_docker
    }

    String sets_raw_table_name = "~{table_prefix}_sets_raw"
    String shifted_keys_table_name = "~{table_prefix}_shifted_keys"
    String key_to_vid_table_name = "~{table_prefix}_key_to_vid"
    String samples_table_name = "~{table_prefix}_samples"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        hail_outputs='~{hail_outputs}'
        bucket=$(echo "${hail_outputs}" | cut -d/ -f3)

        # Name the part files rather than globbing the Parquet directory: a BigQuery wildcard matches across /, and an
        # executor lost mid-write can leave partial files under _temporary/ after the job commits. One wildcard per
        # two-digit part-number prefix keeps the list to about a hundred URIs however many partitions there are.
        # `objects list` exits 0 with no output when nothing matches, and its own * does not cross /.
        part_wildcards() {
            gcloud storage objects list "${hail_outputs}/$1/part-*.parquet" --format='value(name)' |
                sed -E -e 's#(/part-[0-9]{2})[^/]*$#\1*.parquet#' -e "s#^#gs://${bucket}/#" | sort -u | paste -sd, -
        }

        # --parquet_enable_list_inference loads Spark's arrays as REPEATED columns rather than list.element structs.
        # Genome-wide on Foxtrot the sets are estimated at 4.5 TiB, under the 15 TB a single load job may read.
        sets=$(part_wildcards nocall_sets.parquet)
        if [[ -z "${sets}" ]]
        then
            echo "No part files under ${hail_outputs}/nocall_sets.parquet" >&2
            exit 1
        fi
        bq --apilog=false load --project_id=~{project_id} --replace --source_format=PARQUET \
            --parquet_enable_list_inference --clustering_fields=location \
            ~{dataset_name}.~{sets_raw_table_name} "${sets}"

        cat > shifted_keys_schema.json <<'FIN'
        [{"name": "location", "type": "INTEGER"},
         {"name": "vds_key", "type": "STRING"},
         {"name": "carrier_idxs", "type": "INTEGER", "mode": "REPEATED"}]
        FIN

        # A VDS with no key to normalize may leave nothing to load.
        shifted_keys=$(part_wildcards shifted_keys.parquet)
        if [[ -n "${shifted_keys}" ]]
        then
            bq --apilog=false load --project_id=~{project_id} --replace --source_format=PARQUET \
                --parquet_enable_list_inference ~{dataset_name}.~{shifted_keys_table_name} "${shifted_keys}"
        else
            bq --apilog=false rm -f --project_id=~{project_id} ~{dataset_name}.~{shifted_keys_table_name}
            bq --apilog=false mk --project_id=~{project_id} --table \
                ~{dataset_name}.~{shifted_keys_table_name} shifted_keys_schema.json
        fi

        if [[ -s ~{key_to_vid} ]]
        then
            bq --apilog=false load --project_id=~{project_id} --replace --source_format=CSV --field_delimiter='\t' \
                ~{dataset_name}.~{key_to_vid_table_name} ~{key_to_vid} vds_key:STRING,vid:STRING
        else
            bq --apilog=false rm -f --project_id=~{project_id} ~{dataset_name}.~{key_to_vid_table_name}
            bq --apilog=false mk --project_id=~{project_id} --table \
                ~{dataset_name}.~{key_to_vid_table_name} vds_key:STRING,vid:STRING
        fi

        # An explicit schema: AoU sample names are integers, and the table function takes names as strings.
        bq --apilog=false load --project_id=~{project_id} --replace --source_format=CSV --field_delimiter='\t' \
            --skip_leading_rows=1 ~{dataset_name}.~{samples_table_name} "${hail_outputs}/samples.tsv" \
            sample_idx:INTEGER,sample_name:STRING,chrx_nonpar_haploid_calls:INTEGER,chrx_nonpar_diploid_calls:INTEGER,chrx_default_ploidy:INTEGER

        for table in ~{sets_raw_table_name} ~{shifted_keys_table_name} ~{key_to_vid_table_name}; do
            bq --apilog=false update --project_id=~{project_id} --expiration=604800 ~{dataset_name}."${table}"
        done
    >>>

    meta {
        # The staging tables expire, so a relaunch must reload them rather than reuse a cached result naming tables
        # that may be gone.
        volatile: true
    }

    runtime {
        docker: variants_docker
    }

    output {
        String sets_raw_table = sets_raw_table_name
        String shifted_keys_table = shifted_keys_table_name
        String key_to_vid_table = key_to_vid_table_name
        String samples_table = samples_table_name
    }
}

task CreateNoCallTables {
    input {
        String project_id
        String dataset_name
        String table_prefix
        String sets_raw_table
        String shifted_keys_table
        String key_to_vid_table
        String samples_table
        String variants_docker
    }

    String fq_dataset = "~{project_id}.~{dataset_name}"
    String sets_table_name = "~{table_prefix}_sets"
    String detail_table_name = "~{table_prefix}_vid_override_detail"
    String overrides_table_name = "~{table_prefix}_vid_overrides"
    String function_name = "~{table_prefix}_cohort_an"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cat > create_nocall_tables.sql <<'SQL'
        CREATE OR REPLACE TABLE `~{fq_dataset}.~{sets_table_name}`
        CLUSTER BY location AS
        -- A set is sparse (nocall_sample_idxs) unless it has a bitmap. Sample i is bit (128 >> MOD(i, 8)) of byte
        -- DIV(i, 8), and the bitmap is padded to whole 32-bit words.
        SELECT * EXCEPT (nocall_bitmap_hex),
               nocall_bitmap_hex IS NULL AS is_sparse,
               FROM_HEX(nocall_bitmap_hex) AS nocall_bitmap
        FROM `~{fq_dataset}.~{sets_raw_table}`;

        -- One row per VID that at least one key not at its left-aligned position normalizes to.
        CREATE OR REPLACE TABLE `~{fq_dataset}.~{detail_table_name}` AS
        WITH keys AS (
            SELECT k.vid, s.vds_key, s.location AS key_location, s.carrier_idxs,
                   (CASE SPLIT(k.vid, '-')[OFFSET(0)]
                        WHEN 'X' THEN 23
                        WHEN 'Y' THEN 24
                        ELSE CAST(SPLIT(k.vid, '-')[OFFSET(0)] AS INT64) END) * 1000000000000 +
                   CAST(SPLIT(k.vid, '-')[OFFSET(1)] AS INT64) AS vid_location,
                   CONCAT(SPLIT(k.vid, '-')[OFFSET(2)], '-', SPLIT(k.vid, '-')[OFFSET(3)]) AS vid_alleles
            FROM `~{fq_dataset}.~{key_to_vid_table}` AS k
            JOIN `~{fq_dataset}.~{shifted_keys_table}` AS s USING (vds_key)
        ),
        carriers AS (
            SELECT vid, ARRAY_AGG(DISTINCT c ORDER BY c) AS carrier_idxs
            FROM keys, UNNEST(carrier_idxs) AS c
            GROUP BY vid
        ),
        vids AS (
            SELECT vid, vid_location,
                   ANY_VALUE(vid_alleles) AS vid_alleles,
                   COUNT(*) AS n_shifted_keys,
                   ARRAY_AGG(key_location ORDER BY key_location) AS key_locations,
                   IFNULL(ANY_VALUE(c.carrier_idxs), []) AS carrier_idxs
            -- LEFT: a key has no carriers when every genotype carrying it has GQ 0 or fails FT, and its VID still
            -- needs a lookup location.
            FROM keys LEFT JOIN carriers AS c USING (vid)
            GROUP BY vid, vid_location
        )
        SELECT v.vid, v.vid_location, v.n_shifted_keys, v.key_locations,
               ARRAY_LENGTH(v.carrier_idxs) AS n_carriers,
               ns.location IS NOT NULL AS row_at_vid_location,
               -- Whether the VID is also a key of its own at its location, or the row there is another variant's.
               IFNULL(v.vid_alleles IN UNNEST(ns.alt_alleles), FALSE) AS own_key_at_vid_location,
               ns.n_nocall AS n_nocall_at_vid_location,
               -- A sample that carries the VID at any locus is called, and any other sample is called if it is called
               -- at the VID's left-aligned location. With no row there, a single key's own locus stands in, which is
               -- what the VAT reports. With several keys and no row nothing is defined: lookup_location is NULL.
               CASE WHEN ns.location IS NOT NULL THEN v.vid_location
                    WHEN v.n_shifted_keys = 1 THEN v.key_locations[OFFSET(0)] END AS lookup_location,
               -- The carriers who are no-calls at the lookup location, to add back.
               IF(ns.location IS NULL, [], ARRAY(
                   SELECT c FROM UNNEST(v.carrier_idxs) AS c
                   WHERE IF(ns.is_sparse,
                            c IN UNNEST(ns.nocall_sample_idxs),
                            (TO_CODE_POINTS(SUBSTR(ns.nocall_bitmap, DIV(c, 8) + 1, 1))[OFFSET(0)]
                             & (128 >> MOD(c, 8))) != 0)
                   ORDER BY c)) AS rescued_idxs
        FROM vids AS v
        LEFT JOIN `~{fq_dataset}.~{sets_table_name}` AS ns ON ns.location = v.vid_location;

        CREATE OR REPLACE TABLE `~{fq_dataset}.~{overrides_table_name}`
        CLUSTER BY vid AS
        SELECT vid, lookup_location, rescued_idxs FROM `~{fq_dataset}.~{detail_table_name}`;
        SQL

        # bq query --max_rows check: ok, results go to the new tables
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} "$(cat create_nocall_tables.sql)"

        cat > create_cohort_an.sql <<'SQL'
        CREATE OR REPLACE TABLE FUNCTION `~{fq_dataset}.~{function_name}`(
            sample_names ARRAY<STRING>, vids ARRAY<STRING>) AS (
            WITH width AS (
                -- Stored bitmaps are padded to whole 32-bit words, and & needs operands of equal length.
                SELECT 4 * DIV(COUNT(*) + 31, 32) AS bitmap_bytes FROM `~{fq_dataset}.~{samples_table}`
            ),
            ploidy_groups AS (
                -- The cohort split by default ploidy on non-PAR chrX. Everywhere else every group is diploid.
                SELECT chrx_default_ploidy AS chrx_ploidy, COUNT(*) AS n, ARRAY_AGG(sample_idx) AS idxs
                FROM `~{fq_dataset}.~{samples_table}`
                WHERE sample_name IN UNNEST(sample_names)
                GROUP BY chrx_default_ploidy
            ),
            cohort AS (
                SELECT g.chrx_ploidy, g.n, g.idxs,
                       CODE_POINTS_TO_BYTES(ARRAY(
                           SELECT IFNULL(SUM(128 >> MOD(i, 8)), 0)
                           FROM UNNEST(GENERATE_ARRAY(0, w.bitmap_bytes - 1)) AS b
                           LEFT JOIN UNNEST(g.idxs) AS i ON DIV(i, 8) = b
                           GROUP BY b ORDER BY b)) AS bitmap
                FROM ploidy_groups AS g CROSS JOIN width AS w
            ),
            requested AS (
                SELECT vid,
                       -- hg38 PAR1 and PAR2.
                       contig_index = 23 AND NOT (position BETWEEN 10001 AND 2781479
                                                  OR position BETWEEN 155701383 AND 156030895) AS x_nonpar,
                       contig_index * 1000000000000 + position AS own_location
                FROM (
                    SELECT vid,
                           CAST(SPLIT(vid, '-')[OFFSET(1)] AS INT64) AS position,
                           CASE SPLIT(vid, '-')[OFFSET(0)]
                               WHEN 'X' THEN 23
                               WHEN 'Y' THEN 24
                               ELSE CAST(SPLIT(vid, '-')[OFFSET(0)] AS INT64) END AS contig_index
                    FROM UNNEST(vids) AS vid)
            )
            SELECT r.vid,
                   SUM(c.n) AS cohort_size,
                   -- NULL when there is no no-call set to look the VID up in.
                   SUM(IF(ns.location IS NULL, NULL, IF(r.x_nonpar, c.chrx_ploidy, 2) * (c.n
                       - IF(ns.is_sparse,
                            (SELECT COUNT(*) FROM UNNEST(ns.nocall_sample_idxs) AS id WHERE id IN UNNEST(c.idxs)),
                            BIT_COUNT(ns.nocall_bitmap & c.bitmap))
                       -- Carriers of the VID at another locus who are no-calls at this one count as called.
                       + (SELECT COUNT(*) FROM UNNEST(o.rescued_idxs) AS id WHERE id IN UNNEST(c.idxs))))) AS an
            FROM requested AS r
            CROSS JOIN cohort AS c
            LEFT JOIN `~{fq_dataset}.~{overrides_table_name}` AS o ON o.vid = r.vid
            LEFT JOIN `~{fq_dataset}.~{sets_table_name}` AS ns
                -- An override with a NULL lookup_location has no set; do not fall back to the VID's own location.
                ON ns.location = IF(o.vid IS NULL, r.own_location, o.lookup_location)
            GROUP BY r.vid
        );
        SQL

        # bq query --max_rows check: ok, creates a function
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} "$(cat create_cohort_an.sql)"

        # How the overridden VIDs divide among the cases. no_lookup_location counts those the rule does not cover.
        # bq query --max_rows check: ok, a single row of counts
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} --format=json \
            'SELECT COUNT(*) AS overridden_vids,
                    COUNTIF(n_shifted_keys > 1) AS several_keys,
                    COUNTIF(n_carriers = 0) AS no_carriers,
                    COUNTIF(row_at_vid_location AND own_key_at_vid_location) AS own_key_at_vid_location,
                    COUNTIF(row_at_vid_location AND NOT own_key_at_vid_location) AS another_variants_row_at_vid_location,
                    COUNTIF(NOT row_at_vid_location AND n_shifted_keys = 1) AS single_key_no_row,
                    COUNTIF(lookup_location IS NULL) AS no_lookup_location,
                    COUNTIF(lookup_location IS NULL AND n_carriers > 0) AS no_lookup_location_with_carriers,
                    COUNTIF(ARRAY_LENGTH(rescued_idxs) > 0) AS vids_with_rescued_carriers,
                    IFNULL(SUM(ARRAY_LENGTH(rescued_idxs)), 0) AS rescued_carriers
             FROM `~{fq_dataset}.~{detail_table_name}`' | jq '.[0] | map_values(tonumber)' > override_summary.json
        cat override_summary.json
    >>>

    runtime {
        docker: variants_docker
    }

    output {
        File override_summary = "override_summary.json"
        String sets_table = "~{fq_dataset}.~{sets_table_name}"
        String fq_samples_table = "~{fq_dataset}.~{samples_table}"
        String vid_overrides_table = "~{fq_dataset}.~{overrides_table_name}"
        String vid_override_detail_table = "~{fq_dataset}.~{detail_table_name}"
        String cohort_an_function = "~{fq_dataset}.~{function_name}"
    }
}
