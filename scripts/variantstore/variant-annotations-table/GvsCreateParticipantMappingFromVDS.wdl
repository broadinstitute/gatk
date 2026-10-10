version 1.0

# Participant mapping built from the VDS (VS-2013). hail_create_participant_mapping_inputs.py, on a Dataproc cluster,
# writes the carrier pairs and the VDS keys that are not left-aligned; this workflow then gives those keys their VAT
# VIDs and builds, in the given dataset, `<participant_mapping_table_name>_base` (het_ids, hom_ids, hemi_ids per VID),
# its provenance table `<participant_mapping_table_name>_provenance`, and the view `<participant_mapping_table_name>`
# with the delivered (vid, person_ids) shape.
import "../wdl/GvsUtils.wdl" as Utils

workflow GvsCreateParticipantMappingFromVDS {
    input {
        String vds_path
        String project_id
        String dataset_name
        String fq_vat_table
        String participant_mapping_table_name
        String? interval
        String? mapping_inputs_path
        String reference_name = "hg38"

        String region = "us-central1"
        # Compute Engine stockouts are per-zone, and Dataproc's default placement picks one zone and fails rather than
        # moving on. "auto" tries every zone in the region.
        String cluster_zones = "auto"
        Int num_primary_workers = 2
        # Wider clusters have been unstable.
        Int max_secondary_workers = 200
        Int num_local_ssds = 1
        String worker_machine_type = "n1-highmem-8"
        String master_machine_type = "n1-highmem-32"
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
            help: "The VDS the VAT was built from."
        }
        fq_vat_table: {
            help: "The VAT, as project.dataset.table. Every VID the mapping files a carrier under must appear in it."
        }
        interval: {
            help: "Optional single locus interval, e.g. chr21, to restrict the build to it. The acceptance checks restrict the VAT to the contigs the mapping covers, so pass whole contigs."
        }
        mapping_inputs_path: {
            help: "GCS prefix for the Hail job's outputs: keys_to_normalize.vcf.bgz, pairs.parquet and summary.json. Defaults to a path in the workspace bucket named for the dataset and the view. A relaunch with the same VDS, interval and prefix resumes from the job's checkpoint."
        }
        participant_mapping_table_name: {
            help: "Name of the delivered view. The base and provenance tables, and the staging tables, take it as a prefix. None may already exist, apart from the staging tables, which are replaced and expire after seven days."
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

    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    # Deterministic, unlike the hex-suffixed paths elsewhere, so that a relaunch finds the Hail job's resume state.
    String mapping_inputs = sub(select_first([mapping_inputs_path,
        "~{GetToolVersions.workspace_bucket}/participant_mapping_inputs/~{dataset_name}/~{participant_mapping_table_name}"]), "/$", "")

    call Utils.GetHailScripts {
        input:
            variants_docker = effective_variants_docker,
    }

    call Utils.GetReference {
        input:
            reference_name = reference_name,
            basic_docker = effective_basic_docker,
    }

    call CreateMappingInputs {
        input:
            vds_path = vds_path,
            interval = interval,
            output_path = mapping_inputs,
            temp_path = "~{GetToolVersions.workspace_bucket}/hail-temp/participant-mapping-inputs",
            region = region,
            cluster_zones = cluster_zones,
            num_primary_workers = num_primary_workers,
            max_secondary_workers = max_secondary_workers,
            num_local_ssds = num_local_ssds,
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
            hail_create_participant_mapping_inputs_script = GetHailScripts.hail_create_participant_mapping_inputs_script,
            hail_create_vat_inputs_script = GetHailScripts.hail_create_vat_inputs_script,
            create_vat_inputs_script = GetHailScripts.create_vat_inputs_script,
            workspace_project = effective_google_project,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
    }

    call NormalizeKeys {
        input:
            keys_to_normalize_vcf = CreateMappingInputs.keys_to_normalize_vcf,
            ref = GetReference.reference.reference_fasta,
            variants_docker = effective_variants_docker,
    }

    call LoadMappingInputs {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            participant_mapping_table_name = participant_mapping_table_name,
            pairs_parquet_path = CreateMappingInputs.pairs_parquet_path,
            key_to_vid = NormalizeKeys.key_to_vid,
            variants_docker = effective_variants_docker,
    }

    call CreateMappingTables {
        input:
            project_id = project_id,
            dataset_name = dataset_name,
            participant_mapping_table_name = participant_mapping_table_name,
            pairs_table = LoadMappingInputs.pairs_table,
            normalized_pairs_table = LoadMappingInputs.normalized_pairs_table,
            key_to_vid_table = LoadMappingInputs.key_to_vid_table,
            variants_docker = effective_variants_docker,
    }

    call CheckMappingTables {
        input:
            project_id = project_id,
            fq_vat_table = fq_vat_table,
            base_table = CreateMappingTables.base_table,
            provenance_table = CreateMappingTables.provenance_table,
            variants_docker = effective_variants_docker,
    }

    output {
        File summary = CreateMappingInputs.summary
        File key_to_vid = NormalizeKeys.key_to_vid
        File acceptance_checks = CheckMappingTables.acceptance_checks
        String base_table = CreateMappingTables.base_table
        String provenance_table = CreateMappingTables.provenance_table
        String view = CreateMappingTables.view
    }
}

task CreateMappingInputs {
    input {
        String vds_path
        String? interval
        String output_path
        String temp_path

        String region
        String cluster_zones
        Int num_primary_workers
        Int max_secondary_workers
        Int num_local_ssds
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
        cluster_name="participant-mapping-${hex}"

        # run_in_hail_cluster.py renders each key as `--key value`.
        cat > script-arguments.json <<FIN
        {
            "vds": "~{vds_path}",
            "output": "~{output_path}",
            "temp-path": "~{temp_path}"~{', "interval": "' + interval + '"'}
        }
        FIN
        cat script-arguments.json

        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{hail_create_participant_mapping_inputs_script} \
            --secondary-script-path-list ~{hail_create_vat_inputs_script} \
            --secondary-script-path-list ~{create_vat_inputs_script} \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            ~{true='--use-tiny-dataproc-cluster' false='' use_tiny_dataproc_cluster} \
            --num-primary-workers ~{num_primary_workers} \
            --max-secondary-workers ~{max_secondary_workers} \
            --region ~{region} \
            --zones ~{cluster_zones} \
            --num-local-ssds ~{num_local_ssds} \
            --worker-machine-type ~{worker_machine_type} \
            --master-machine-type ~{master_machine_type} \
            --workspace-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end}

        echo "~{output_path}/keys_to_normalize.vcf.bgz" > keys_to_normalize_vcf.txt
        echo "~{output_path}/summary.json" > summary.txt
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
        File keys_to_normalize_vcf = read_string("keys_to_normalize_vcf.txt")
        File summary = read_string("summary.txt")
        String pairs_parquet_path = "~{output_path}/pairs.parquet"
    }
}

task NormalizeKeys {
    input {
        File keys_to_normalize_vcf
        File ref
        String variants_docker
    }

    Int disk_size = ceil(size(ref, "GB") * 2 + size(keys_to_normalize_vcf, "GB") * 5) + 20

    # Writes key_to_vid.tsv, one headerless `vds_key<TAB>vid` row per input key.
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # The same normalization RemoveDuplicatesFromSitesOnlyVCF in GvsCreateVATfromVDS.wdl applies, against the same
        # reference, so each key gets the VID the VAT gave its allele. ID passes through untouched.
        bcftools norm -m- --check-ref w -f ~{ref} ~{keys_to_normalize_vcf} -O v -o normalized.vcf

        # VID format: contig without the chr prefix.
        bcftools query -f '%ID\t%CHROM-%POS-%REF-%ALT\n' normalized.vcf \
            | awk -F'\t' -v OFS='\t' '{ sub(/^chr/, "", $2); print }' > key_to_vid.tsv

        # Every key is biallelic, so norm -m- neither splits nor drops; anything else means the input is not what the
        # Hail script writes.
        input_keys=$(bcftools view -H ~{keys_to_normalize_vcf} | wc -l)
        output_rows=$(wc -l < key_to_vid.tsv)
        distinct_keys=$(cut -f1 key_to_vid.tsv | sort -u | wc -l)
        if [[ ${input_keys} -ne ${output_rows} || ${input_keys} -ne ${distinct_keys} ]]
        then
            echo "Expected one row per key: ${input_keys} keys in, ${output_rows} rows and ${distinct_keys} distinct keys out." >&2
            exit 1
        fi

        # The Hail script flags a key only when its indel can move left, so every one must come out at a new VID.
        unmoved=$(awk -F'\t' '$1 == $2' key_to_vid.tsv | wc -l)
        if [[ ${unmoved} -ne 0 ]]
        then
            echo "${unmoved} flagged keys normalized to themselves, e.g.:" >&2
            awk -F'\t' '$1 == $2' key_to_vid.tsv | head -5 >&2
            exit 1
        fi
    >>>

    runtime {
        docker: variants_docker
        maxRetries: 3
        memory: "4 GB"
        preemptible: 3
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File key_to_vid = "key_to_vid.tsv"
    }
}

task LoadMappingInputs {
    input {
        String project_id
        String dataset_name
        String participant_mapping_table_name
        String pairs_parquet_path
        File key_to_vid
        String variants_docker
    }

    String pairs_table_name = "~{participant_mapping_table_name}_pairs"
    String normalized_pairs_table_name = "~{participant_mapping_table_name}_normalized_pairs"
    String key_to_vid_table_name = "~{participant_mapping_table_name}_key_to_vid"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Staging tables, replaced on a rerun and expiring after a week. Load jobs are free; the pairs are on the order
        # of 2.5e12 rows and 11 TB of Parquet genome-wide, close to the 15 TB a single load job may read, so they load
        # in batches. Each batch is a set of wildcards over the part files, one per two-digit part-number prefix in
        # either flag directory: at most 200 URIs, so no batch argument approaches the kernel's per-argument limit,
        # and Hail writes partitions in locus order, so each batch holds a run of nearby keys. Name the part files
        # rather than globbing pairs.parquet: a BigQuery wildcard matches across /, and an executor lost mid-write can
        # leave partial files under _temporary/ after the job commits. Spark writes no directory for a flag value with
        # no rows, which an interval with no flagged carriers can produce; listing the part files handles that.
        # Without hive partitioning the flag, which lives only in the directory names, is not loaded.
        # Clustered on vds_key so CreateMappingTables can read one contig at a time.
        gcloud storage ls '~{pairs_parquet_path}/needs_normalization=*/part-*.parquet' |
            sed -E 's#(/part-[0-9]{2})[^/]*$#\1*.parquet#' | sort -u > sources.txt
        wc -l sources.txt

        replace=--replace
        while mapfile -t -n 20 batch && (( ${#batch[@]} > 0 )); do
            bq --apilog=false load --project_id=~{project_id} "${replace}" --source_format=PARQUET \
                --clustering_fields=vds_key ~{dataset_name}.~{pairs_table_name} "$(IFS=,; echo "${batch[*]}")" < /dev/null
            replace=--noreplace
        done < sources.txt

        normalized='~{pairs_parquet_path}/needs_normalization=true'

        # The pairs on keys bcftools norm moved, again, as their own small table: the provenance table reads these
        # rather than scanning the pairs a second time.
        if gcloud storage ls "${normalized}/" > /dev/null 2>&1; then
            bq --apilog=false load --project_id=~{project_id} --replace --source_format=PARQUET \
                ~{dataset_name}.~{normalized_pairs_table_name} "${normalized}/*.parquet"
        else
            bq --apilog=false rm -f --project_id=~{project_id} ~{dataset_name}.~{normalized_pairs_table_name}
            bq --apilog=false mk --project_id=~{project_id} --table \
                ~{dataset_name}.~{normalized_pairs_table_name} vds_key:STRING,person_id:INTEGER,zygosity:STRING
        fi

        bq --apilog=false load --project_id=~{project_id} --replace --source_format=CSV --field_delimiter='\t' \
            ~{dataset_name}.~{key_to_vid_table_name} ~{key_to_vid} vds_key:STRING,vid:STRING

        # Genome-wide the pairs table is around 77 TB, so do not leave it to accrue storage.
        for table in ~{pairs_table_name} ~{normalized_pairs_table_name} ~{key_to_vid_table_name}; do
            bq --apilog=false update --project_id=~{project_id} --expiration=604800 ~{dataset_name}."${table}"
        done
    >>>

    meta {
        # The tables expire, so a relaunch must reload them rather than reuse a cached result naming tables that may
        # be gone.
        volatile: true
    }

    runtime {
        docker: variants_docker
    }

    output {
        String pairs_table = pairs_table_name
        String normalized_pairs_table = normalized_pairs_table_name
        String key_to_vid_table = key_to_vid_table_name
    }
}

task CreateMappingTables {
    input {
        String project_id
        String dataset_name
        String participant_mapping_table_name
        String pairs_table
        String normalized_pairs_table
        String key_to_vid_table
        String variants_docker
    }

    String base_table_name = "~{participant_mapping_table_name}_base"
    String provenance_table_name = "~{participant_mapping_table_name}_provenance"
    String fq_dataset = "~{project_id}.~{dataset_name}"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cat > create_mapping_tables.sql <<'SQL'
        DECLARE unnormalized_keys INT64;
        DECLARE unnormalized_examples ARRAY<STRING>;

        CREATE TEMP FUNCTION vidToLocation(vid STRING)
        RETURNS INT64
        AS (
            (CASE SPLIT(vid, '-')[OFFSET(0)]
                WHEN 'X' THEN 23
                WHEN 'Y' THEN 24
                ELSE CAST(SPLIT(vid, '-')[OFFSET(0)] AS INT64) END) * 1000000000000 +
            CAST(SPLIT(vid, '-')[OFFSET(1)] AS INT64)
        );

        -- Every flagged key with carriers must have a VID from bcftools norm. One without would be filed under its own,
        -- non-VAT key. Reads only the small tables, so a mismatched set of inputs fails before the pairs are scanned.
        SET (unnormalized_keys, unnormalized_examples) = (
            SELECT AS STRUCT COUNT(*), ARRAY_AGG(n.vds_key ORDER BY n.vds_key LIMIT 10)
            FROM (SELECT DISTINCT vds_key FROM `~{fq_dataset}.~{normalized_pairs_table}`) AS n
            LEFT JOIN `~{fq_dataset}.~{key_to_vid_table}` AS k USING (vds_key)
            WHERE k.vid IS NULL
        );
        IF unnormalized_keys > 0 THEN
            RAISE USING MESSAGE = FORMAT('%d keys flagged for normalization have no VID, e.g. %s',
                                         unnormalized_keys, ARRAY_TO_STRING(unnormalized_examples, ', '));
        END IF;

        -- Filled one contig at a time below.
        CREATE TABLE `~{fq_dataset}.~{base_table_name}` (
            vid STRING,
            het_ids ARRAY<INT64>,
            hom_ids ARRAY<INT64>,
            hemi_ids ARRAY<INT64>
        )
        -- Clustered at creation: a bq update after the fact does not take. No ORDER BY; see VS-2013.
        CLUSTER BY vid;

        -- One row per normalized key that has carriers. input_location, input_ref and input_alt join to alt_allele on
        -- (location, ref, allele), which explains a participant filed under a VID whose own location has no alt_allele
        -- row for them.
        CREATE TABLE `~{fq_dataset}.~{provenance_table_name}`
        CLUSTER BY vid
        AS
        SELECT k.vid,
               vidToLocation(n.vds_key) AS input_location,
               SPLIT(n.vds_key, '-')[OFFSET(2)] AS input_ref,
               SPLIT(n.vds_key, '-')[OFFSET(3)] AS input_alt,
               n.n_het, n.n_hom, n.n_hemi
        FROM (
            SELECT vds_key,
                   COUNTIF(zygosity = 'het') AS n_het,
                   COUNTIF(zygosity = 'hom') AS n_hom,
                   COUNTIF(zygosity = 'hemi') AS n_hemi
            FROM `~{fq_dataset}.~{normalized_pairs_table}`
            GROUP BY vds_key
        ) AS n
        JOIN `~{fq_dataset}.~{key_to_vid_table}` AS k USING (vds_key);
        SQL

        # bq query --max_rows check: ok, results go to the new tables
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} "$(cat create_mapping_tables.sql)"

        # The scan of the pairs, one contig at a time: genome-wide, a single query would group some 2.5e12 rows, near
        # BigQuery's six-hour limit. The pairs table is clustered on vds_key, so each query reads only its contig's
        # rows. Normalization never moves a key to another contig, so no VID's carriers span two queries. The list is
        # every contig GVS can hold: its locations encode only chromosomes 1-22, X and Y. Whether every VID here is a
        # VAT VID is checked against the finished table, by mapped_not_in_vat in CheckMappingTables.
        cat > insert_contig.sql <<'SQL'
        INSERT INTO `~{fq_dataset}.~{base_table_name}` (vid, het_ids, hom_ids, hemi_ids)
        WITH person_zygosity AS (
            -- One row per (vid, person). A person appears once per key, so only VIDs that more than one key normalizes
            -- to can repeat a person here. Such a person keeps the strongest zygosity any of their keys gives them,
            -- so no person lands in two of the three arrays.
            SELECT COALESCE(k.vid, p.vds_key) AS vid,
                   p.person_id,
                   MAX(CASE p.zygosity WHEN 'het' THEN 1 WHEN 'hemi' THEN 2 WHEN 'hom' THEN 3 END) AS zygosity_rank
            FROM `~{fq_dataset}.~{pairs_table}` AS p
            LEFT JOIN `~{fq_dataset}.~{key_to_vid_table}` AS k USING (vds_key)
            -- '-' sorts just below '.', so this is exactly the keys on CONTIG, and as constants they prune clusters.
            WHERE p.vds_key >= 'CONTIG-' AND p.vds_key < 'CONTIG.'
            -- By position: an unqualified vid here could resolve to k.vid, which is NULL for nearly every key.
            GROUP BY 1, 2
        )
        SELECT vid,
               -- ARRAY_AGG over nothing but NULLs is NULL, and ARRAY_CONCAT with a NULL is NULL.
               IFNULL(ARRAY_AGG(IF(zygosity_rank = 1, person_id, NULL) IGNORE NULLS), []) AS het_ids,
               IFNULL(ARRAY_AGG(IF(zygosity_rank = 3, person_id, NULL) IGNORE NULLS), []) AS hom_ids,
               IFNULL(ARRAY_AGG(IF(zygosity_rank = 2, person_id, NULL) IGNORE NULLS), []) AS hemi_ids
        FROM person_zygosity
        GROUP BY vid
        SQL

        for contig in $(seq 1 22) X Y; do
            # bq query --max_rows check: ok, results go to the base table
            bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} \
                "$(sed "s/CONTIG/${contig}/g" insert_contig.sql)"
        done

        # A query of its own: BigQuery will not create a view in a script that declares a temporary function, whether
        # or not the view uses it.
        # bq query --max_rows check: ok, creates a view
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} \
            'CREATE VIEW `~{fq_dataset}.~{participant_mapping_table_name}` AS
             SELECT vid, ARRAY_CONCAT(het_ids, hom_ids, hemi_ids) AS person_ids
             FROM `~{fq_dataset}.~{base_table_name}`'
    >>>

    runtime {
        docker: variants_docker
    }

    output {
        String base_table = "~{fq_dataset}.~{base_table_name}"
        String provenance_table = "~{fq_dataset}.~{provenance_table_name}"
        String view = "~{fq_dataset}.~{participant_mapping_table_name}"
    }
}

# The VS-2013 acceptance criteria, against the tables as built. Fails the workflow, leaving the tables in place for
# inspection, if any is not met. One scan of the base table's arrays; everything else reads vid columns.
task CheckMappingTables {
    input {
        String project_id
        String fq_vat_table
        String base_table
        String provenance_table
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cat > acceptance_checks.sql <<'SQL'
        WITH contigs AS (
            SELECT DISTINCT SPLIT(vid, '-')[OFFSET(0)] AS contig FROM `~{base_table}`
        ),
        -- Only the VAT VIDs on contigs the mapping covers, so a run restricted to whole contigs is checked as a whole.
        -- Not the first line of the query: bq would parse a query starting with -- as a flag.
        vat AS (
            SELECT vid, ANY_VALUE(gvs_all_ac) AS ac, ANY_VALUE(gvs_all_sc) AS sc
            FROM `~{fq_vat_table}`
            WHERE SPLIT(vid, '-')[OFFSET(0)] IN (SELECT contig FROM contigs)
            GROUP BY vid
        ),
        mapping AS (
            SELECT vid,
                   ARRAY_LENGTH(het_ids) + ARRAY_LENGTH(hom_ids) + ARRAY_LENGTH(hemi_ids) AS n,
                   ARRAY_LENGTH(hom_ids) AS n_hom,
                   (SELECT COUNT(DISTINCT person_id)
                    FROM UNNEST(ARRAY_CONCAT(het_ids, hom_ids, hemi_ids)) AS person_id) AS n_distinct
            FROM `~{base_table}`
        ),
        provenance AS (
            SELECT DISTINCT vid FROM `~{provenance_table}`
        ),
        j AS (
            SELECT m.vid IS NOT NULL AS in_mapping,
                   v.vid IS NOT NULL AS in_vat,
                   p.vid IS NOT NULL AS in_provenance,
                   m.n, m.n_hom, m.n_distinct, v.ac, v.sc
            FROM mapping AS m
            FULL JOIN vat AS v ON v.vid = m.vid
            LEFT JOIN provenance AS p ON p.vid = COALESCE(m.vid, v.vid)
        )
        SELECT
            COUNTIF(in_mapping) AS mapped_vids,
            COUNTIF(in_provenance) AS provenance_vids,
            -- Criterion 1. Outside the provenance VIDs the mapping count is exactly gvs_all_sc; inside, where synonym
            -- carriers are unioned in, it may exceed it but never fall below.
            COUNTIF(in_mapping AND in_vat AND NOT in_provenance AND n = sc) AS exact_outside_provenance,
            COUNTIF(in_mapping AND in_vat AND NOT in_provenance AND n != sc) AS mismatch_outside_provenance,
            COUNTIF(in_provenance AND n > sc) AS above_inside_provenance,
            COUNTIF(in_provenance AND n < sc) AS below_inside_provenance,
            -- The split. The VAT has no zygosity columns, but a hom carrier adds 2 to gvs_all_ac and a het or hemi
            -- carrier 1, so hom = ac - sc. With the total exact, that also fixes het + hemi; het and hemi are not
            -- separable from the VAT. Inside the provenance VIDs a person on two keys keeps the stronger zygosity,
            -- so there is nothing exact to compare.
            COUNTIF(in_mapping AND in_vat AND NOT in_provenance AND n_hom != ac - sc) AS hom_mismatch_outside_provenance,
            -- Criterion 2. Every mapped VID is a VAT VID, and every VAT VID with carriers has a mapping row.
            COUNTIF(in_mapping AND NOT in_vat) AS mapped_not_in_vat,
            COUNTIF(in_vat AND NOT in_mapping AND sc > 0) AS vat_carriers_without_row,
            -- Every provenance VID has a mapping row. Counted from provenance, as j drops a VID in neither table.
            (SELECT COUNT(*) FROM provenance WHERE vid NOT IN (SELECT vid FROM mapping)) AS provenance_not_mapped,
            -- Criterion 3. No person appears twice for one VID, within an array or across them.
            COUNTIF(n != n_distinct) AS vids_with_duplicate_person
        FROM j
        SQL

        # bq query --max_rows check: ok, a single row of counts
        bq --apilog=false query --nouse_legacy_sql --project_id=~{project_id} --format=json \
            "$(cat acceptance_checks.sql)" | jq '.[0] | map_values(tonumber)' > acceptance_checks.json
        cat acceptance_checks.json

        failures=$(jq -r 'to_entries[]
            | select(.key | IN("mismatch_outside_provenance", "below_inside_provenance", "hom_mismatch_outside_provenance",
                               "mapped_not_in_vat",
                               "vat_carriers_without_row", "provenance_not_mapped", "vids_with_duplicate_person"))
            | select(.value != 0) | "\(.key)=\(.value)"' acceptance_checks.json)
        if [[ -n "${failures}" ]]
        then
            echo "Acceptance checks failed: ${failures//$'\n'/, }" >&2
            exit 1
        fi
    >>>

    runtime {
        docker: variants_docker
    }

    output {
        File acceptance_checks = "acceptance_checks.json"
    }
}
