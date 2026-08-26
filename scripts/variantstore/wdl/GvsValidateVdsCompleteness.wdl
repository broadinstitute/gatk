version 1.0

# Screen a GVS VDS for "rectangle" data dropouts: a contiguous genomic window in which one
# GVS superpartition has little or no data while every other superpartition has the usual
# amount. This is the shape produced when a single Avro export shard is lost, truncated, or
# never read, because the `EXPORT DATA ... ORDER BY location` in GvsExtractAvroFilesForHail.wdl
# means each numbered output file holds a contiguous location range for exactly one
# superpartition.
#
# Every sample is screened. There is no sampling, of samples or of loci: measured on Foxtrot
# r2 (535,662 samples, 119,189 variant_data partitions), a genome-wide variant scan runs in
# about an hour and a reference scan in a couple of hours at full autoscaling width, which
# leaves nothing worth buying by screening a subset. An earlier design materialized a
# downsampled copy first; that only made sense while a full-width pass was assumed expensive.
#
# Run one action per invocation. Typical sequence for a callset:
#   1. action = "scan", mode = "variants"    -> summary, candidate rectangles, adjudication SQL
#   2. action = "scan", mode = "references"  -> the same for reference coverage
#   3. action = "full-depth"                 -> per-sample detail for whatever step 1 or 2 flagged
#
# `action = "probe"` is optional, for sizing an unfamiliar callset before committing to a full
# scan. Give it a bounded interval wide enough to span many partitions: a partition is the unit
# of work, so an interval resolving to a handful of them measures near-serial throughput and
# takes about as long as a fully parallel genome-wide scan.
#
# The sample map is generated automatically when bq_project_id and bq_dataset_name are supplied.
# Pass sample_map_path to supply one yourself, or see
# scripts/variantstore/bq/vds_dropout_sample_map.sql to build it by hand.
#
# The BigQuery adjudication SQL that `scan` emits is deliberately not executed here. The screen
# produces candidates; running the queries against the callset dataset is a separate, deliberate
# step.

import "GvsUtils.wdl" as Utils


workflow GvsValidateVdsCompleteness {
    input {
        # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go = true

        String action = "scan"
        String vds_path
        String output_prefix
        String mode = "variants"

        String? sample_map_path
        # Supply bq_project_id and bq_dataset_name and the sample map is generated for you;
        # the query is a two-line projection of sample_info.
        String bq_sample_table = "sample_info"

        Int superpartition_size = 4000
        Int bin_size = 50000

        String? contigs
        String? intervals
        String? target_superpartitions
        # Fully callable 10 Mb window on chr20, clear of the centromere gap
        # (chr20:26,386,232-30,088,349), the telomere, and the smaller gaps near 31 Mb.
        String probe_interval = "chr20:37000000-47000000"

        # Detection thresholds. Cheap to revisit: the detector runs against the small summary
        # table, so re-judging an existing scan costs nothing.
        Float? ratio_threshold
        Float? score_threshold
        Float? min_expected
        Float? scale_threshold
        Float? baseline_quantile

        # Adjudication SQL generation. Omitting the project or dataset skips it.
        String? bq_project_id
        String? bq_dataset_name
        # Which ref_ranges schema the dataset uses. Normally left unset: it is detected from
        # the dataset with Utils.IsUsingCompressedReferences, the same way
        # GvsExtractAvroFilesForHail.wdl decides how to export. Set it only to override that.
        String? reference_schema

        Boolean use_tiny_dataproc_cluster = false
        # Two rather than four. Primary workers must all be provisioned before the cluster
        # is usable, and Dataproc's floor is 2 datanodes, so asking for 4 doubles the
        # up-front capacity ask for no benefit -- a us-central1-b create failed having
        # provisioned only 1 of 4. Secondary workers are added by the autoscaling policy
        # afterwards and are preemptible, so peak parallelism is unaffected.
        Int num_primary_workers = 2
        Int max_secondary_workers = 300
        String cluster_prefix = "vds-completeness"
        # Compute Engine stockouts are per-zone, and Dataproc's default placement picks one
        # zone and fails rather than moving on. "auto" tries every zone in the region.
        String cluster_zones = "auto"
        # Cores per worker, used only for the probe's concurrency estimate.
        Int worker_cores = 4
        # Kept at 1, matching every other GVS Hail workflow. In principle this screen is a
        # single streaming pass with nothing to spill, and dropping local SSD widens the
        # zones able to serve a request -- but that is now handled by cluster_zones retrying
        # a stockout elsewhere, so there is no reason to also diverge from the configuration
        # that is known to work. Worth revisiting for cost once the pipeline is proven.
        Int num_local_ssds = 1
        String worker_machine_type = "n1-highmem-8"
        # Also left matching the other GVS workflows. The driver only collects a summary on
        # the order of 62,000 x 134 numbers, so this is oversized on paper, but a smaller
        # master is a variable worth removing while the pipeline is still being brought up.
        String master_machine_type = "n1-highmem-32"
        String? hail_temp_path
        String region = "us-central1"
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Boolean leave_cluster_running_at_end = false
        Float? master_memory_fraction

        String? git_branch_or_tag
        String? hail_version
        File? hail_wheel
        String? basic_docker
        String? variants_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? workspace_bucket
        String? workspace_project
    }

    parameter_meta {
        action: {
            help: "One of scan, probe, full-depth. See the header comment for the usual sequence."
        }
        vds_path: {
            help: "VDS to read."
        }
        output_prefix: {
            help: "GCS prefix under which outputs are written (sample map, summary, report, adjudication SQL, logs)."
        }
        mode: {
            help: "variants counts variant_data entries; references sums reference-block coverage."
        }
        sample_map_path: {
            help: "TSV of sample_name and sample_id. Optional: if omitted and bq_project_id and bq_dataset_name are set, it is generated from sample_info automatically."
        }
        bq_sample_table: {
            help: "Sample table or view the generated map reads from. Defaults to sample_info; use a view such as sample_info_new_to_foxtrot to restrict the sample universe."
        }
        intervals: {
            help: "Comma-separated Hail locus intervals. Overrides contigs. Required for full-depth."
        }
        bq_project_id: {
            help: "BigQuery project used only to generate adjudication SQL. No queries are run by this workflow."
        }
        cluster_zones: {
            help: "Comma-separated zones to try in order, or auto for every zone in the region. A stockout is a per-zone condition, so a create that fails for capacity is retried in the next zone rather than failing the workflow."
        }
        num_local_ssds: {
            help: "Local SSDs per master and worker. Kept at 1 to match the other GVS Hail workflows; lowering it widens the zones able to serve a request, at the cost of diverging from a known-good configuration."
        }
        probe_interval: {
            help: "Interval used when action is probe. Defaults to a fully callable 10 Mb window on chr20; pick a gap-free interval, since one overlapping a centromere reads faster than average and skews the extrapolation."
        }
        reference_schema: {
            help: "Override for the ref_ranges schema, either compressed or uncompressed. Leave unset to detect it from the dataset. AoU callsets are compressed, where reference adjudication filters on packed_ref_data, the clustering field."
        }
    }

    if (!defined(variants_docker) || !defined(basic_docker) || !defined(cloud_sdk_slim_docker) ||
        !defined(cloud_sdk_docker) || !defined(workspace_bucket) || !defined(workspace_project) ||
        !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_google_project = select_first([workspace_project, GetToolVersions.google_project])
    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])

    if (defined(hail_version) && defined(hail_wheel)) {
        call Utils.TerminateWorkflow as BothHailVersionAndHailWheelDefined {
            input:
                message = "Cannot define both `hail_version` and `hail_wheel`, exiting.",
                basic_docker = effective_basic_docker,
        }
    }

    if (action != "scan" && action != "probe" && action != "full-depth") {
        call Utils.TerminateWorkflow as UnrecognizedAction {
            input:
                message = "`action` must be one of scan, probe, full-depth; got '" + action + "'.",
                basic_docker = effective_basic_docker,
        }
    }

    if (defined(reference_schema) && reference_schema != "compressed" && reference_schema != "uncompressed") {
        call Utils.TerminateWorkflow as UnrecognizedReferenceSchema {
            input:
                message = "`reference_schema` must be either compressed or uncompressed; got '" + select_first([reference_schema, ""]) + "'.",
                basic_docker = effective_basic_docker,
        }
    }

    if (mode != "variants" && mode != "references") {
        call Utils.TerminateWorkflow as UnrecognizedMode {
            input:
                message = "`mode` must be either variants or references; got '" + mode + "'.",
                basic_docker = effective_basic_docker,
        }
    }

    # Every action needs the sample map: superpartition membership is a function of
    # sample_id, which a VDS does not carry.
    Boolean can_generate_sample_map = defined(bq_project_id) && defined(bq_dataset_name)

    if (!defined(sample_map_path) && !can_generate_sample_map) {
        call Utils.TerminateWorkflow as NeedsSampleMap {
            input:
                message = "Provide `sample_map_path`, or `bq_project_id` and `bq_dataset_name` so the map can be generated from sample_info.",
                basic_docker = effective_basic_docker,
        }
    }

    if (action == "full-depth" && !defined(intervals)) {
        call Utils.TerminateWorkflow as FullDepthNeedsIntervals {
            input:
                message = "`intervals` is required when action is full-depth.",
                basic_docker = effective_basic_docker,
        }
    }

    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    # The ref_ranges schema is a property of the dataset, not a choice, so detect it rather
    # than trusting an input. Only needed when reference adjudication SQL will be generated,
    # which is the only thing the schema affects.
    if (mode == "references" && !defined(reference_schema) &&
        defined(bq_project_id) && defined(bq_dataset_name)) {
        call Utils.GetBQTableLastModifiedDatetime as RefTableDatetimeCheck {
            input:
                project_id = select_first([bq_project_id]),
                fq_table = select_first([bq_project_id]) + "." + select_first([bq_dataset_name]) + ".ref_ranges_001",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call Utils.IsUsingCompressedReferences {
            input:
                query_project_id = select_first([bq_project_id]),
                dest_project_id = select_first([bq_project_id]),
                dataset_name = select_first([bq_dataset_name]),
                ref_table_timestamp = RefTableDatetimeCheck.last_modified_timestamp,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }
    }

    # Explicit override wins; otherwise use what was detected. The fallback only applies when
    # no detection ran, in which case no reference SQL is generated and the value is unused.
    String effective_reference_schema =
        if defined(reference_schema) then select_first([reference_schema])
        else if defined(IsUsingCompressedReferences.is_using_compressed_references) then
            (if select_first([IsUsingCompressedReferences.is_using_compressed_references])
             then "compressed" else "uncompressed")
        else "compressed"

    if (!defined(sample_map_path) && can_generate_sample_map) {
        call GenerateSampleMap {
            input:
                bq_project_id = select_first([bq_project_id]),
                bq_dataset_name = select_first([bq_dataset_name]),
                sample_table = bq_sample_table,
                output_prefix = output_prefix,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }
    }

    String? effective_sample_map_path =
        if defined(sample_map_path) then sample_map_path else GenerateSampleMap.sample_map_path

    call Utils.GetHailScripts {
        input:
            variants_docker = effective_variants_docker,
    }

    call ScanVdsForDropouts {
        input:
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
            vds_dropout_scan_script = GetHailScripts.vds_dropout_scan_script,
            vds_dropout_detect_script = GetHailScripts.vds_dropout_detect_script,
            action = action,
            vds_path = vds_path,
            output_prefix = output_prefix,
            mode = mode,
            sample_map_path = effective_sample_map_path,
            superpartition_size = superpartition_size,
            bin_size = bin_size,
            contigs = contigs,
            intervals = intervals,
            target_superpartitions = target_superpartitions,
            probe_interval = probe_interval,
            ratio_threshold = ratio_threshold,
            score_threshold = score_threshold,
            min_expected = min_expected,
            scale_threshold = scale_threshold,
            baseline_quantile = baseline_quantile,
            bq_project_id = bq_project_id,
            bq_dataset_name = bq_dataset_name,
            reference_schema = effective_reference_schema,
            prefix = cluster_prefix,
            cluster_zones = cluster_zones,
            worker_cores = worker_cores,
            num_local_ssds = num_local_ssds,
            worker_machine_type = worker_machine_type,
            master_machine_type = master_machine_type,
            use_tiny_dataproc_cluster = use_tiny_dataproc_cluster,
            num_primary_workers = num_primary_workers,
            max_secondary_workers = max_secondary_workers,
            hail_version = effective_hail_version,
            hail_wheel = hail_wheel,
            hail_temp_path = hail_temp_path,
            workspace_project = effective_google_project,
            workspace_bucket = effective_workspace_bucket,
            region = region,
            leave_cluster_running_at_end = leave_cluster_running_at_end,
            cluster_max_idle_minutes = cluster_max_idle_minutes,
            cluster_max_age_minutes = cluster_max_age_minutes,
            master_memory_fraction = master_memory_fraction,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
    }

    output {
        String? generated_sample_map = GenerateSampleMap.sample_map_path
        String cluster_name = ScanVdsForDropouts.cluster_name
        String scan_log = ScanVdsForDropouts.scan_log
        File report = ScanVdsForDropouts.report
        File adjudication_sql = ScanVdsForDropouts.adjudication_sql
        Boolean done = true
    }
}


task GenerateSampleMap {
    meta {
        description: "Export sample_name and sample_id from sample_info as the --sample-map-path input."
        # Not cached: sample_info changes as samples are withdrawn, and a stale map would
        # be silently wrong. Regenerating costs seconds. A map that has fallen behind the
        # VDS is caught loudly anyway -- vds_dropout_scan.py fails on a VDS sample it
        # cannot place in a superpartition rather than screening a biased subset.
        volatile: true
    }

    input {
        String bq_project_id
        String bq_dataset_name
        String sample_table
        String output_prefix
        String cloud_sdk_docker
    }

    parameter_meta {
        sample_table: {
            help: "Table or view to read. sample_info for the whole callset; a view to restrict the sample universe."
        }
    }

    String map_prefix = output_prefix + "/samples/sample_map"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # The same filter GvsExtractAvroFilesForHail.wdl applies when exporting Avro, so the
        # sample universe matches the one the VDS was built from. Kept in step with
        # scripts/variantstore/bq/vds_dropout_sample_map.sql, which is the hand-run copy.
        #
        # EXPORT DATA rather than `bq query > file`: it writes tab-delimited with a header
        # straight to GCS, which is the format the scan expects, and it does not depend on
        # paging half a million rows through the CLI.
        bq query \
            --project_id=~{bq_project_id} \
            --use_legacy_sql=false \
            --nouse_cache \
            --format=none \
            "EXPORT DATA OPTIONS(
                 uri='~{map_prefix}_*.tsv',
                 format='CSV',
                 field_delimiter='\t',
                 header=true,
                 overwrite=true) AS
             SELECT sample_name, sample_id
             FROM \`~{bq_project_id}.~{bq_dataset_name}.~{sample_table}\`
             WHERE withdrawn IS NULL
               AND is_control = false
             ORDER BY sample_id"

        # EXPORT DATA requires a wildcard and shards above 1GB. A sample map is a few MB, so
        # exactly one file is expected; more than one means the downstream reader would
        # silently see only part of the callset, so fail rather than guess.
        gsutil ls "~{map_prefix}_*.tsv" > exported_files.txt
        file_count=$(wc -l < exported_files.txt)
        if [[ "${file_count}" -ne 1 ]]
        then
            echo "Expected exactly one exported sample map file, found ${file_count}:" >&2
            cat exported_files.txt >&2
            exit 1
        fi
        tr -d '\n' < exported_files.txt > sample_map_path.txt

        # Report the row count so a suspiciously small map is obvious in the log.
        gsutil cat "$(cat sample_map_path.txt)" | tail -n +2 | wc -l > sample_count.txt
        echo "Exported $(cat sample_count.txt) samples to $(cat sample_map_path.txt)"
    >>>

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        disks: "local-disk 50 HDD"
        cpu: 1
        preemptible: 0
    }

    output {
        String sample_map_path = read_string("sample_map_path.txt")
        Int sample_count = read_int("sample_count.txt")
    }
}


task ScanVdsForDropouts {
    input {
        File run_in_hail_cluster_script
        File vds_dropout_scan_script
        File vds_dropout_detect_script

        String action
        String vds_path
        String output_prefix
        String mode

        String? sample_map_path

        Int superpartition_size
        Int bin_size

        String? contigs
        String? intervals
        String? target_superpartitions
        String probe_interval

        Float? ratio_threshold
        Float? score_threshold
        Float? min_expected
        Float? scale_threshold
        Float? baseline_quantile

        String? bq_project_id
        String? bq_dataset_name
        String reference_schema

        String prefix
        String cluster_zones
        Int worker_cores
        Int num_local_ssds
        String worker_machine_type
        String master_machine_type
        Boolean use_tiny_dataproc_cluster
        Int num_primary_workers
        Int max_secondary_workers
        String? hail_version
        File? hail_wheel
        String? hail_temp_path
        String workspace_project
        String workspace_bucket
        String region
        Boolean leave_cluster_running_at_end
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction
        String cloud_sdk_slim_docker
    }

    meta {
        # Always re-run: the point is to observe the current state of a VDS.
        volatile: true
    }

    String summary_path = output_prefix + "/summary_" + mode + ".tsv"
    String superpartitions_path = output_prefix + "/superpartitions_" + mode + ".tsv"
    String full_depth_path = output_prefix + "/full_depth_" + mode + ".tsv"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        apt-get update
        apt install --assume-yes python3.11-venv
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

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        if [[ -z "~{hail_temp_path}" ]]
        then
            hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"
        else
            hail_temp_path="~{hail_temp_path}"
        fi

        # Declare the report outputs up front so they always exist. Cromwell resolves task
        # outputs regardless of which branch ran, and only the scan action produces these.
        echo "action '~{action}' does not produce a report; see the scan action." > report.tsv
        cp report.tsv adjudicate.sql

        # Build the arguments JSON for the script that will run inside the Hail cluster.
        # run_in_hail_cluster.py renders each key as `--key value`, so every key must be
        # kebab-case and must carry a value; boolean flags are spelled out explicitly and
        # optional keys are omitted entirely rather than passed empty.
        python3 - "${hail_temp_path}" > script-arguments.json <<'PYTHON_HEREDOC'
        import json
        import sys

        action = "~{action}"

        arguments = {
            "action": action,
            "vds-path": "~{vds_path}",
            "mode": "~{mode}",
            "superpartition-size": ~{superpartition_size},
            "bin-size": ~{bin_size},
            "probe-interval": "~{probe_interval}",
            # Only used to annotate the probe report with the concurrency a full run could
            # reach; without them the report can report throughput but not project from it.
            "max-secondary-workers": ~{max_secondary_workers},
            "worker-cores": ~{worker_cores},
            "temp-path": sys.argv[1],
        }

        if action == "scan":
            arguments["summary-path"] = "~{summary_path}"
            arguments["superpartitions-path"] = "~{superpartitions_path}"
        elif action == "full-depth":
            arguments["full-depth-path"] = "~{full_depth_path}"

        for key, value in [
            ("sample-map-path", "~{default='' sample_map_path}"),
            ("contigs", "~{default='' contigs}"),
            ("intervals", "~{default='' intervals}"),
            ("target-superpartitions", "~{default='' target_superpartitions}"),
        ]:
            if value:
                arguments[key] = value

        print(json.dumps(arguments, indent=2))
        PYTHON_HEREDOC

        cat script-arguments.json

        # GetHailScripts copies /app/*.py out of the Variants image, so the helper running here
        # is whatever was baked into the current variants_docker tag, not what is in the repo.
        # Check for the flags this workflow relies on before invoking it: without this the
        # failure is an argparse usage dump listing the old arguments, which says nothing about
        # the actual cause.
        for required_flag in --zones --num-local-ssds
        do
            if ! python3 ~{run_in_hail_cluster_script} --help 2>&1 | grep -q -- "${required_flag}"
            then
                echo "ERROR: run_in_hail_cluster.py in this variants_docker image does not support" >&2
                echo "  ${required_flag}, which this workflow passes." >&2
                echo "Rebuild the Variants image and point GetToolVersions at the new tag:" >&2
                echo "  1. scripts/variantstore/scripts/build_docker.sh" >&2
                echo "  2. set variants_docker in GvsUtils.wdl GetToolVersions to the printed tag" >&2
                exit 1
            fi
        done

        # vds_dropout_detect.py rides along as a secondary file so the Hail job can import it,
        # and so the same module that judges the summary on the cluster is the one CI tested.
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{vds_dropout_scan_script} \
            --secondary-script-path-list ~{vds_dropout_detect_script} \
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
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end} \
            2>&1 | tee scan.log

        # Judge the summary locally. This is deliberately not done on the cluster: the summary
        # is small, the judging is pure Python, and running it here means thresholds can be
        # re-applied to an existing scan without another cluster.
        if [[ "~{action}" == "scan" ]]
        then
            gsutil cp "~{summary_path}" ./summary.tsv
            gsutil cp "~{superpartitions_path}" ./superpartitions.tsv

            detect_args=(
                --summary ./summary.tsv
                --superpartitions ./superpartitions.tsv
                --mode ~{mode}
                --report-path ./report.tsv
            )
            ~{'detect_args+=(--ratio-threshold ' + ratio_threshold + ')'}
            ~{'detect_args+=(--score-threshold ' + score_threshold + ')'}
            ~{'detect_args+=(--min-expected ' + min_expected + ')'}
            ~{'detect_args+=(--scale-threshold ' + scale_threshold + ')'}
            ~{'detect_args+=(--baseline-quantile ' + baseline_quantile + ')'}

            if [[ -n "~{default='' bq_project_id}" && -n "~{default='' bq_dataset_name}" ]]
            then
                detect_args+=(--sql-path ./adjudicate.sql)
                detect_args+=(--project-id "~{default='' bq_project_id}")
                detect_args+=(--dataset-name "~{default='' bq_dataset_name}")
                detect_args+=(--reference-schema ~{reference_schema})
            fi

            python3 ~{vds_dropout_detect_script} "${detect_args[@]}" | tee -a scan.log

            gsutil cp ./report.tsv "~{output_prefix}/report_~{mode}.tsv"
            if [[ -f ./adjudicate.sql ]]
            then
                gsutil cp ./adjudicate.sql "~{output_prefix}/adjudicate_~{mode}.sql"
            fi
        fi

        gsutil cp scan.log "~{output_prefix}/scan_~{action}_~{mode}.log"
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
        String cluster_name = read_string("cluster_name.txt")
        String scan_log = output_prefix + "/scan_" + action + "_" + mode + ".log"
        # Always written, so these resolve for every action; see the command body.
        File report = "report.tsv"
        File adjudication_sql = "adjudicate.sql"
    }
}
