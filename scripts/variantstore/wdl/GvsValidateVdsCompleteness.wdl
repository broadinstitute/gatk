version 1.0

# Screen a GVS VDS for "rectangle" data dropouts: a contiguous genomic window in which one
# GVS superpartition has little or no data while every other superpartition has the usual
# amount. This is the shape produced when a single Avro export shard is lost, truncated, or
# never read, because the `EXPORT DATA ... ORDER BY location` in GvsExtractAvroFilesForHail.wdl
# means each numbered output file holds a contiguous location range for exactly one
# superpartition.
#
# The work is split into actions rather than run end to end, because the cost is very unevenly
# distributed. `materialize` pays one full-width read of the callset to write a downsampled
# narrow copy, and is the single large bill. Every other action runs against that narrow copy
# for a small fraction of the cost, so calibration and threshold tuning are cheap and
# repeatable. Run one action per invocation.
#
# Typical sequence for a callset:
#   1. action = "materialize"  on the source VDS      -> narrow copy + sample list
#   2. action = "probe"        on the narrow copy     -> measured throughput, cost estimate
#   3. action = "scan"         on the narrow copy     -> summary + candidate rectangles + SQL
#   4. action = "full-depth"   on the source VDS      -> per-sample detail for flagged windows
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
        String? sample_list_path
        String? narrow_vds_path

        Int samples_per_superpartition = 100
        String random_seed = "vs-1998"
        Int superpartition_size = 4000
        Int bin_size = 50000
        Int stride = 1

        String? contigs
        String? intervals
        String? target_superpartitions
        String probe_interval = "chr20:30000000-40000000"
        Boolean overwrite = false

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
        Int num_primary_workers = 4
        Int max_secondary_workers = 300
        String cluster_prefix = "vds-completeness"
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
            help: "One of materialize, scan, probe, full-depth. See the header comment for the usual sequence."
        }
        vds_path: {
            help: "VDS to read. The source VDS for materialize and full-depth; normally the narrow copy for scan and probe."
        }
        output_prefix: {
            help: "GCS prefix under which outputs are written (narrow VDS, sample list, summary, report, SQL)."
        }
        mode: {
            help: "variants counts variant_data entries; references sums reference-block coverage."
        }
        sample_map_path: {
            help: "TSV of sample_name and sample_id from sample_info. Required for materialize."
        }
        sample_list_path: {
            help: "Chosen-sample TSV written by materialize. Preferred by later actions so every VDS is screened on exactly the same samples, which is what makes figures comparable across classic, r2 and r3."
        }
        narrow_vds_path: {
            help: "Destination for the narrow copy when action is materialize. Defaults to <output_prefix>/narrow.vds."
        }
        samples_per_superpartition: {
            help: "Stratified sampling depth. 100 retains essentially all detection power for shard-scale dropouts at a fortieth of the columns; it cannot see dropouts affecting only a handful of individual samples."
        }
        random_seed: {
            help: "Seed for sample selection. Keep fixed across VDSes or cross-VDS comparisons become meaningless."
        }
        stride: {
            help: "Read one genomic bin in every N. Default 1 reads everything. Larger values cost roughly 1/N and guarantee detection of any dropout wider than stride * bin_size, but say nothing about unscanned bins."
        }
        intervals: {
            help: "Comma-separated Hail locus intervals. Overrides contigs. Required for full-depth."
        }
        bq_project_id: {
            help: "BigQuery project used only to generate adjudication SQL. No queries are run by this workflow."
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

    if (action != "materialize" && action != "scan" && action != "probe" && action != "full-depth") {
        call Utils.TerminateWorkflow as UnrecognizedAction {
            input:
                message = "`action` must be one of materialize, scan, probe, full-depth; got '" + action + "'.",
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

    # materialize is the only action that derives the sample selection from scratch; every other
    # action should reuse the recorded list so results stay comparable between VDSes.
    if (action == "materialize" && !defined(sample_map_path)) {
        call Utils.TerminateWorkflow as MaterializeNeedsSampleMap {
            input:
                message = "`sample_map_path` is required when action is materialize.",
                basic_docker = effective_basic_docker,
        }
    }

    if (action != "materialize" && !defined(sample_map_path) && !defined(sample_list_path)) {
        call Utils.TerminateWorkflow as ScanNeedsSamples {
            input:
                message = "action '" + action + "' requires either `sample_list_path` (preferred) or `sample_map_path`.",
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
            sample_map_path = sample_map_path,
            sample_list_path = sample_list_path,
            narrow_vds_path = narrow_vds_path,
            samples_per_superpartition = samples_per_superpartition,
            random_seed = random_seed,
            superpartition_size = superpartition_size,
            bin_size = bin_size,
            stride = stride,
            contigs = contigs,
            intervals = intervals,
            target_superpartitions = target_superpartitions,
            probe_interval = probe_interval,
            overwrite = overwrite,
            ratio_threshold = ratio_threshold,
            score_threshold = score_threshold,
            min_expected = min_expected,
            scale_threshold = scale_threshold,
            baseline_quantile = baseline_quantile,
            bq_project_id = bq_project_id,
            bq_dataset_name = bq_dataset_name,
            reference_schema = effective_reference_schema,
            prefix = cluster_prefix,
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
        String cluster_name = ScanVdsForDropouts.cluster_name
        String scan_log = ScanVdsForDropouts.scan_log
        File report = ScanVdsForDropouts.report
        File adjudication_sql = ScanVdsForDropouts.adjudication_sql
        Boolean done = true
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
        String? sample_list_path
        String? narrow_vds_path

        Int samples_per_superpartition
        String random_seed
        Int superpartition_size
        Int bin_size
        Int stride

        String? contigs
        String? intervals
        String? target_superpartitions
        String probe_interval
        Boolean overwrite

        Float? ratio_threshold
        Float? score_threshold
        Float? min_expected
        Float? scale_threshold
        Float? baseline_quantile

        String? bq_project_id
        String? bq_dataset_name
        String reference_schema

        String prefix
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

    String effective_narrow_vds_path = select_first([narrow_vds_path, output_prefix + "/narrow.vds"])
    String summary_path = output_prefix + "/summary_" + mode + ".tsv"
    String superpartitions_path = output_prefix + "/superpartitions_" + mode + ".tsv"
    String effective_sample_list_path = select_first([sample_list_path, output_prefix + "/sample_list.tsv"])
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
            "samples-per-superpartition": ~{samples_per_superpartition},
            "random-seed": "~{random_seed}",
            "superpartition-size": ~{superpartition_size},
            "bin-size": ~{bin_size},
            "stride": ~{stride},
            "probe-interval": "~{probe_interval}",
            "temp-path": sys.argv[1],
        }

        if action == "materialize":
            arguments["output-path"] = "~{effective_narrow_vds_path}"
            arguments["sample-list-path"] = "~{effective_sample_list_path}"
            arguments["superpartitions-path"] = "~{superpartitions_path}"
            arguments["overwrite"] = "~{true='true' false='false' overwrite}"
        elif action == "scan":
            arguments["summary-path"] = "~{summary_path}"
            arguments["superpartitions-path"] = "~{superpartitions_path}"
        elif action == "full-depth":
            arguments["full-depth-path"] = "~{full_depth_path}"

        # A recorded sample list is preferred wherever one exists, so that every VDS is
        # screened on identical samples; the raw map is the fallback for materialize.
        for key, value in [
            ("sample-map-path", "~{default='' sample_map_path}"),
            ("sample-list-path", "~{default='' sample_list_path}"),
            ("contigs", "~{default='' contigs}"),
            ("intervals", "~{default='' intervals}"),
            ("target-superpartitions", "~{default='' target_superpartitions}"),
        ]:
            if value:
                arguments[key] = value

        print(json.dumps(arguments, indent=2))
        PYTHON_HEREDOC

        cat script-arguments.json

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
