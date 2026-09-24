version 1.0

# Screen a GVS VDS for "rectangle" data dropouts: a contiguous genomic window in which one
# GVS superpartition has little or no data while every other superpartition has the usual
# amount. This is the shape produced when a single Avro export shard is lost, truncated, or
# never read, because the `EXPORT DATA ... ORDER BY location` in GvsExtractAvroFilesForHail.wdl
# means each numbered output file holds a contiguous location range for exactly one
# superpartition.
#
# Every sample is screened. There is no sampling, of samples or of loci. Measured on Foxtrot
# r2 (535,662 samples, 119,189 variant_data partitions) at full autoscaling width, the Hail
# aggregation took 5 h 01 m for variants and about 10 hours for references. The shard merge adds
# about 8 minutes at the 10 kb default and the tasks either side of it about twenty, so budget
# roughly 5.5 h and 10.5 h end to end. That is an overnight job whether or not a subset is screened, so
# sampling would not change how this is run -- while costing the exhaustive answer it exists to
# give. Size a run from those figures, not from the ~1 h and ~2 h this project first projected
# off a single-contig probe; those were low by roughly 5-6x.
#
# Which VDSes this applies to. The screen judges each superpartition against its peers, so it
# needs enough of them to have peers at all. Superpartitions hold 4,000 samples each
# (floorDiv(sample_id - 1, 4000) + 1), which puts the useful floor at more than 8,000 samples
# and the recommended width at more than 20,000. Below that, sensitivity to a partially
# depleted window falls off, and a single-superpartition VDS -- any callset of 4,000 samples
# or fewer, which today means everything except AoU -- cannot be screened this way at all:
# the baseline would be the superpartition's own rate, so nothing could ever be flagged.
# vds_dropout_detect.py refuses that case rather than reporting it clean. See
# MIN_SUPERPARTITIONS there for the arithmetic.
#
# What a clean run does not establish. This screens for one shape: data present in the other
# superpartitions and missing from one of them over a contiguous window. It counts variant_data
# entries and reference-block coverage, and never looks at filters, scores, globals or allele
# representation -- so a scan that flags nothing says the data is present at full width, not that
# it is correct. Score and AC/AN/AF correctness belong to the VDS tieout (vds_validation.validate
# and the rescoring in merge_and_rescore_vdses.py), which this does not substitute for.
#
# Run one action per invocation. Typical sequence for a callset:
#   1. action = "scan", mode = "variants"    -> summary, candidate rectangles, adjudication SQL
#   2. action = "scan", mode = "references"  -> the same for reference coverage
#   3. action = "full-depth"                 -> per-sample detail for whatever step 1 or 2 flagged
#
# Step 3 is diagnostic and optional; steps 1 and 2 are what detect a dropout. Where it does earn
# its keep is after a repair that rebuilt a superpartition's samples: the aggregate cannot see a
# handful of samples that survive as columns but carry no data, since one sample is a few
# ten-thousandths of a cell and moves that superpartition's median by the same factor. Such a
# failure would be genome-wide for the affected sample, so one narrow interval with
# target_superpartitions set to the rebuilt superpartitions, checking that n_zero is zero,
# closes that gap in minutes.
#
# Contigs are checkpointed, so an interrupted scan resumes rather than restarting, and the first
# contig to finish reports its own partition count and duration -- which is how a run is sized.
# Partition geometry itself is metadata: `hl.vds.read_vds(path).variant_data.n_partitions()`.
# To rehearse the pipeline cheaply, run a scan over a bounded `intervals` value; it produces a
# real summary the detector can consume.
#
# The sample map is generated automatically when bq_project_id and bq_dataset_name are supplied.
# Pass sample_map_path to supply one yourself; the format is two columns, sample_name and
# sample_id, with a header, and the GenerateSampleMap task below shows the query.
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
        # 10 kb rather than the 50 kb the VS-1946 analysis used. Bin size fixes the
        # detection floor -- clearing the 0.5 ratio gate needs more than half a bin
        # depleted, so 50 kb only sees a gap wider than ~25 kb -- and it is baked into the
        # Hail pass, so it cannot be re-chosen offline the way the thresholds can. The
        # finer bin costs no cluster time, since cost is per-partition and the bin is only
        # a grouping key; what grows is the summary file and the passes that read it.
        Int bin_size = 10000

        String? contigs
        String? intervals
        String? target_superpartitions

        String? inject_dropout

        # Detection thresholds. Cheap to revisit: the detector runs against the small summary
        # table, so re-judging an existing scan costs nothing.
        Float? ratio_threshold
        Float? score_threshold
        Float? min_expected
        Float? min_coverage_fraction
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
        # Two rather than the 4 GvsValidateVDS.wdl uses. Primary workers must all be
        # provisioned before the cluster is usable, and Dataproc's floor is 2 datanodes, so
        # asking for 4 doubles the up-front capacity ask for no benefit -- a us-central1-b
        # create failed having provisioned only 1 of 4. Secondary workers are added by the
        # autoscaling policy afterwards and are preemptible, so peak parallelism is
        # unaffected.
        Int num_primary_workers = 2
        Int max_secondary_workers = 300
        String cluster_prefix = "vds-completeness"
        # Compute Engine stockouts are per-zone, and Dataproc's default placement picks one
        # zone and fails rather than moving on. "auto" tries every zone in the region.
        String cluster_zones = "auto"
        # Kept at 1, matching every other GVS Hail workflow. In principle this screen is a
        # single streaming pass with nothing to spill, and dropping local SSD widens the
        # zones able to serve a request -- but that is now handled by cluster_zones retrying
        # a stockout elsewhere, so there is no reason to also diverge from the configuration
        # that is known to work. Worth revisiting for cost once the pipeline is proven.
        Int num_local_ssds = 1
        String worker_machine_type = "n1-highmem-8"
        # Also left matching the other GVS workflows. The driver only collects a summary on
        # the order of 25,000 x 134 numbers per contig, so this is oversized on paper, but a smaller
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
            help: "Either scan or full-depth. See the header comment for the usual sequence."
        }
        vds_path: {
            help: "VDS to read. Must span more than one GVS superpartition -- see the header comment on which VDSes this screen applies to."
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
        inject_dropout: {
            help: "DIAGNOSTIC ONLY, as contig:start-end:superpartition, e.g. chr20:1000000-1100000:83. Removes that superpartition's data over that window before summarizing, so a known dropout can be watched through the scan and the detector end to end. One window per run. A run using this does not describe the VDS it reads, and its shards are marked so a later clean run refuses to resume from them."
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
        min_coverage_fraction: {
            help: "References mode only: skip bins whose baseline covers less than this fraction of the bin per sample. Defaults to 0.05. Dead sequence -- centromeres, satellite arrays, assembly gaps -- otherwise gets judged on the ratio between two near-zero numbers, which is where every false positive in the Foxtrot r2 reference scan came from. Bins it excludes are listed in sparse_bins_references.tsv. Has no effect in variants mode, where an entry count has no such denominator; use min_expected there."
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

    if (action != "scan" && action != "full-depth") {
        call Utils.TerminateWorkflow as UnrecognizedAction {
            input:
                message = "`action` must be either scan or full-depth; got '" + action + "'.",
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

    # Normalize the prefix, because one carrying a doubled slash cannot be read back. Every
    # path here is built by concatenation, so `gs://b/p/` yields `gs://b/p//summary_x.tsv`.
    # The scan writes through `hl.hadoop_open`, and `org.apache.hadoop.fs.Path` collapses
    # consecutive slashes, so the object lands at the single-slash name; gsutil does not
    # collapse them, because a GCS object name is a literal string in which `//` is simply
    # two characters. The write then succeeds, the scan log looks perfectly healthy, and the
    # `gsutil cp` that follows reports `No URLs matched` for a file sitting right there
    # under a name one character shorter. Worse, `GenerateSampleMap` fails whichever way
    # BigQuery resolves the export URI -- either its `gsutil ls` matches nothing, or it
    # hands the doubled path to a Hail read that normalizes it away -- so the same input
    # can also kill a run before it starts.
    #
    # Collapsing every run of slashes and then repairing the scheme, rather than a
    # lookbehind, because Cromwell's `sub` is specified against POSIX ERE and an engine
    # that rejected `(?<!:)` would fail the workflow outright -- worse than the bug. This
    # also catches an interior `gs://b/a//c`, which a trailing-slash strip does not.
    String collapsed_output_prefix = sub(sub(output_prefix, "/+", "/"), "/$", "")
    String clean_output_prefix = sub(collapsed_output_prefix, "^gs:/", "gs://")

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
                output_prefix = clean_output_prefix,
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
            output_prefix = clean_output_prefix,
            mode = mode,
            sample_map_path = effective_sample_map_path,
            superpartition_size = superpartition_size,
            bin_size = bin_size,
            contigs = contigs,
            intervals = intervals,
            target_superpartitions = target_superpartitions,
            inject_dropout = inject_dropout,
            ratio_threshold = ratio_threshold,
            score_threshold = score_threshold,
            min_expected = min_expected,
            min_coverage_fraction = min_coverage_fraction,
            scale_threshold = scale_threshold,
            baseline_quantile = baseline_quantile,
            bq_project_id = bq_project_id,
            bq_dataset_name = bq_dataset_name,
            reference_schema = effective_reference_schema,
            prefix = cluster_prefix,
            cluster_zones = cluster_zones,
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
        File sparse_bins = ScanVdsForDropouts.sparse_bins
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
        # sample universe matches the one the VDS was built from. This is the only copy of
        # the query; vds_dropout_scan.py documents the format a hand-built map must match.
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

        String? inject_dropout

        Float? ratio_threshold
        Float? score_threshold
        Float? min_expected
        Float? min_coverage_fraction
        Float? scale_threshold
        Float? baseline_quantile

        String? bq_project_id
        String? bq_dataset_name
        String reference_schema

        String prefix
        String cluster_zones
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

        # Both output files must exist because Cromwell resolves task outputs regardless of
        # which branch ran, so each states what it is and why it is empty. These are
        # defaults, overwritten by the detect step below when action is scan. Each file's
        # text has to name that file rather than the other, and has to read correctly for
        # both actions -- the scan case is the only one in which anyone actually reads
        # these, where a placeholder means the task died before detect ran.
        #
        # Appended rather than written as one multi-line string: Cromwell dedents a command
        # block by its common leading whitespace, so a continuation line starting at column
        # zero drops that common prefix to nothing. Nothing is then stripped, the
        # PYTHON_HEREDOC terminator below keeps its indentation, and bash never finds it --
        # `unexpected end of file`, reported at the bottom of the script and pointing
        # nowhere near the string that caused it.
        if [[ "~{action}" == "scan" ]]
        then
            placeholder="the scan did not get as far as the detect step, so this file"
            placeholder="${placeholder} was never written. The task failed earlier -- see"
            placeholder="${placeholder} scan_~{action}_~{mode}.log and the task's stderr."
            placeholder="${placeholder} This file's presence is the symptom, not the cause."
        else
            placeholder="action '~{action}' does not produce this file; only the scan action does."
        fi
        echo "# No findings report: ${placeholder}" > report.tsv
        echo "-- No adjudication SQL: ${placeholder}" > adjudicate.sql
        echo "# No excluded-bin list: ${placeholder}" > sparse_bins.tsv

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
            "temp-path": sys.argv[1],
        }

        # Per-action arguments only, so no command line carries flags another action owns.
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
            ("inject-dropout", "~{default='' inject_dropout}"),
        ]:
            if value:
                arguments[key] = value

        print(json.dumps(arguments, indent=2))
        PYTHON_HEREDOC

        cat script-arguments.json

        if [[ -n "~{default='' inject_dropout}" ]]
        then
            echo "*** INJECTED DROPOUT ~{default='' inject_dropout} -- this run is a diagnostic" >&2
            echo "*** and its outputs do NOT describe ~{vds_path} ***" >&2
        fi

        # Upload the log however this task exits, not just when it succeeds. Copying it at
        # the end of the script instead would make the one artifact worth having after an
        # hours-long failure the one artifact a failure guarantees you do not get: errexit
        # aborts at the first failing step, and every step that can fail comes before the
        # upload, leaving the log in Cromwell's stdout only. `|| true` so a failed upload
        # cannot overwrite the exit code that explains the failure.
        trap 'if [[ -f scan.log ]]
              then
                  gsutil cp scan.log "~{output_prefix}/scan_~{action}_~{mode}.log" || true
              fi' EXIT

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
            # Check both before copying either. A scan that resumed with every contig
            # already checkpointed complete does not rewrite the superpartition table --
            # vds_dropout_scan.py leaves the original in place, on the grounds that it is
            # still correct -- so this path is reachable with a perfectly healthy scan log,
            # where a bare `gsutil cp` failure would name neither the file nor the reason.
            missing=()
            for required_object in "~{summary_path}" "~{superpartitions_path}"
            do
                gsutil -q stat "${required_object}" || missing+=("${required_object}")
            done
            if (( ${#missing[@]} ))
            then
                echo "ERROR: the scan finished but detect cannot run, because these objects" >&2
                echo "are not present:" >&2
                printf '  %s\n' "${missing[@]}" >&2
                echo "" >&2
                echo "The usual cause is a scan that resumed with every contig already" >&2
                echo "complete: it leaves the superpartition table from the original run in" >&2
                echo "place rather than rewriting it, so if that run never wrote one -- or" >&2
                echo "wrote it under a different output_prefix -- nothing here will." >&2
                echo "" >&2
                echo "No need to re-scan. Judging is pure Python and needs no cluster:" >&2
                echo "  python3 vds_dropout_detect.py --summary <summary.tsv> \\" >&2
                echo "    --superpartitions <superpartitions.tsv> --mode ~{mode}" >&2
                echo "Or re-run the scan under a fresh output_prefix to rebuild both." >&2
                exit 1
            fi

            gsutil cp "~{summary_path}" ./summary.tsv
            gsutil cp "~{superpartitions_path}" ./superpartitions.tsv

            detect_args=(
                --summary ./summary.tsv
                --superpartitions ./superpartitions.tsv
                --mode ~{mode}
                --report-path ./report.tsv
                --sparse-bins-path ./sparse_bins.tsv
            )
            ~{'detect_args+=(--ratio-threshold ' + ratio_threshold + ')'}
            ~{'detect_args+=(--score-threshold ' + score_threshold + ')'}
            ~{'detect_args+=(--min-expected ' + min_expected + ')'}
            ~{'detect_args+=(--min-coverage-fraction ' + min_coverage_fraction + ')'}
            ~{'detect_args+=(--scale-threshold ' + scale_threshold + ')'}
            ~{'detect_args+=(--baseline-quantile ' + baseline_quantile + ')'}

            if [[ -n "~{default='' bq_project_id}" && -n "~{default='' bq_dataset_name}" ]]
            then
                detect_args+=(--sql-path ./adjudicate.sql)
                detect_args+=(--project-id "~{default='' bq_project_id}")
                detect_args+=(--dataset-name "~{default='' bq_dataset_name}")
                detect_args+=(--reference-schema ~{reference_schema})
            else
                # Adjudication is half the method: the screen produces candidates and
                # BigQuery decides whether each is real. Skipping it silently would leave
                # unproven findings looking finished, so say so in the file and in the log.
                cat > adjudicate.sql <<'NO_SQL_GENERATED'
        -- No adjudication SQL was generated, because bq_project_id and bq_dataset_name were not
        -- supplied to this workflow.
        --
        -- The screen only identifies candidates; BigQuery is what establishes whether the data is
        -- genuinely missing rather than merely anomalous. Findings are unproven without it.
        --
        -- No need to re-run the scan. The summary it wrote is enough, and judging is pure Python:
        --
        --   python3 vds_dropout_detect.py \
        --     --summary <summary.tsv> --superpartitions <superpartitions.tsv> --mode <mode> \
        --     --project-id <project> --dataset-name <dataset> --sql-path adjudicate.sql
        --
        -- Or pass bq_project_id and bq_dataset_name next time; they can be supplied alongside
        -- sample_map_path, in which case the map is reused and only the SQL is generated.
        NO_SQL_GENERATED
                echo "WARNING: bq_project_id/bq_dataset_name were not supplied, so no" >&2
                echo "         adjudication SQL was generated. The candidates in" >&2
                echo "         report_~{mode}.tsv are unproven until checked against" >&2
                echo "         BigQuery. See adjudicate_~{mode}.sql for how to generate it" >&2
                echo "         without re-running the scan." >&2
            fi

            python3 ~{vds_dropout_detect_script} "${detect_args[@]}" | tee -a scan.log

            gsutil cp ./report.tsv "~{output_prefix}/report_~{mode}.tsv"
            # Unconditional: the placeholder written above guarantees the path exists, and
            # whichever version of the file is there says something true about itself.
            gsutil cp ./adjudicate.sql "~{output_prefix}/adjudicate_~{mode}.sql"
            # In references mode the coverage floor drops bins too sparsely covered to
            # judge, and a clean report is only as strong as the list of places the screen
            # declined to look. Empty but for its header in variants mode, which has no
            # such floor.
            gsutil cp ./sparse_bins.tsv "~{output_prefix}/sparse_bins_~{mode}.tsv"
        fi

        # scan.log is uploaded by the EXIT trap installed above, so there is no copy here.
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
        File sparse_bins = "sparse_bins.tsv"
    }
}
