version 1.0

task PreprocessIntervals {
    input {
      File? intervals
      File? blacklist_intervals
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      Int? padding
      Int? bin_length
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # Determine output filename
    String filename = select_first([intervals, "wgs"])
    String base_filename = basename(filename, ".interval_list")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" PreprocessIntervals \
            ~{"-L " + intervals} \
            ~{"-XL " + blacklist_intervals} \
            --reference ~{ref_fasta} \
            --padding ~{default="250" padding} \
            --bin-length ~{default="1000" bin_length} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.preprocessed.interval_list
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File preprocessed_intervals = "~{base_filename}.preprocessed.interval_list"
    }
}

task AnnotateIntervals {
    input {
      File intervals
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      File? mappability_track_bed
      File? mappability_track_bed_idx
      File? segmental_duplication_track_bed
      File? segmental_duplication_track_bed_idx
      Int? feature_query_lookahead
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500
    
    # Determine output filename
    String base_filename = basename(intervals, ".interval_list")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" AnnotateIntervals \
            -L ~{intervals} \
            --reference ~{ref_fasta} \
            ~{"--mappability-track " + mappability_track_bed} \
            ~{"--segmental-duplication-track " + segmental_duplication_track_bed} \
            --feature-query-lookahead ~{default=1000000 feature_query_lookahead} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.annotated.tsv
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File annotated_intervals = "~{base_filename}.annotated.tsv"
    }
}

task FilterIntervals {
    input {
      File intervals
      File? blacklist_intervals
      File? annotated_intervals
      Array[File]? read_count_files
      Float? minimum_gc_content
      Float? maximum_gc_content
      Float? minimum_mappability
      Float? maximum_mappability
      Float? minimum_segmental_duplication_content
      Float? maximum_segmental_duplication_content
      Int? low_count_filter_count_threshold
      Float? low_count_filter_percentage_of_samples
      Float? extreme_count_filter_minimum_percentile
      Float? extreme_count_filter_maximum_percentile
      Float? extreme_count_filter_percentage_of_samples
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500
    
    # Determine output filename
    String base_filename = basename(intervals, ".interval_list")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" FilterIntervals \
            -L ~{intervals} \
            ~{"-XL " + blacklist_intervals} \
            ~{"--annotated-intervals " + annotated_intervals} \
            ~{if defined(read_count_files) then "--input " else ""} ~{sep=" --input " read_count_files} \
            --minimum-gc-content ~{default="0.1" minimum_gc_content} \
            --maximum-gc-content ~{default="0.9" maximum_gc_content} \
            --minimum-mappability ~{default="0.9" minimum_mappability} \
            --maximum-mappability ~{default="1.0" maximum_mappability} \
            --minimum-segmental-duplication-content ~{default="0.0" minimum_segmental_duplication_content} \
            --maximum-segmental-duplication-content ~{default="0.5" maximum_segmental_duplication_content} \
            --low-count-filter-count-threshold ~{default="5" low_count_filter_count_threshold} \
            --low-count-filter-percentage-of-samples ~{default="90.0" low_count_filter_percentage_of_samples} \
            --extreme-count-filter-minimum-percentile ~{default="1.0" extreme_count_filter_minimum_percentile} \
            --extreme-count-filter-maximum-percentile ~{default="99.0" extreme_count_filter_maximum_percentile} \
            --extreme-count-filter-percentage-of-samples ~{default="90.0" extreme_count_filter_percentage_of_samples} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.filtered.interval_list
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File filtered_intervals = "~{base_filename}.filtered.interval_list"
    }
}

task CollectCounts {
    input {
      File intervals
      File bam
      File bam_idx
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      Array[String]? disabled_read_filters
      Boolean? enable_indexing
      String? format
      File? gatk4_jar_override
      String? gcs_project_for_requester_pays

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    parameter_meta {
      bam: {
        localization_optional: true
      }
      bam_idx: {
        localization_optional: true
      }
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    Boolean enable_indexing_ = select_first([enable_indexing, false])

    Array[String] disabled_read_filters_arr = if defined(disabled_read_filters) then prefix("--disable-read-filter ", select_first([disabled_read_filters])) else []

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")
    String format_ = select_first([format, "HDF5"])
    String hdf5_or_tsv_or_null_format =
        if format_ == "HDF5" then "HDF5" else
        (if format_ == "TSV" then "TSV" else
        (if format_ == "TSV_GZ" then "TSV" else "null")) # until we can write TSV_GZ in CollectReadCounts, we write TSV and use bgzip
    String counts_filename_extension =
        if format_ == "HDF5" then "counts.hdf5" else
        (if format_ == "TSV" then "counts.tsv" else
        (if format_ == "TSV_GZ" then "counts.tsv.gz" else "null"))
    String counts_index_filename_extension =
        if format_ == "HDF5" then "null" else
        (if format_ == "TSV" then "counts.tsv.idx" else
        (if format_ == "TSV_GZ" then "counts.tsv.gz.tbi" else "null"))
    Boolean do_block_compression =
        if format_ == "HDF5" then false else
        (if format_ == "TSV" then false else
        (if format_ == "TSV_GZ" then true else false))
    String counts_filename = "~{base_filename}.~{counts_filename_extension}"
    String counts_filename_for_collect_read_counts = basename(counts_filename, ".gz")

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        case ~{format_} in
            HDF5 | TSV | TSV_GZ)
                ;;
            *)
                echo "ERROR: Unknown format specified. Format must be one of HDF5, TSV, or TSV_GZ."
                exit 1
                ;;
        esac

        if [ ~{format_} = "HDF5" ] && [ ~{enable_indexing_} = "true" ]; then
            echo "ERROR: Incompatible WDL parameters. Cannot have format = HDF5 and enable_indexing = true."
            exit 1
        fi

        if [ ~{hdf5_or_tsv_or_null_format} = "null" ]; then
            echo "ERROR: Should never reach here."
            exit 1
        fi

        gatk --java-options "-Xmx~{command_mem_mb}m" CollectReadCounts \
            -L ~{intervals} \
            --input ~{bam} \
            --read-index ~{bam_idx} \
            --reference ~{ref_fasta} \
            --format ~{default="HDF5" hdf5_or_tsv_or_null_format} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{counts_filename_for_collect_read_counts} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays} \
            ~{sep=' ' disabled_read_filters_arr}

        if [ ~{do_block_compression} = "true" ]; then
            bgzip ~{counts_filename_for_collect_read_counts}
        fi

        if [ ~{enable_indexing_} = "true" ]; then
            gatk --java-options "-Xmx~{command_mem_mb}m" IndexFeatureFile \
                -I ~{counts_filename}
        fi
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        String entity_id = base_filename
        File counts = counts_filename
    }
}

task CollectAllelicCounts {
    input {
      File common_sites
      File bam
      File bam_idx
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      Int? minimum_base_quality
      File? gatk4_jar_override
      String? gcs_project_for_requester_pays

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    parameter_meta {
      bam: {
        localization_optional: true
      }
      bam_idx: {
        localization_optional: true
      }
    }

    Int machine_mem_mb = select_first([mem_gb, 13]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")

    String allelic_counts_filename = "~{base_filename}.allelicCounts.tsv"

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" CollectAllelicCounts \
            -L ~{common_sites} \
            --input ~{bam} \
            --read-index ~{bam_idx} \
            --reference ~{ref_fasta} \
            --minimum-base-quality ~{default="20" minimum_base_quality} \
            --output ~{allelic_counts_filename} \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        String entity_id = base_filename
        File allelic_counts = allelic_counts_filename
    }
}

# Users should consult the IntervalListTools documentation and/or manually inspect the results of this task
# to ensure that the number of intervals in each shard is as desired, as the logic IntervalListTools uses
# for dividing intervals can yield shards that are unexpectedly larger than num_intervals_per_scatter.
# Depending on their use case, users may want to modify this task to instead use the SCATTER_COUNT option of
# IntervalListTools, which allows the number of shards to be directly specified.
task ScatterIntervals {
    input {
      File interval_list
      Int num_intervals_per_scatter
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out";
    String output_dir_ = select_first([output_dir, "out"])

    String base_filename = basename(interval_list, ".interval_list")

    command <<<
        set -eu
        # IntervalListTools will fail if the output directory does not exist, so we create it
        mkdir ~{output_dir_}
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        # IntervalListTools behaves differently when scattering to a single or multiple shards, so we do some handling in bash

        # IntervalListTools tries to equally divide intervals across shards to give at least INTERVAL_COUNT in each and
        # puts remainder intervals in the last shard, so integer division gives the number of shards
        # (unless NUM_INTERVALS < num_intervals_per_scatter and NUM_SCATTERS = 0, in which case we still want a single shard)
        NUM_INTERVALS=$(grep -v '@' ~{interval_list} | wc -l)
        NUM_SCATTERS=$(echo $((NUM_INTERVALS / ~{num_intervals_per_scatter})))

        if [ $NUM_SCATTERS -le 1 ]; then
            # if only a single shard is required, then we can just rename the original interval list
            >&2 echo "Not running IntervalListTools because only a single shard is required. Copying original interval list..."
            cp ~{interval_list} ~{output_dir_}/~{base_filename}.scattered.0001.interval_list
        else
            gatk --java-options "-Xmx~{command_mem_mb}m" IntervalListTools \
                --INPUT ~{interval_list} \
                --SUBDIVISION_MODE INTERVAL_COUNT \
                --SCATTER_CONTENT ~{num_intervals_per_scatter} \
                --OUTPUT ~{output_dir_}

            # output files are named output_dir_/temp_0001_of_N/scattered.interval_list, etc. (N = number of scatters);
            # we rename them as output_dir_/base_filename.scattered.0001.interval_list, etc.
            ls -v ~{output_dir_}/*/scattered.interval_list | \
                cat -n | \
                while read n filename; do mv $filename ~{output_dir_}/~{base_filename}.scattered.$(printf "%04d" $n).interval_list; done
            rm -rf ~{output_dir_}/temp_*_of_*
        fi
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] scattered_interval_lists = glob("~{output_dir_}/~{base_filename}.scattered.*.interval_list")
    }
}

task PostprocessGermlineCNVCalls {
    input {
      String entity_id
      Array[File] gcnv_calls_tars
      Array[File] gcnv_model_tars
      Array[File] calling_configs
      Array[File] denoising_configs
      Array[File] gcnvkernel_version
      Array[File] sharded_interval_lists
      File contig_ploidy_calls_tar
      Array[String]? allosomal_contigs
      Int ref_copy_number_autosomal_contigs
      Int sample_index
      Int maximum_number_events
      Int maximum_number_pass_events
      File? intervals_vcf
      File? intervals_vcf_index
      File? clustered_vcf
      File? clustered_vcf_index
      File? reference_fasta
      File? reference_fasta_fai
      File? reference_dict
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    String genotyped_intervals_vcf_filename = "genotyped-intervals-~{entity_id}.vcf.gz"
    String genotyped_segments_vcf_filename = "genotyped-segments-~{entity_id}.vcf.gz"
    String denoised_copy_ratios_filename = "denoised_copy_ratios-~{entity_id}.tsv"
    String qc_status_filename = "~{entity_id}.qcStatus.txt"

    Array[String] allosomal_contigs_args = if defined(allosomal_contigs) then prefix("--allosomal-contig ", select_first([allosomal_contigs])) else []

    command <<<
        set -eu
        ~{"export GATK_LOCAL_JAR=" + gatk4_jar_override}

        sharded_interval_lists_array=(~{sep=" " sharded_interval_lists})

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        # also copy over shard config and interval files
        gcnv_calls_tar_array=(~{sep=" " gcnv_calls_tars})
        calling_configs_array=(~{sep=" " calling_configs})
        denoising_configs_array=(~{sep=" " denoising_configs})
        gcnvkernel_version_array=(~{sep=" " gcnvkernel_version})
        sharded_interval_lists_array=(~{sep=" " sharded_interval_lists})
        calls_args=""
        for index in ${!gcnv_calls_tar_array[@]}; do
            gcnv_calls_tar=${gcnv_calls_tar_array[$index]}
            mkdir -p CALLS_$index/SAMPLE_~{sample_index}
            tar xzf $gcnv_calls_tar -C CALLS_$index/SAMPLE_~{sample_index}
            cp ${calling_configs_array[$index]} CALLS_$index/
            cp ${denoising_configs_array[$index]} CALLS_$index/
            cp ${gcnvkernel_version_array[$index]} CALLS_$index/
            cp ${sharded_interval_lists_array[$index]} CALLS_$index/
            calls_args="$calls_args --calls-shard-path CALLS_$index"
        done

        # untar models to MODEL_0, MODEL_1, etc directories and build the command line
        gcnv_model_tar_array=(~{sep=" " gcnv_model_tars})
        model_args=""
        for index in ${!gcnv_model_tar_array[@]}; do
            gcnv_model_tar=${gcnv_model_tar_array[$index]}
            mkdir MODEL_$index
            tar xzf $gcnv_model_tar -C MODEL_$index
            model_args="$model_args --model-shard-path MODEL_$index"
        done

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        gatk --java-options "-Xmx~{command_mem_mb}m" PostprocessGermlineCNVCalls \
            $calls_args \
            $model_args \
            ~{sep=" " allosomal_contigs_args} \
            --autosomal-ref-copy-number ~{ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls contig-ploidy-calls \
            --sample-index ~{sample_index} \
            --output-genotyped-intervals ~{genotyped_intervals_vcf_filename} \
            --output-genotyped-segments ~{genotyped_segments_vcf_filename} \
            --output-denoised-copy-ratios ~{denoised_copy_ratios_filename} \
            ~{"--input-intervals-vcf " + intervals_vcf} \
            ~{"--clustered-breakpoints " + clustered_vcf} \
            ~{"-R " + reference_fasta}

        #use wc instead of grep -c so zero count isn't non-zero exit
        #use grep -P to recognize tab character
        NUM_SEGMENTS=$(zgrep '^[^#]' ~{genotyped_segments_vcf_filename} | grep -v '0/0' | grep -v -P '\t0:1:' | grep '' | wc -l)
        NUM_PASS_SEGMENTS=$(zgrep '^[^#]' ~{genotyped_segments_vcf_filename} | grep -v '0/0' | grep -v -P '\t0:1:' | grep 'PASS' | wc -l)
        if [ $NUM_SEGMENTS -lt ~{maximum_number_events} ]; then
            if [ $NUM_PASS_SEGMENTS -lt ~{maximum_number_pass_events} ]; then
              echo "PASS" >> ~{qc_status_filename}
            else
              echo "EXCESSIVE_NUMBER_OF_PASS_EVENTS" >> ~{qc_status_filename}
            fi
        else
            echo "EXCESSIVE_NUMBER_OF_EVENTS" >> ~{qc_status_filename}
        fi

        rm -rf CALLS_*
        rm -rf MODEL_*
        rm -rf contig-ploidy-calls
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File genotyped_intervals_vcf = genotyped_intervals_vcf_filename
        File genotyped_intervals_vcf_index = genotyped_intervals_vcf_filename + ".tbi"
        File genotyped_segments_vcf = genotyped_segments_vcf_filename
        File genotyped_segments_vcf_index = genotyped_segments_vcf_filename + ".tbi"
        File denoised_copy_ratios = denoised_copy_ratios_filename
        File qc_status_file = qc_status_filename
        String qc_status_string = read_string(qc_status_filename)
    }
}

task CollectModelQualityMetrics {
    input {
      Array[File] gcnv_model_tars

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 1]) * 1000

    command <<<
        sed -e 
        qc_status="PASS"

        gcnv_model_tar_array=(~{sep=" " gcnv_model_tars})
        for index in ${!gcnv_model_tar_array[@]}; do
            gcnv_model_tar=${gcnv_model_tar_array[$index]}
            mkdir MODEL_$index
            tar xzf $gcnv_model_tar -C MODEL_$index
            ard_file=MODEL_$index/mu_ard_u_log__.tsv

            #check whether all values for ARD components are negative
            NUM_POSITIVE_VALUES=$(awk '{ if (index($0, "@") == 0) {if ($1 > 0.0) {print $1} }}' MODEL_$index/mu_ard_u_log__.tsv | wc -l)
            if [ $NUM_POSITIVE_VALUES -eq 0 ]; then
                qc_status="ALL_PRINCIPAL_COMPONENTS_USED"
                break
            fi
        done
        echo $qc_status >> qcStatus.txt
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File qc_status_file = "qcStatus.txt"
        String qc_status_string = read_string("qcStatus.txt")
    }
}

task SplitInputArray {
    input {
      Array[String] input_array
      Int num_inputs_in_scatter_block
      String gatk_docker

      Int machine_mem_mb = 4000
      Int disk_space_gb = 20
      Int cpu = 1
      Int? preemptible_attempts
      Boolean use_ssd = false
    }

    File input_array_file = write_lines(input_array)

    # This tasks takes as input an array of strings and number of columns (num_inputs_in_scatter_block)
    # and outputs a 2-dimensional reshaped array with same contents and with width equal to num_inputs_in_scatter_block
    # (with last row potentially having a smaller length than others)
    command <<<
        python <<CODE
        import math
        with open("~{input_array_file}", "r") as input_array_file:
            input_array = input_array_file.read().splitlines()
        num = ~{num_inputs_in_scatter_block}
        values_to_write = [input_array[num*i:num*i+min(num, len(input_array)-num*i)] for i in range(int(math.ceil(len(input_array)/num)))]
        with open('input_array_split.tsv', 'w') as outfile:
            for i in range(len(values_to_write)):
                current_sub_array = values_to_write[i]
                for j in range(len(current_sub_array)):
                    outfile.write(current_sub_array[j] + "\t")
                outfile.write("\n")
        CODE
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + disk_space_gb + if use_ssd then " SSD" else " HDD"
        cpu: cpu
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[Array[String]] split_array = read_tsv("input_array_split.tsv")
    }
}

task ScatterPloidyCallsBySample {
    input {
      File contig_ploidy_calls_tar
      Array[String] samples

      # Runtime parameters
      String docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int num_samples = length(samples)
    String out_dir = "calls_renamed"

    command <<<
      set -eu

      # Extract ploidy calls
      mkdir calls
      tar xzf ~{contig_ploidy_calls_tar} -C calls/

      # Archive call files by sample, renaming so they will be glob'd in order
      sample_ids=(~{sep=" " samples})
      num_samples=~{num_samples}
      num_digits=${#num_samples}
      for (( i=0; i<~{num_samples}; i++ ))
      do
        sample_id=${sample_ids[$i]}
        padded_sample_index=$(printf "%0${num_digits}d" $i)
        tar -czf sample_${padded_sample_index}.${sample_id}.contig_ploidy_calls.tar.gz -C calls/SAMPLE_${i} .
      done
    >>>

    runtime {
        docker: docker
        memory: select_first([mem_gb, 2]) + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, 10]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] sample_contig_ploidy_calls_tar = glob("sample_*.contig_ploidy_calls.tar.gz")
    }
}