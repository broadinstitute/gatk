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
    Array[String] disabled_read_filters_arr = if(defined(disabled_read_filters))
      then
        prefix(
          "--disable-read-filter ",
          select_first([disabled_read_filters])
        )
      else
        []

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
            --reference ~{ref_fasta} \
            --format ~{default="HDF5" hdf5_or_tsv_or_null_format} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{counts_filename_for_collect_read_counts} \
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
            --reference ~{ref_fasta} \
            --minimum-base-quality ~{default="20" minimum_base_quality} \
            --output ~{allelic_counts_filename}
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
        File bundled_gcnv_outputs
        String entity_id
        File contig_ploidy_calls_tar
        Array[String]? allosomal_contigs
        Int ref_copy_number_autosomal_contigs
        Int sample_index
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

    Float bundled_gcnv_outputs_size = size(bundled_gcnv_outputs, "GiB")
    Float disk_overhead = 20.0
    Float tar_disk_factor= 5.0
    Int vm_disk_size = ceil(tar_disk_factor * bundled_gcnv_outputs_size + disk_overhead)

    String genotyped_intervals_vcf_filename = "genotyped-intervals-~{entity_id}.vcf.gz"
    String genotyped_segments_vcf_filename = "genotyped-segments-~{entity_id}.vcf.gz"
    String denoised_copy_ratios_filename = "denoised_copy_ratios-~{entity_id}.tsv"

    Array[String] allosomal_contigs_args = if defined(allosomal_contigs) then prefix("--allosomal-contig ", select_first([allosomal_contigs])) else []

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        # also copy over shard config and interval files
        tar xzf ~{bundled_gcnv_outputs}
        rm ~{bundled_gcnv_outputs}
        number_of_shards=`find . -name 'CALLS_*' | wc -l`

        touch calls_and_model_args.txt
        for i in $(seq 0 `expr $number_of_shards - 1`); do
            echo "--calls-shard-path CALLS_$i" >> calls_and_model_args.txt
            echo "--model-shard-path MODEL_$i" >> calls_and_model_args.txt
        done

        mkdir -p extracted-contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C extracted-contig-ploidy-calls
        rm ~{contig_ploidy_calls_tar}

        gatk --java-options "-Xmx~{command_mem_mb}m" PostprocessGermlineCNVCalls \
             --arguments_file calls_and_model_args.txt \
            ~{sep=" " allosomal_contigs_args} \
            --autosomal-ref-copy-number ~{ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls extracted-contig-ploidy-calls \
            --sample-index ~{sample_index} \
            --output-genotyped-intervals ~{genotyped_intervals_vcf_filename} \
            --output-genotyped-segments ~{genotyped_segments_vcf_filename} \
            --output-denoised-copy-ratios ~{denoised_copy_ratios_filename}

        rm -rf CALLS_*
        rm -rf MODEL_*
        rm -rf extracted-contig-ploidy-calls
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, vm_disk_size]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
        maxRetries: 1
    }

    output {
        File genotyped_intervals_vcf = genotyped_intervals_vcf_filename
        File genotyped_segments_vcf = genotyped_segments_vcf_filename
        File denoised_copy_ratios = denoised_copy_ratios_filename
    }
}

task CollectSampleQualityMetrics {
    input {
      File genotyped_segments_vcf
      String entity_id
      Int maximum_number_events

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
        set -eu
        NUM_SEGMENTS=$(gunzip -c ~{genotyped_segments_vcf} | grep -v '#' | wc -l)
        if [ $NUM_SEGMENTS -lt ~{maximum_number_events} ]; then
            echo "PASS" >> ~{entity_id}.qcStatus.txt
        else 
            echo "EXCESSIVE_NUMBER_OF_EVENTS" >> ~{entity_id}.qcStatus.txt
        fi
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 20]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File qc_status_file = "~{entity_id}.qcStatus.txt"
        String qc_status_string = read_string("~{entity_id}.qcStatus.txt")
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

task BundleCallerOutputs {
    input {
      Array[File] calls_tars
      Array[File] model_tars
      Array[File] calling_configs
      Array[File] denoising_configs
      Array[File] gcnvkernel_version
      Array[File] sharded_interval_lists

      # Runtime parameters
      String docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    command <<<
        set -euo pipefail
        mkdir -p out

        calls_files_tar_list=~{write_lines(calls_tars)}
        model_files_tar_list=~{write_lines(model_tars)}

        calling_configs_list=~{write_lines(calling_configs)}
        denoising_configs_list=~{write_lines(denoising_configs)}
        gcnvkernel_version_list=~{write_lines(gcnvkernel_version)}
        sharded_interval_lists_list=~{write_lines(sharded_interval_lists)}

        cat $calls_files_tar_list | sort -V > calls_files_tar_list.sorted
        cat $model_files_tar_list | sort -V > model_files_tar_list.sorted

        cat $calling_configs_list | sort -V > calling_configs_list.sorted
        cat $denoising_configs_list | sort -V > denoising_configs_list.sorted
        cat $gcnvkernel_version_list | sort -V > gcnvkernel_version_list.sorted
        cat $sharded_interval_lists_list | sort -V > sharded_interval_lists_list.sorted

        paste calls_files_tar_list.sorted model_files_tar_list.sorted calling_configs_list.sorted denoising_configs_list.sorted gcnvkernel_version_list.sorted sharded_interval_lists_list.sorted |\
              awk '{print (NR-1)"\t"$0}' > file_sets.sorted
        OIFS=$IFS
        IFS=$'\t'
        while read index calls_tar model_tar call_config denoise version intervals; do
            mkdir -p out/CALLS_$index
            mkdir -p out/MODEL_$index
            tar xzf $calls_tar -C out/CALLS_$index
            tar xzf $model_tar -C out/MODEL_$index
            cp $call_config out/CALLS_$index
            cp $denoise out/CALLS_$index
            cp $version out/CALLS_$index
            cp $intervals out/CALLS_$index
            rm $calls_tar $model_tar $call_config $denoise $version $intervals

        done < file_sets.sorted
        IFS=$OIFS

        tar c -C out . | gzip -1 > case-gcnv-postprocessing-invariants.tar.gz
        rm -Rf out
    >>>

    runtime {
        docker: docker
        memory: select_first([mem_gb, 2]) + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File bundle_tar = "case-gcnv-postprocessing-invariants.tar.gz"
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
      for (( i=0; i<~{num_samples}; i++ ))
      do
        sample_id=${sample_ids[$i]}
        sample_no=`printf %04d $i`
        tar -czf sample_${sample_no}.${sample_id}.contig_ploidy_calls.tar.gz -C calls/SAMPLE_${i} .
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