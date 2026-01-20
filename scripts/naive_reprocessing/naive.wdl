version 1.0

workflow ReprocessAndValidate {
  input {
    String rlg="1000000"
    File umil="gs://fc-d26a03c8-5f34-4452-93ff-b3fc13cc2950/naive/k100.Unique.Mappability.ucsc.bed"
    File   input_bam="gs://fc-d26a03c8-5f34-4452-93ff-b3fc13cc2950/naive/NA12878.cram"
    File   input_bam_idx="gs://fc-d26a03c8-5f34-4452-93ff-b3fc13cc2950/naive/NA12878.cram.crai"
		File ref_fasta="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
		File ref_fasta_fai="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    File ref_dict="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    File gatk_jar = "gs://fc-d26a03c8-5f34-4452-93ff-b3fc13cc2950/naive/gatkdev_naive.jar"
    Int java_mem = 4
    String output_prefix="test"
    
    # Parallelization option - if true, splits by chromosome groups for faster processing
    # Recommended: true (enables preemptible VMs for cost savings)
    Boolean parallelize = true
    
    # Two-stage optimization - if true, uses memory-efficient two-stage processing:
    # Stage 1: Generate uncertain CRAM + .fcl file
    # Stage 2: Filter input CRAM using .fcl to create naive CRAM
    # This uses 85% less memory (1.5 GB vs 10.5 GB for 38M read names)
    # Note: Not compatible with parallelize=true
    Boolean two_stage = false
    
    # Parallel processing options (only used when parallelize=true)
    # 8 shards optimized for ~equal runtime and preemptible VM safety (~1.1 hours each)
    Array[Array[String]] chromosome_groups = [
      ["chr1"],           # Largest chromosome
      ["chr2"],           # Second largest
      ["chr3", "chr4"],   # Medium chromosomes
      ["chr5", "chr6"],   # Medium chromosomes
      ["chr7", "chr8", "chr9"],      # Smaller chromosomes
      ["chr10", "chr11", "chr12"],   # Smaller chromosomes
      ["chr13", "chr14", "chr15", "chr16", "chr17", "chr18"],  # Small chromosomes
      ["chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]  # Smallest + sex + mito
    ]
    Array[String] shard_names = [
      "shard1_chr1", 
      "shard2_chr2", 
      "shard3_chr3_4", 
      "shard4_chr5_6",
      "shard5_chr7_9",
      "shard6_chr10_12",
      "shard7_chr13_18",
      "shard8_chr19_M"
    ]
    Array[Int] expected_uncertain_per_shard = [
      10000000,  # chr1
      9000000,   # chr2
      8000000,   # chr3-4
      8000000,   # chr5-6
      7000000,   # chr7-9
      7000000,   # chr10-12
      6000000,   # chr13-18
      5000000    # chr19-M
    ]
  }

  # ============================================================================
  # TWO-STAGE OPTIMIZED PATH (two_stage=true) - DEFAULT
  # Uses memory-efficient two-stage processing with front-coded binary lists
  # ============================================================================
  if (two_stage) {
    call naive_processing_two_stage {
      input:
        rlg = rlg,
        umil = umil,
        input_bam = input_bam,
        input_bam_idx = input_bam_idx,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_dict = ref_dict,
        gatk_jar = gatk_jar,
        output_prefix = output_prefix,
        java_mem = java_mem
    }
  }

  # ============================================================================
  # SEQUENTIAL PROCESSING PATH (parallelize=false, two_stage=false)
  # ============================================================================
  if (!parallelize && !two_stage) {
    call naive_processing {
      input:
        rlg = rlg,
        umil = umil,
        input_bam = input_bam,
        input_bam_idx = input_bam_idx,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_dict = ref_dict,
        gatk_jar = gatk_jar,
        output_prefix = output_prefix,
        java_mem = java_mem
    }
  }

  # ============================================================================
  # PARALLEL PROCESSING PATH (parallelize=true, two_stage=false)
  # ============================================================================
  if (parallelize && !two_stage) {
    # Scatter across chromosome groups
    scatter (idx in range(length(chromosome_groups))) {
      call naive_processing_shard {
        input:
          rlg = rlg,
          umil = umil,
          input_bam = input_bam,
          input_bam_idx = input_bam_idx,
          ref_fasta = ref_fasta,
          ref_fasta_fai = ref_fasta_fai,
          ref_dict = ref_dict,
          gatk_jar = gatk_jar,
          java_mem = java_mem,
          output_prefix = output_prefix,
          shard_name = shard_names[idx],
          chromosomes = chromosome_groups[idx],
          expected_uncertain_reads = expected_uncertain_per_shard[idx]
      }
    }

    # Gather: merge naive CRAMs
    call merge_crams as merge_naive {
      input:
        input_crams = naive_processing_shard.naive_bam,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_dict = ref_dict,
        gatk_jar = gatk_jar,
        java_mem = java_mem,
        output_prefix = output_prefix,
        output_suffix = "naive"
    }

    # Gather: merge uncertain CRAMs
    call merge_crams as merge_uncertain {
      input:
        input_crams = naive_processing_shard.uncertain_bam,
        ref_fasta = ref_fasta,
        ref_fasta_fai = ref_fasta_fai,
        ref_dict = ref_dict,
        gatk_jar = gatk_jar,
        java_mem = java_mem,
        output_prefix = output_prefix,
        output_suffix = "uncertain"
    }
  }

  # ============================================================================
  # SELECT OUTPUTS FROM ANY PATH
  # ============================================================================
  File final_naive_bam = select_first([
    naive_processing_two_stage.naive_bam,
    merge_naive.merged_cram, 
    naive_processing.naive_bam
  ])
  File final_uncertain_bam = select_first([
    naive_processing_two_stage.uncertain_bam,
    merge_uncertain.merged_cram, 
    naive_processing.uncertain_bam
  ])

  output {
    File naive_bam                 = final_naive_bam
    File uncertain_bam             = final_uncertain_bam
    
    # Parallel processing debug outputs (only populated when parallelize=true)
    Array[File]? shard_naive_bams     = naive_processing_shard.naive_bam
    Array[File]? shard_uncertain_bams = naive_processing_shard.uncertain_bam
    
    # Two-stage mode output (only populated when two_stage=true)
    File? uncertain_read_names = naive_processing_two_stage.uncertain_read_names
  }
}

task naive_processing {
  input {
    String rlg
    String umil
    File   input_bam
    File   input_bam_idx
    File   ref_fasta
		File ref_fasta_fai
    File ref_dict
    File gatk_jar
    Int java_mem
    String output_prefix
  }

  command <<<
    set -e
    echo "NAIVE PROCESSING ============"
    #java -jar -Xmx~{java_mem}G ~{gatk_jar} NaiveReprocessReads \
    #  -RLG ~{rlg} \
    #  -I ~{input_bam} \
    #  -O ~{output_prefix}.naive.cram \
    #  -Ou ~{output_prefix}.uncertain.cram \
    #  -umil ~{umil} \
    #  -R ~{ref_fasta}
    java -jar -Xmx~{java_mem}G ~{gatk_jar} SplitReadsByRealignmentDifficulty \
      -RLG ~{rlg} \
      -I ~{input_bam} \
      -O ~{output_prefix}.naive.cram \
      -Ou ~{output_prefix}.uncertain.cram \
      -umil ~{umil} \
      -R ~{ref_fasta}
  >>>

  output {
    File naive_bam     = "~{output_prefix}.naive.cram"
    File uncertain_bam = "~{output_prefix}.uncertain.cram"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 2
    memory: "~{java_mem + 2}G"
    disks: "local-disk 1500 SSD"
  }
}

# ============================================================================
# UNCERTAIN-ONLY PROCESSING (optimized I/O - only writes uncertain reads)
# ============================================================================
# This task writes only uncertain reads during GATK traversal, then uses
# samtools to create the naive CRAM by subtracting uncertain read names
# from the original input. This dramatically reduces I/O time since
# uncertain reads are typically <10% of total.

task naive_processing_uncertain_only {
  input {
    String rlg
    File   umil
    File   input_bam
    File   input_bam_idx
    File   ref_fasta
    File   ref_fasta_fai
    File   ref_dict
    File   gatk_jar
    Int    java_mem
    String output_prefix
    
    # Bloom filter parameters
    Int expected_uncertain_reads = 100000000
    Float bloom_filter_fpr = 0.001
  }

  command <<<
    set -e
    echo "NAIVE PROCESSING (UNCERTAIN-ONLY MODE) ============"
    echo "Step 1: Write only uncertain reads + read names file..."
    
    java -jar -Xmx~{java_mem}G ~{gatk_jar} SplitReadsByRealignmentDifficulty \
      -RLG ~{rlg} \
      -I ~{input_bam} \
      -O /dev/null \
      -Ou ~{output_prefix}.uncertain.cram \
      -umil ~{umil} \
      -R ~{ref_fasta} \
      -UO true \
      -URN ~{output_prefix}.uncertain_read_names.txt \
      -EUR ~{expected_uncertain_reads} \
      -FPR ~{bloom_filter_fpr}
    
    echo "Step 1 complete. Uncertain CRAM and read names file created."
    echo "Uncertain read names count: $(wc -l < ~{output_prefix}.uncertain_read_names.txt)"
    
    echo "Step 2: Create naive CRAM by excluding uncertain read names from input..."
    # Use GATK FilterSamReads to create naive CRAM by excluding uncertain reads
    gatk FilterSamReads \
      -I ~{input_bam} \
      -O ~{output_prefix}.naive.cram \
      -R ~{ref_fasta} \
      --FILTER excludeReadList \
      --READ_LIST_FILE ~{output_prefix}.uncertain_read_names.txt
    
    echo "Step 2 complete. Naive CRAM created."
    ls -lh ~{output_prefix}.*.cram
    echo "Done!"
  >>>

  output {
    File naive_bam     = "~{output_prefix}.naive.cram"
    File uncertain_bam = "~{output_prefix}.uncertain.cram"
    File uncertain_read_names = "~{output_prefix}.uncertain_read_names.txt"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 4
    memory: "~{java_mem + 4}G"
    disks: "local-disk 1500 SSD"
  }
}

# ============================================================================
# TWO-STAGE PROCESSING: Generate uncertain CRAM + .fcl, then filter to naive
# ============================================================================
# This task uses the memory-efficient two-stage approach:
# Stage 1: SplitReadsByRealignmentDifficulty writes uncertain CRAM + .fcl file
# Stage 2: FilterReadsByNameList uses .fcl to filter input -> naive CRAM
# The .fcl format uses ~85% less memory than HashSet (1.5 GB vs 10.5 GB for 38M names)

task naive_processing_two_stage {
  input {
    String rlg
    File   umil
    File   input_bam
    File   input_bam_idx
    File   ref_fasta
    File   ref_fasta_fai
    File   ref_dict
    File   gatk_jar
    Int    java_mem
    String output_prefix
    
    # Bloom filter parameters for Stage 1
    Int expected_uncertain_reads = 100000000
    Float bloom_filter_fpr = 0.001
  }

  command <<<
    set -e
    echo "TWO-STAGE NAIVE PROCESSING ============"
    echo "Single-pass: Generate both naive and uncertain CRAMs + .fcl file..."
    echo ""
    
    java -jar -Xmx~{java_mem}G ~{gatk_jar} SplitReadsByRealignmentDifficulty \
      -RLG ~{rlg} \
      -I ~{input_bam} \
      -O ~{output_prefix}.naive.cram \
      -Ou ~{output_prefix}.uncertain.cram \
      -umil ~{umil} \
      -R ~{ref_fasta} \
      -UO true \
      -URN ~{output_prefix}.uncertain_read_names.fcl \
      -EUR ~{expected_uncertain_reads} \
      -FPR ~{bloom_filter_fpr}
    
    echo ""
    echo "============================================"
    echo "PROCESSING COMPLETE"
    echo "Output files:"
    ls -lh ~{output_prefix}.*.cram ~{output_prefix}.*.fcl
    echo "============================================"
  >>>

  output {
    File naive_bam     = "~{output_prefix}.naive.cram"
    File uncertain_bam = "~{output_prefix}.uncertain.cram"
    File uncertain_read_names = "~{output_prefix}.uncertain_read_names.fcl"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 2
    memory: "~{java_mem + 2}G"
    disks: "local-disk 1500 SSD"
  }
}

# ============================================================================
# PARALLEL PROCESSING TASKS
# ============================================================================

# Process a single chromosome group shard
task naive_processing_shard {
  input {
    String rlg
    File   umil  # Full BED file - we'll filter it
    File   input_bam
    File   input_bam_idx
    File   ref_fasta
    File   ref_fasta_fai
    File   ref_dict
    File   gatk_jar
    Int    java_mem
    String output_prefix
    String shard_name
    Array[String] chromosomes  # Chromosomes to process in this shard
    
    # Bloom filter parameters
    Int expected_uncertain_reads = 50000000  # Lower per-shard estimate
    Float bloom_filter_fpr = 0.001
  }

  command <<<
    set -e
    echo "NAIVE PROCESSING SHARD: ~{shard_name} ============"
    echo "Processing chromosomes: ~{sep=' ' chromosomes}"
    
    # Create intervals argument for GATK
    INTERVALS=""
    for chr in ~{sep=' ' chromosomes}; do
      INTERVALS="$INTERVALS -L $chr"
    done
    echo "Intervals: $INTERVALS"
    
    # Filter BED file to only include regions in our chromosomes
    # This speeds up overlap detection significantly and reduces memory usage
    echo "Filtering BED file for shard chromosomes..."
    echo "Original BED file size: $(wc -l < ~{umil}) lines"
    
    CHROM_PATTERN=$(echo "~{sep='|' chromosomes}" | sed 's/|/\\|/g')
    grep -E "^($CHROM_PATTERN)\t" ~{umil} > filtered_umil.bed || true
    
    # If no BED regions for these chromosomes, create empty file
    if [ ! -s filtered_umil.bed ]; then
      echo "Warning: No BED regions found for chromosomes ~{sep=' ' chromosomes}"
      touch filtered_umil.bed
    fi
    
    FILTERED_LINES=$(wc -l < filtered_umil.bed)
    echo "Filtered BED lines: $FILTERED_LINES"
    if [ "$(wc -l < ~{umil})" -gt 0 ]; then
      REDUCTION=$(awk "BEGIN {printf \"%.1f\", (1 - $FILTERED_LINES / $(wc -l < ~{umil})) * 100}")
      echo "BED file reduced by ${REDUCTION}% for this shard"
    fi
    
    # Run SplitReadsByRealignmentDifficulty on this shard
    java -jar -Xmx4G ~{gatk_jar} SplitReadsByRealignmentDifficulty \
      -RLG ~{rlg} \
      -I ~{input_bam} \
      $INTERVALS \
      -O ~{output_prefix}.~{shard_name}.naive.cram \
      -Ou ~{output_prefix}.~{shard_name}.uncertain.cram \
      -umil filtered_umil.bed \
      -R ~{ref_fasta} \
      -EUR ~{expected_uncertain_reads} \
      -FPR ~{bloom_filter_fpr}
    
    echo "Shard ~{shard_name} complete"
    ls -lh ~{output_prefix}.~{shard_name}.*.cram
  >>>

  output {
    File naive_bam     = "~{output_prefix}.~{shard_name}.naive.cram"
    File uncertain_bam = "~{output_prefix}.~{shard_name}.uncertain.cram"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 1
    memory: "6G"
    disks: "local-disk 400 SSD"
    preemptible: 3
  }
}

# Merge multiple CRAM files into one
task merge_crams {
  input {
    Array[File] input_crams
    File   ref_fasta
    File   ref_fasta_fai
    File   ref_dict
    File   gatk_jar
    Int    java_mem
    String output_prefix
    String output_suffix  # "naive" or "uncertain"
  }

  command <<<
    set -e
    echo "MERGING ~{output_suffix} CRAMs ============"
    echo "Input files: ~{sep=' ' input_crams}"
    
    # Use samtools merge for speed (it's in the GATK docker)
    samtools merge \
      -@ 4 \
      --reference ~{ref_fasta} \
      -O CRAM \
      ~{output_prefix}.~{output_suffix}.merged.cram \
      ~{sep=' ' input_crams}
    
    echo "Merge complete"
    ls -lh ~{output_prefix}.~{output_suffix}.merged.cram
  >>>

  output {
    File merged_cram = "~{output_prefix}.~{output_suffix}.merged.cram"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 1
    memory: "6G"
    disks: "local-disk 1000 SSD"
    preemptible: 3
  }
}
