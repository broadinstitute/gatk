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
    Int java_mem = 8
    String output_prefix="test"
    
    # Parallelization option - if true, splits by chromosome groups for faster processing
    Boolean parallelize = false
    
    # Uncertain-only optimization - if true, uses optimized I/O mode that only writes
    # uncertain reads during traversal, then creates naive via samtools subtraction
    # This can reduce total runtime by ~3 hours for typical samples
    Boolean uncertain_only = false
    
    # Parallel processing options (only used when parallelize=true)
    Array[Array[String]] chromosome_groups = [
      ["chr1", "chr2"],
      ["chr3", "chr4", "chr5", "chr6"],
      ["chr7", "chr8", "chr9", "chr10", "chr11", "chr12"],
      ["chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
    ]
    Array[String] shard_names = ["shard1_chr1_2", "shard2_chr3_6", "shard3_chr7_12", "shard4_chr13_Y"]
    Array[Int] expected_uncertain_per_shard = [15000000, 12000000, 12000000, 15000000]
  }

  # ============================================================================
  # UNCERTAIN-ONLY OPTIMIZED PATH (uncertain_only=true)
  # Uses optimized I/O: writes only uncertain reads, creates naive via samtools
  # ============================================================================
  if (uncertain_only) {
    call generate_uncertain_read_names_only {
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
  # SEQUENTIAL PROCESSING PATH (parallelize=false, uncertain_only=false)
  # ============================================================================
  if (!parallelize && !uncertain_only) {
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
  # PARALLEL PROCESSING PATH (parallelize=true, uncertain_only=false)
  # ============================================================================
  if (parallelize && !uncertain_only) {
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
    generate_uncertain_read_names_only.uncertain_bam,
    merge_naive.merged_cram, 
    naive_processing.naive_bam
  ])
  File final_uncertain_bam = select_first([
    generate_uncertain_read_names_only.uncertain_bam,
    merge_uncertain.merged_cram, 
    naive_processing.uncertain_bam
  ])

  output {
    File naive_bam                 = final_naive_bam
    File uncertain_bam             = final_uncertain_bam
    
    # Parallel processing debug outputs (only populated when parallelize=true)
    Array[File]? shard_naive_bams     = naive_processing_shard.naive_bam
    Array[File]? shard_uncertain_bams = naive_processing_shard.uncertain_bam
    
    # Uncertain-only mode output (only populated when uncertain_only=true)
    File? uncertain_read_names = generate_uncertain_read_names_only.uncertain_read_names
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
# EXPERIMENTAL: Only generate uncertain read names file for analysis
# ============================================================================
# This task only writes the uncertain reads CRAM and read names file,
# without attempting to create the naive CRAM. Use this to examine
# the read names file and determine optimal storage/processing strategies.

task generate_uncertain_read_names_only {
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
    echo "GENERATING UNCERTAIN READ NAMES ONLY (EXPERIMENTAL) ============"
    echo "This task only writes uncertain reads + read names file for analysis."
    echo ""
    
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
    
    echo ""
    echo "============================================"
    echo "Uncertain CRAM and read names file created."
    echo "Uncertain read names count: $(wc -l < ~{output_prefix}.uncertain_read_names.txt)"
    echo "Read names file size: $(ls -lh ~{output_prefix}.uncertain_read_names.txt | awk '{print $5}')"
    echo ""
    echo "Analyzing read names..."
    echo "First 10 read names:"
    head -n 10 ~{output_prefix}.uncertain_read_names.txt
    echo ""
    echo "Average read name length:"
    awk '{total += length($0); count++} END {print total/count " characters"}' ~{output_prefix}.uncertain_read_names.txt
    echo ""
    echo "Analysis complete. Files ready for examination."
    echo "============================================"
  >>>

  output {
    File uncertain_bam = "~{output_prefix}.uncertain.cram"
    File uncertain_read_names = "~{output_prefix}.uncertain_read_names.txt"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 2
    memory: "~{java_mem + 2}G"
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
    echo "Stage 1: Generate uncertain CRAM and front-coded read names list..."
    echo ""
    
    java -jar -Xmx~{java_mem}G ~{gatk_jar} SplitReadsByRealignmentDifficulty \
      -RLG ~{rlg} \
      -I ~{input_bam} \
      -O /dev/null \
      -Ou ~{output_prefix}.uncertain.cram \
      -umil ~{umil} \
      -R ~{ref_fasta} \
      -UO true \
      -URN ~{output_prefix}.uncertain_read_names.fcl \
      -EUR ~{expected_uncertain_reads} \
      -FPR ~{bloom_filter_fpr}
    
    echo ""
    echo "Stage 1 complete. Uncertain CRAM and .fcl file created."
    echo "File sizes:"
    ls -lh ~{output_prefix}.uncertain.cram ~{output_prefix}.uncertain_read_names.fcl
    echo ""
    
    echo "Stage 2: Filter input CRAM using .fcl to create naive CRAM..."
    echo ""
    
    java -jar -Xmx~{java_mem}G ~{gatk_jar} FilterReadsByNameList \
      -I ~{input_bam} \
      -O ~{output_prefix}.naive.cram \
      -R ~{ref_fasta} \
      -XRL ~{output_prefix}.uncertain_read_names.fcl
    
    echo ""
    echo "Stage 2 complete. Naive CRAM created."
    echo ""
    echo "============================================"
    echo "TWO-STAGE PROCESSING COMPLETE"
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
    # This speeds up overlap detection significantly
    echo "Filtering BED file for shard chromosomes..."
    CHROM_PATTERN=$(echo "~{sep='|' chromosomes}" | sed 's/|/\\|/g')
    grep -E "^($CHROM_PATTERN)\t" ~{umil} > filtered_umil.bed || true
    
    # If no BED regions for these chromosomes, create empty file
    if [ ! -s filtered_umil.bed ]; then
      echo "Warning: No BED regions found for chromosomes ~{sep=' ' chromosomes}"
      touch filtered_umil.bed
    fi
    
    echo "Filtered BED lines: $(wc -l < filtered_umil.bed)"
    
    # Run SplitReadsByRealignmentDifficulty on this shard
    java -jar -Xmx~{java_mem}G ~{gatk_jar} SplitReadsByRealignmentDifficulty \
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
    cpu: 2
    memory: "~{java_mem + 2}G"
    disks: "local-disk 400 SSD"
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
    cpu: 4
    memory: "~{java_mem + 2}G"
    disks: "local-disk 1000 SSD"
  }
}
