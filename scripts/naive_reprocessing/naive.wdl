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
    Int java_mem = 10
    String output_prefix="test"
  }

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

  call validate_1 {
    input:
      naive_bam = naive_processing.naive_bam,
      uncertain_bam = naive_processing.uncertain_bam,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_dict = ref_dict,
      gatk_jar = gatk_jar,
      output_prefix = output_prefix,
      java_mem = java_mem
  }

  call sort_block as naive_sort_block {
    input:
      bam = naive_processing.naive_bam,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_dict = ref_dict,
      gatk_jar = gatk_jar,
      output_prefix = output_prefix,
      output_prefix_append = "naive",
      java_mem = java_mem
  }

  call sort_block as uncertain_sort_block {
    input:
      bam = naive_processing.uncertain_bam,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_dict = ref_dict,
      gatk_jar = gatk_jar,
      output_prefix = output_prefix,
      output_prefix_append = "uncertain",
      java_mem = java_mem
  }

  call validate_2 {
    input:
      naive_sorted_bam = naive_sort_block.sorted_bam,
      uncertain_sorted_bam = uncertain_sort_block.sorted_bam,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_dict = ref_dict,
      gatk_jar = gatk_jar,
      output_prefix = output_prefix,
      java_mem = java_mem
  }

  call validate_input_block {
    input:
      input_bam = input_bam,
      input_bam_idx = input_bam_idx,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_dict = ref_dict,
      gatk_jar = gatk_jar,
      output_prefix = output_prefix,
      java_mem = java_mem
  }

  output {
    File naive_bam                 = naive_processing.naive_bam
    File uncertain_bam             = naive_processing.uncertain_bam

    File naive_sorted_bam          = naive_sort_block.sorted_bam
    File naive_sorted_bai          = naive_sort_block.sorted_bai
    File uncertain_sorted_bam      = uncertain_sort_block.sorted_bam
    File uncertain_sorted_bai      = uncertain_sort_block.sorted_bai
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

task validate_1 {
  input {
    File naive_bam
    File uncertain_bam
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File gatk_jar
    Int java_mem
    String output_prefix
  }

  command <<<
    set -e
    echo "VALIDATE 1 ========================"
    java -jar -Xmx~{java_mem}G ~{gatk_jar} ValidateSamFile --IGNORE INVALID_TAG_NM --IGNORE MISSING_TAG_NM -M SUMMARY -R ~{ref_fasta} -I ~{naive_bam}     > ~{output_prefix}.validate1.naive.summary.txt     || true
    java -jar -Xmx~{java_mem}G ~{gatk_jar} ValidateSamFile --IGNORE INVALID_TAG_NM --IGNORE MISSING_TAG_NM -M SUMMARY -R ~{ref_fasta} -I ~{uncertain_bam} > ~{output_prefix}.validate1.uncertain.summary.txt || true
  >>>

  output {
    File naive_summary     = "~{output_prefix}.validate1.naive.summary.txt"
    File uncertain_summary = "~{output_prefix}.validate1.uncertain.summary.txt"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 2
    memory: "~{java_mem + 2}G"
    disks: "local-disk 100 SSD"
  }
}

task sort_block {
  input {
    File bam
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File gatk_jar
    Int java_mem
    String output_prefix
    String output_prefix_append
  }

  command <<<
    set -e
    echo "SORT ========================"
    java -jar -Xmx~{java_mem}G ~{gatk_jar} SortSam -I ~{bam} -O ~{output_prefix}.~{output_prefix_append}.sorted.cram  -R ~{ref_fasta} --CREATE_INDEX true -SO coordinate
    echo "DONE SORT ========================"
    ls -altrh
    echo "That's all folks!"
  >>>

  output {
    File sorted_bam     = "~{output_prefix}.~{output_prefix_append}.sorted.cram"
    File sorted_bai     = "~{output_prefix}.~{output_prefix_append}.sorted.cram.bai"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 2
    memory: "~{java_mem + 2}G"
    disks: "local-disk 500 SSD"
  }
}

task validate_2 {
  input {
    File naive_sorted_bam
    File uncertain_sorted_bam
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File gatk_jar
    Int java_mem
    String output_prefix
  }

  command <<<
    set -e
    echo "VALIDATE 2 ========================"
    java -jar -Xmx~{java_mem}G ~{gatk_jar} ValidateSamFile --IGNORE INVALID_TAG_NM --IGNORE MISSING_TAG_NM -M SUMMARY -R ~{ref_fasta} -I ~{naive_sorted_bam}     > ~{output_prefix}.validate2.naive.sorted.summary.txt     || true
    java -jar -Xmx~{java_mem}G ~{gatk_jar} ValidateSamFile --IGNORE INVALID_TAG_NM --IGNORE MISSING_TAG_NM -M SUMMARY -R ~{ref_fasta} -I ~{uncertain_sorted_bam} > ~{output_prefix}.validate2.uncertain.sorted.summary.txt || true
  >>>

  output {
    File naive_sorted_summary     = "~{output_prefix}.validate2.naive.sorted.summary.txt"
    File uncertain_sorted_summary = "~{output_prefix}.validate2.uncertain.sorted.summary.txt"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 2
    memory: "~{java_mem + 2}G"
    disks: "local-disk 100 SSD"
  }
}

task validate_input_block {
  input {
    File input_bam
    File   input_bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File gatk_jar
    Int java_mem
    String output_prefix
  }

  command <<<
    set -e
    echo "VALIDATE INPUT ===================="
    java -jar -Xmx~{java_mem}G ~{gatk_jar} ValidateSamFile --IGNORE INVALID_TAG_NM --IGNORE MISSING_TAG_NM -M SUMMARY -R ~{ref_fasta} -I ~{input_bam} > ~{output_prefix}.validate.input.summary.txt || true
  >>>

  output {
    File input_summary = "~{output_prefix}.validate.input.summary.txt"
  }

  runtime {
    docker: "broadinstitute/gatk:latest"
    cpu: 2
    memory: "~{java_mem + 2}G"
    disks: "local-disk 100 HDD"
  }
}
