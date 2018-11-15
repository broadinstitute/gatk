
workflow HaplotypeCallerSubWorkflow {

  String picard_docker

  Int num_intervals
  Array[File] intervals

  String bam  
  String base_file_name

  Int disk_size

  File ref_dict = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File ref_fasta = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File ref_fasta_index = "gs://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  # Call variants in parallel over WGS calling intervals
  scatter (index in range(num_intervals)) {
    # Generate GVCF by interval
    call HaplotypeCallerGATK3 {
      input:
        input_bam = bam,
        interval_list = intervals[index],
        gvcf_basename = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = disk_size
     }

    call HaplotypeCallerGATK4 {
      input:
        input_bam = bam,
        interval_list = intervals[index],
        gvcf_basename = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = disk_size
     }

     call MakeVCFIndex {
          input:
            vcf = HaplotypeCallerGATK4.output_gvcf,
            sample = base_file_name,
            disk_size = disk_size,
      }
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs as MergeVCFsGATK3 {
    input:
      input_vcfs = HaplotypeCallerGATK3.output_gvcf,
      input_vcfs_indexes = HaplotypeCallerGATK3.output_gvcf_index,
      output_vcf_name = base_file_name + ".g.vcf.gz",
      disk_size = disk_size,
      picard_docker = picard_docker
  }

  call MergeVCFs as MergeVCFsGATK4 {
    input:
      input_vcfs = HaplotypeCallerGATK4.output_gvcf,
      input_vcfs_indexes = MakeVCFIndex.output_vcf_index,
      output_vcf_name = base_file_name + ".g.vcf.gz",
      disk_size = disk_size,
      picard_docker = picard_docker
  }

  call GenotypeGVCF as GenotypeGVCF3 {
    input:
      gvcf = MergeVCFsGATK3.output_vcf,
      gvcf_index = MergeVCFsGATK3.output_vcf_index,
      output_vcf_name = base_file_name + ".vcf.gz",
      disk_size = disk_size,
      ref = ref_fasta,
      ref_index = ref_fasta_index,
      ref_dict = ref_dict
  }
  call GenotypeGVCF as GenotypeGVCF4 {
    input:
      gvcf = MergeVCFsGATK4.output_vcf,
      gvcf_index = MergeVCFsGATK4.output_vcf_index,
      output_vcf_name = base_file_name + ".vcf.gz",
      disk_size = disk_size,
      ref = ref_fasta,
      ref_index = ref_fasta_index,
      ref_dict = ref_dict
  }

  output {
    File gatk3_vcf = GenotypeGVCF3.output_vcf
    File gatk3_vcf_index = GenotypeGVCF3.output_vcf_index
    File gatk4_vcf = GenotypeGVCF4.output_vcf
    File gatk4_vcf_index = GenotypeGVCF4.output_vcf_index
  }
}


# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGATK3 {
  String input_bam
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int disk_size

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # Using PrintReads is a temporary solution until we update HaploypeCaller to use GATK4. Once that is done,
  # HaplotypeCaller can stream the required intervals directly from the cloud.
  command {
    /gatk/gatk-launch --javaOptions "-Xms2g" \
      PrintReads \
      -I ${input_bam} \
      --interval-padding 500 \
      -L ${interval_list} \
      -O local.sharded.bam \
    && \
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
      -jar /gatk/gatk3/GenomeAnalysisTK_3.8-1-1-gdde23f56a6.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I local.sharded.bam \
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ${default=0 contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}


# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCallerGATK4 {
  String input_bam
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int disk_size


  command {
    /gatk/gatk-launch --javaOptions "-Xms2g" \
      PrintReads \
      -I ${input_bam} \
      --interval-padding 500 \
      -L ${interval_list} \
      -R ${ref_fasta} \
      -O local.sharded.bam \
    && \
    /gatk/gatk-launch --javaOptions "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m" HaplotypeCaller \
      -R ${ref_fasta} \
      -O ${gvcf_basename}.vcf.gz \
      -I local.sharded.bam \
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -contamination ${default=0 contamination} \
      --read-filter "OverclippedReadFilter"
  }

  runtime {
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    # File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi" I think it creates an index now but I don't want to rerun it..
  }
}

task MakeVCFIndex {
    File vcf
    Int disk_size
    String sample
    
    command {
        /gatk/gatk-launch IndexFeatureFile -F ${vcf} -O ${sample}.vcf.gz.tbi
    }

    runtime {
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf_index = "${sample}.vcf.gz.tbi"
    }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  Int disk_size
  String picard_docker

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_vcf_name}
  }
  runtime {
    docker: picard_docker
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task GenotypeGVCF {
  File gvcf
  File gvcf_index
  String output_vcf_name
  Int disk_size
  File ref
  File ref_index
  File ref_dict

  command {
        /gatk/gatk-launch GenotypeGVCFs \
        -R ${ref} \
        --variant ${gvcf} \
        -O ${output_vcf_name}
    }

    runtime {
         memory: "3 GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf = "${output_vcf_name}"
        File output_vcf_index = "${output_vcf_name}.tbi"
    }
}