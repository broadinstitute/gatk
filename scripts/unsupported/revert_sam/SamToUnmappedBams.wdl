workflow CramToUnmappedBams {
  File? input_bam
  File? input_cram
  File? ref_fasta
  File? ref_fasta_index
  String dir_pattern = "gs://.*/"
  Int? cram_to_bam_disk_size
  Int revert_sam_disk_size
  Int sort_sam_disk_size
  Int validate_sam_file_disk_size

  if(defined(input_cram)) {
    String input_cram_path = select_first([input_cram])
    String bam_from_cram_name = basename(input_cram_path, ".cram") + ".bam"

    call CramToBam {
      input:
        input_cram = select_first([input_cram]),
        output_bam_name = bam_from_cram_name,
        ref_fasta = select_first([ref_fasta]),
        ref_fasta_index = select_first([ref_fasta_index]),
        disk_size = select_first([cram_to_bam_disk_size])
    }
  }

  File input_data = select_first([CramToBam.output_bam, input_bam])

  call GenerateOutputMap {
    input:
      input_bam = input_data,
      disk_size = revert_sam_disk_size
  }

  call RevertSam {
    input:
      input_bam = input_data,
      output_map = GenerateOutputMap.output_map,
      disk_size = revert_sam_disk_size
  }

  scatter (unmapped_bam in RevertSam.unmapped_bams) {
    String output_basename = sub(sub(unmapped_bam, dir_pattern, ""), ".coord.sorted.unmapped_bam$", "")

    call SortSam {
      input:
        input_bam = unmapped_bam,
        sorted_bam_name = output_basename + ".unmapped.bam",
        disk_size = sort_sam_disk_size
    }
    call ValidateSamFile {
      input:
        input_bam = SortSam.sorted_bam,
        report_filename = output_basename + ".validation_report",
        disk_size = validate_sam_file_disk_size
    }
  }

  output {
    Array[File] sortsam_out = SortSam.sorted_bam
    Array[File] validatesam_out = ValidateSamFile.report
  }
}

task GenerateOutputMap {
  File input_bam
  Int disk_size

  command {
    samtools view -H ${input_bam} | grep @RG | cut -f2 | tr -d ID: > readgroups.txt

    echo -e "READ_GROUP_ID\tOUTPUT" > output_map.tsv
    
    for rg in `cat readgroups.txt`; do
      echo -e "$rg\t$rg.unmapped.bam" >> output_map.tsv
    done
  }

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    disks: "local-disk " + disk_size + " HDD"
    memory: "1 GB"
  }
  output {
    File output_map = "output_map.tsv"
  }
}

task CramToBam {
  File input_cram
  String output_bam_name
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  command {
    samtools view -1 -h -T ${ref_fasta} -o ${output_bam_name} ${input_cram}
  }

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    disks: "local-disk " + disk_size + " HDD"
    memory: "1 GB"
  }
  output {
    File output_bam = "${output_bam_name}"
  }
}

task RevertSam {
  File input_bam
  File output_map
  Int disk_size

  command {
    java -Xmx1000m -jar /usr/gitc/picard.jar \
    RevertSam \
    INPUT=${input_bam} \
    OUTPUT_MAP=${output_map} \
    OUTPUT_BY_READGROUP=true \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    SORT_ORDER=coordinate
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    disks: "local-disk " + disk_size + " HDD"
    memory: "1200 MB"
  }
  output {
    Array[File] unmapped_bams = glob("*.bam")
  }
}

task SortSam {
  File input_bam
  String sorted_bam_name
  Int disk_size

  command {
    java -Xmx3000m -jar /usr/gitc/picard.jar \
    SortSam \
    INPUT=${input_bam} \
    OUTPUT=${sorted_bam_name} \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=1000000
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MB"
    preemptible: 3
  }
  output {
    File sorted_bam = "${sorted_bam_name}"
  }
}

task ValidateSamFile {
  File input_bam
  String report_filename
  Int disk_size

  command {
    java -Xmx3000m -jar /usr/gitc/picard.jar \
      ValidateSamFile \
      INPUT=${input_bam} \
      OUTPUT=${report_filename} \
      MODE=VERBOSE \
      IS_BISULFITE_SEQUENCED=false 
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.2.3-1469027018"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MB"
    preemptible: 3
  }
  output {
    File report = "${report_filename}"
  }
}