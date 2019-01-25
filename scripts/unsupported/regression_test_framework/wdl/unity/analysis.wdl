workflow analysis {
  String message_string

  call run_analysis {
    input:
      message = message_string
  }

  output {
    File analysis_output = run_analysis.analysis_output
  }
}


task run_analysis {
  String message

  command {
    echo "${message}" > output_file.txt
    echo "Aligning fastqs and calculating expression" >> output_file.txt
    echo "Created output.transcript.bam and output_matrix.txt" >> output_file.txt
    echo "Creating qc matrix" >> output_file.txt
    echo "Wrote qc matrix" >> output_file.txt
  }
  runtime {
    docker: "ubuntu"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
  output {
    File analysis_output = "output_file.txt"
  }
}