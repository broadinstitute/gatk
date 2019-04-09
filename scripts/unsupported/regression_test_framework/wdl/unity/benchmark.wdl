workflow benchmark {
  File reference
  File analysis

  call run_benchmark {
    input:
      reference_input = reference,
      analysis_input = analysis
  }

  output {
    File benchmark_results = run_benchmark.benchmark_results
  }
}

task run_benchmark {
  File reference_input
  File analysis_input

  command {
    diff ${reference_input} ${analysis_input} >> benchmark.txt || :
  }
  runtime {
    docker: "ubuntu"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
  output {
    File benchmark_results = "benchmark.txt"
  }
}


