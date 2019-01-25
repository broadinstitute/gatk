import "https://api.firecloud.org/ga4gh/v1/tools/unity-benchmark-development:test-analysis/versions/5/plain-WDL/descriptor" as analysis_workflow
import "https://api.firecloud.org/ga4gh/v1/tools/unity-benchmark-development:test-benchmark/versions/3/plain-WDL/descriptor" as benchmark_workflow

workflow orchestration {
  String message_string
  File reference_dataset

  call analysis_workflow.analysis as aw {
    input:
      message_string = message_string
  }

  call benchmark_workflow.benchmark as bw {
    input:
      reference = reference_dataset,
      analysis = aw.analysis_output
  }

  output {
    File benchmark_data = bw.benchmark_results
  }
}
