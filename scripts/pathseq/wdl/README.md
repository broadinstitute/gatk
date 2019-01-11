## Running the PathSeq WDL

### Setting up parameter json file for a run

To get started, *copy* the ``pathseq_pipeline_template.json`` for the workflow and modify the parameters accordingly.
DO NOT commit your json file to this repo. This file has reasonable default parameters.

PathSeq reference files are available in the [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle) (located [here](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/beta/PathSeq)).

*Please note that there are optional parameters that do not appear in the template files, since we do not want to specify them by default*

### Docker image
- "broadinstitute/gatk:4.0.0.0"

### Parameter descriptions

Recommended default values (where possible) are found in ``pathseq_pipeline_template.json``

- ``PathSeqPipelineWorkflow.sample_name`` -- sample ID
- ``PathSeqPipelineWorkflow.input_bam`` -- sample BAM file
- ``PathSeqPipelineWorkflow.is_host_aligned`` -- Set to true if the input has already been aligned to a host reference. *NOTE: common human references (e.g. GrCh38) contain decoy sequences such as the Epstein-Barr Virus genome. If this flag is set to true, reads aligning to these decoys will be misidentified as "host" and filtered out.*
- ``PathSeqPipelineWorkflow.filter_bwa_image`` -- Path to host BWA index image. This corresponds to `pathseq_host.fa.img` in the Resource Bundle.
- ``PathSeqPipelineWorkflow.kmer_file`` -- Path to host k-mer file. This corresponds to `pathseq_host.bfi` in the Resource Bundle.
- ``PathSeqPipelineWorkflow.microbe_bwa_image`` -- Path to microbe BWA index image. This corresponds to `pathseq_microbe.fa.img` in the Resource Bundle.
- ``PathSeqPipelineWorkflow.microbe_fasta`` -- Path to microbe reference FASTA file. This corresponds to `pathseq_microbe.fa` in the Resource Bundle.
- ``PathSeqPipelineWorkflow.microbe_fasta_dict`` -- Path to microbe reference dictionary file.This corresponds to `pathseq_microbe.dict` in the Resource Bundle.
- ``PathSeqPipelineWorkflow.taxonomy_file`` -- Path to PathSeq taxonomy file. This corresponds to `pathseq_taxonomy.db` in the Resource Bundle.
- ``PathSeqPipelineWorkflow.gatk_docker`` -- GATK docker image

Optional parameters:

- ``PathSeqPipelineWorkflow.min_clipped_read_length`` -- Minimum read length after quality trimming. You may need to reduce this if your input reads are shorter than the default value (default 60)
- ``PathSeqPipelineWorkflow.min_score_identity`` -- Fraction of bases aligned to count a microbe alignment in abundance scoring (default 0.9)
- ``PathSeqPipelineWorkflow.identity_margin`` -- If a read maps to more than one microbe, the best alignment is counted. Any additional alignments are also counted if they are within this margin of the best alignment (default 0.02)
- ``PathSeqPipelineWorkflow.divide_by_genome_length`` -- If true, abundance scores are normalized by genome length (default false)
- ``PathSeqPipelineWorkflow.filter_duplicates`` -- If true, reads with identical sequences will be filtered (default true)
- ``PathSeqPipelineWorkflow.skip_quality_filters`` -- If true, skips quality filtering steps (default false)
- ``PathSeqPipelineWorkflow.skip_pre_bwa_repartition`` -- Advanced tuning parameter; recommended if the sample contains a large number (e.g. >10M) microbial reads (default false)
- ``PathSeqPipelineWorkflow.filter_bwa_seed_length`` -- Sets host alignment MEM length. Reducing this valueincreases host read detection sensitivity but is slower (default 19) 
- ``PathSeqPipelineWorkflow.host_min_identity`` -- Minimum identity score for reads to be classified as host (default 30)
- ``PathSeqPipelineWorkflow.bam_partition_size`` -- Advanced tuning parameter; set size of Spark read partitions. Reducing this value increases parallelism but may also increase overhead and affect host BWA alignment (default 4000000)
- ``PathSeqPipelineWorkflow.mem_gb`` -- Virtual machine (VM) memory in GB (default 208 GB)
- ``PathSeqPipelineWorkflow.preemptible_attempts`` -- Minimum number of attempts on preemtible VMs before the job will fail (default 3)
- ``PathSeqPipelineWorkflow.disk_space_gb`` -- VM disk space in GB (automatically sized by default)
- ``PathSeqPipelineWorkflow.cpu`` -- Number of VM virtual processor cores (default 32)
