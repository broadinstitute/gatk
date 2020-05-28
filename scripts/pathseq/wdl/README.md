## Running the PathSeq WDL

### Setting up parameter json file for a run

To get started, *copy* the ``pathseq_pipeline_template.json`` for the workflow and modify the parameters accordingly.
DO NOT commit your json file to this repo.

PathSeq reference files are available in the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) at `gs://gatk-best-practices/pathseq`.

*Please note that there are optional parameters that do not appear in the template files, since we do not want to specify them by default*

### Parameter descriptions

Example values (where possible) are found in ``pathseq_pipeline_template.json``

- ``PathSeqPipeline.sample_name`` -- Sample ID
- ``PathSeqPipeline.input_bam_or_cram`` -- Sample BAM or CRAM file
- ``PathSeqPipeline.input_bam_or_cram_index`` -- Sample BAM or CRAM index
- ``PathSeqPipeline.downsample`` -- Set true to downsample input. It is strongly recommended that there are <10M total microbial reads.
- ``PathSeqPipeline.Downsample.reads_after_downsampling`` - Target number of reads if downsampling is enabled
- ``PathSeqPipeline.PathSeqFilter.is_host_aligned`` -- Set to true if the input has already been aligned to a host reference. *NOTE: some human references (e.g. hg38) contain decoy sequences such as Epstein-Barr virus. If this flag is set to true, use `ignored_alignment_contigs` as necessary to prevent reads aligning to these decoys from being misidentified as "host".*
- ``PathSeqPipeline.PathSeqFilter.filter_bwa_image`` -- Path to host BWA index image (may be omitted to skip). This corresponds to `pathseq_host.fa.img` in the Resource Bundle.
- ``PathSeqPipeline.PathSeqFilter.kmer_file`` -- Path to host k-mer file (may be omitted to skip). This corresponds to `pathseq_host.bfi` in the Resource Bundle.
- ``PathSeqPipeline.PathSeqAlign.microbe_bwa_image`` -- Path to microbe BWA index image. This corresponds to `pathseq_microbe.fa.img` in the Resource Bundle.
- ``PathSeqPipeline.PathSeqAlign.microbe_dict`` -- Path to microbe reference dictionary file. This corresponds to `pathseq_microbe.dict` in the Resource Bundle.
- ``PathSeqPipeline.PathSeqScore.taxonomy_file`` -- Path to PathSeq taxonomy file. This corresponds to `pathseq_taxonomy.db` in the Resource Bundle.
- ``PathSeqPipeline.gatk_docker`` -- GATK docker image ( > v4.1.7.0)

Optional parameters:

- ``PathSeqPipeline.ignored_alignment_contigs`` -- Array of contig names that will be ignored when filtering pre-aligned reads (see `is_host_aligned` above), e.g. "chrEBV".
- ``PathSeqPipeline.PathSeqFilter.min_clipped_read_length`` -- Minimum read length after quality trimming. Increasing may increase microbial classification specificity but may reduce sensitivity (default 31)
- ``PathSeqPipeline.PathSeqFilter.filter_duplicates`` -- If true, reads with identical sequences will be filtered (default true)
- ``PathSeqPipeline.PathSeqFilter.skip_quality_filters`` -- If true, skips quality filtering steps (default false)
- ``PathSeqPipeline.PathSeqFilter.skip_pre_bwa_repartition`` -- Advanced tuning parameter; recommended if the sample contains a large number (e.g. >10M) microbial reads (default false)
- ``PathSeqPipeline.PathSeqFilter.filter_bwa_seed_length`` -- Sets host alignment MEM length. Reducing this value increases host read detection sensitivity but is slower (default 19) 
- ``PathSeqPipeline.PathSeqFilter.host_min_identity`` -- Minimum identity score for reads to be classified as host (default 30)
- ``PathSeqPipeline.PathSeqScore.min_score_identity`` -- Fraction of bases aligned to count a microbe alignment in abundance scoring (default 0.9)
- ``PathSeqPipeline.PathSeqScore.identity_margin`` -- If a read maps to more than one microbe, the best alignment is counted. Any additional alignments are also counted if they are within this margin of the best alignment (default 0.02)
- ``PathSeqPipeline.PathSeqScore.divide_by_genome_length`` -- If true, abundance scores are normalized by genome length (default false)
- ``PathSeqPipeline.PathSeqScore.not_normalized_by_kingdom`` -- If true, normalized scores are tallied by superkingdom (default false)
- ``PathSeqPipeline.*_mem_gb`` -- Virtual machine (VM) memory in GB (default 208 GB)
- ``PathSeqPipeline.*_preemptible_attempts`` -- Minimum number of attempts on preemtible VMs before the job will fail (default 3)
- ``PathSeqPipeline.*_additional_disk_gb`` -- VM disk space padding in GB
- ``PathSeqPipeline.*_cpu`` -- Number of VM virtual processor cores (default 32)
- ``PathSeqPipeline.*_ssd`` -- Use solid-state disks (default false)
