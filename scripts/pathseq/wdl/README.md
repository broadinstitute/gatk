## Running the PathSeq WDL

### Setting up parameter json file for a run

To get started, *copy* the ``*_template.json`` for the workflow you wish to run and modify the parameters accordingly.
DO NOT commit your json file to this repo. This file has reasonable default parameters.

PathSeq reference files are available in the [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle) (direct link [here](ftp://ftp.broadinstitute.org/bundle/beta/PathSeq)).

*Please note that there are optional parameters that do not appear in the template files, since we do not want to specify, by default*

### Docker image
- "broadinstitute/gatk:4.0.0.0"

### Parameter descriptions

Recommended default values (where possible) are found in ``pathseq_pipeline_template.json``

- ``PathSeqPipelineWorkflow.sample_name`` -- sample ID
- ``PathSeqPipelineWorkflow.input_bam`` -- sample BAM file
- ``PathSeqPipelineWorkflow.is_host_aligned`` -- Set to true if the input has already been aligned to a host reference. *NOTE: common human references (e.g. GrCh38) contain decoy sequences such as the Epstein-Barr Virus genome. If this flag is set to true, reads aligning to these decoys will be misidentified as "host" and filtered out.*
- ``PathSeqPipelineWorkflow.min_clipped_read_length`` -- Minimum read length after quality trimming. You may need to reduce this if your input reads are shorter than the default value.
- ``PathSeqPipelineWorkflow.filter_duplicates`` -- If true, reads with identical sequences will be filtered.  
- ``PathSeqPipelineWorkflow.min_score_identity`` -- Fraction of bases aligned to count a microbe alignment in abundance scoring.
- ``PathSeqPipelineWorkflow.identity_margin`` -- If a read maps to more than one microbe, the best alignment is counted. Any additional alignments are also counted if they are within this margin of the best alignment. 
- ``PathSeqPipelineWorkflow.divide_by_genome_length`` -- If true, abundance scores are normalized by genome length.
- ``PathSeqPipelineWorkflow.filter_bwa_image`` -- Path to host BWA index image
- ``PathSeqPipelineWorkflow.kmer_file`` -- Path to host k-mer file
- ``PathSeqPipelineWorkflow.microbe_bwa_image`` -- Path to microbe BWA index image
- ``PathSeqPipelineWorkflow.microbe_fasta`` -- Path to microbe reference FASTA file
- ``PathSeqPipelineWorkflow.microbe_fasta_dict`` -- Path to microbe reference dictionary file
- ``PathSeqPipelineWorkflow.taxonomy_file`` -- Path to PathSeq taxonomy file
- ``PathSeqPipelineWorkflow.gatk_docker`` -- GATK docker image
