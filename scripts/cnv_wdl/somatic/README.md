## Running the Somatic CNV WDL

### Which WDL should you use?

- Building a panel of normals (PoN): ``cnv_somatic_panel_workflow.wdl``
- Running a matched pair: ``cnv_somatic_pair_workflow.wdl``

#### Setting up parameter json file for a run

To get started, create the json template (using ``java -jar wdltool.jar inputs <workflow>``) for the workflow you wish to run and adjust parameters accordingly.

*Please note that there are optional workflow-level and task-level parameters that do not appear in the template file.  These are set to reasonable values by default, but can also be adjusted if desired.*

#### Required parameters in the somatic panel workflow

The reference used must be the same between PoN and case samples.

- ``CNVSomaticPanelWorkflow.gatk_docker`` -- GATK Docker image (e.g., ``broadinstitute/gatk:latest``).
- ``CNVSomaticPanelWorkflow.intervals`` -- Picard or GATK-style interval list.  For WGS, this should typically only include the autosomal chromosomes.
- ``CNVSomaticPanelWorkflow.normal_bais`` -- List of BAI files.  This list must correspond to `normal_bams`.  For example, `["Sample1.bai", "Sample2.bai"]`.
- ``CNVSomaticPanelWorkflow.normal_bams`` -- List of BAM files.  This list must correspond to `normal_bais`.  For example, `["Sample1.bam", "Sample2.bam"]`.
- ``CNVSomaticPanelWorkflow.pon_entity_id`` -- Name of the final PoN file.
- ``CNVSomaticPanelWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVSomaticPanelWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVSomaticPanelWorkflow.ref_fasta`` -- Path to reference fasta file.

Additionally, there are optional workflow-level and task-level parameters that may be set by advanced users; for example:

- ``CNVSomaticPanelWorkflow.do_explicit_gc_correction`` -- (optional) If true, perform explicit GC-bias correction when creating PoN and in subsequent denoising of case samples.  If false, rely on PCA-based denoising to correct for GC bias.
- ``CNVSomaticPanelWorkflow.PreprocessIntervals.bin_length`` -- Size of bins (in bp) for coverage collection.  Typically, this is set to `0` for capture samples.  *This must be the same value used for all case samples.*
- ``CNVSomaticPanelWorkflow.PreprocessIntervals.padding`` -- Amount of padding (in bp) to add to both sides of targets for WES coverage collection.  *This must be the same value used for all case samples.*

Further explanation of other task-level parameters may be found by invoking the ``--help`` documentation available in the gatk.jar for each tool.

#### Required parameters in the somatic pair workflow

The reference and bins (if specified) must be the same between PoN and case samples.

- ``CNVSomaticPairWorkflow.common_sites`` -- Picard or GATK-style interval list of common sites to use for collecting allelic counts.
- ``CNVSomaticPairWorkflow.gatk_docker`` -- GATK Docker image (e.g., ``broadinstitute/gatk:latest``).
- ``CNVSomaticPairWorkflow.intervals`` -- Picard or GATK-style interval list.  For WGS, this should typically only include the autosomal chromosomes.
- ``CNVSomaticPairWorkflow.normal_bam`` -- (optional, but recommended) Path to normal BAM file.  Do not specify in order to run the tumor sample without a matched normal.
- ``CNVSomaticPairWorkflow.normal_bam_idx`` -- (optional, but recommended) Path to normal BAM file index.  Do not specify in order to run the tumor sample without a matched normal.
- ``CNVSomaticPairWorkflow.read_count_pon`` -- Path to read-count PoN created by the panel workflow. 
- ``CNVSomaticPairWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVSomaticPairWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVSomaticPairWorkflow.ref_fasta`` -- Path to reference fasta file.
- ``CNVSomaticPairWorkflow.tumor_bam`` -- Path to tumor BAM file.
- ``CNVSomaticPairWorkflow.tumor_bam_idx`` -- Path to tumor BAM file index.

Additionally, there are several task-level parameters that may be set by advanced users as above.

To invoke Oncotator on the called tumor copy-ratio segments:

- ``CNVSomaticPairWorkflow.is_run_oncotator`` -- (optional) If true, run Oncotator on the called copy-ratio segments.  This will generate both a simple TSV and a gene list.


Further explanation of these task-level parameters may be found by invoking the ``--help`` documentation available in the gatk.jar for each tool.
