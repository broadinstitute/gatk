## Running the Somatic CNV WDL

### Which WDL should you use?
- Building a panel of normals (PoN): ``cnv_somatic_panel_workflow.wdl``
- Running one pair (or one sample in tumor-only mode): ``cnv_somatic_pair_workflow.wdl``
- Create a method (or method configuration) in Workbench or FireCloud: ``cnv_somatic_pair_workflow.wdl`` and set it for entity type of pair (for pairs) or sample (for tumor-only).
- (Advanced) Single-BAM copy-ratio subworkflow: ``cnv_somatic_copy_ratio_bam_workflow.wdl``
- (Advanced) Pair (or sample in tumor-only mode) allele-fraction subworkflow: ``cnv_somatic_allele_fraction_pair_workflow.wdl``

#### Setting up parameter json file for a run

To get started, copy the relevant ``*.template.json`` for the workflow you wish to run and adjust parameters accordingly.  
- Values starting with ``$__`` *must* be replaced with values for your run.  Please note that python [templates](https://docs.python.org/2/library/string.html#template-strings) can be useful for replacing these values.

*Please note that there are task-level parameters that do not appear in the template files.  These are set to reasonable values by default, but can also be adjusted if desired.*

#### Explanation of fields in the somatic panel workflow

The reference used must be the same between PoN and case samples.

- ``CNVSomaticPanelWorkflow.gatk_jar`` -- Absolute path to gatk-protected.jar.
- ``CNVSomaticPanelWorkflow.normal_bams_list`` -- TSV file consisting of corresponding bam and corresponding index files as described in cnv_somatic_panel_workflow.wdl.
- ``CNVSomaticPanelWorkflow.pon_entity_id`` -- Name of the final PoN file.
- ``CNVSomaticPanelWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVSomaticPanelWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVSomaticPanelWorkflow.ref_fasta`` -- Path to reference fasta file.
- ``CNVSomaticPanelWorkflow.targets`` -- (optional) Target file (NOT in bed format) that was used to describe the baits in capture (exome) samples.  Please run ``ConvertBedToTargetFile`` to convert a BED file to a target file.  If provided, then WES workflow will be run; otherwise, WGS workflow will be run.

In additional, there are several task-level parameters that may be set by advanced users; for example:

- ``CNVSomaticPanelWorkflow.CollectCoverage.wgs_bin_size`` -- Size of bins (in bp) for WGS coverage collection.  *This must be the same value used for all case samples.*  Ignored if not running WGS.
- ``CNVSomaticPanelWorkflow.PadTargets.padding`` -- Amount of padding (in bp) to add to both sides of targets for WES coverage collection.  *This must be the same value used for all case samples.*  Ignored if not running WES.

Further explanation of other task-level parameters may be found by invoking the ``--help`` documentation available in the gatk-protected.jar for each tool.  

#### Explanation of fields in the somatic pair workflow

The reference used must be the same between PoN and case samples.

- ``CNVSomaticPairWorkflow.cnv_panel_of_normals`` -- Path to copy-ratio PoN created by the panel workflow. 
- ``CNVSomaticPairWorkflow.common_sites`` -- (optional) List of common SNP sites to use in ``GetBayesianHetCoverage``.  If not provided, the allele-fraction subworkflow will not be run.
- ``CNVSomaticPairWorkflow.gatk_jar`` -- Absolute path to gatk-protected.jar.
- ``CNVSomaticPairWorkflow.normal_bam_idx`` -- (optional, but required if ``normal_bam`` is provided)  File path or storage location (depending on backend) of the normal BAM file index.
- ``CNVSomaticPairWorkflow.normal_bam`` -- (optional, but required if ``normal_bam_index``  is provided)  File path or storage location (depending on backend) of the normal BAM file.
- ``CNVSomaticPanelWorkflow.ref_fasta_dict`` -- Path to reference dict file.
- ``CNVSomaticPanelWorkflow.ref_fasta_fai`` -- Path to reference fasta fai file.
- ``CNVSomaticPanelWorkflow.ref_fasta`` -- Path to reference fasta file.
- ``CNVSomaticPairWorkflow.targets`` -- (optional) Target file (NOT in bed format) that was used to describe the baits in capture (exome) samples.  Please run ``ConvertBedToTargetFile`` to convert a BED file to a target file.  If provided, then WES workflow will be run; otherwise, WGS workflow will be run.
- ``CNVSomaticPairWorkflow.tumor_bam_idx`` -- File path or storage location (depending on backend) of the tumor BAM file index.
- ``CNVSomaticPairWorkflow.tumor_bam`` -- File path or storage location (depending on backend) of the tumor BAM index.

In additional, there are several task-level parameters (for tasks in the subworkflows ``CNVSomaticPairWorkflow.NormalCopyRatioWorkflow``, ``CNVSomaticPairWorkflow.TumorCopyRatioWorkflow``, and ``CNVSomaticPairWorkflow.TumorAlleleFractionWorkflow``) that may be set by advanced users as above.  Be sure to set copy-ratio parameters identically for both the normal and tumor when appropriate.

Further explanation of these task-level parameters may be found by invoking the ``--help`` documentation available in the gatk-protected.jar for each tool.

*For the ``PerformSegmentation`` task-level parameters, see pgs. 11-12 of https://bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf*