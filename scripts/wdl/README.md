This directory contains WDL scripts for running GATK CNV workflows.

For more information on WDL and running WDL workflows with cromwell:

https://github.com/broadinstitute/wdl

https://github.com/broadinstitute/cromwell

#### Setting up parameter json file for a run

To get started, copy the relevant ``*.template.json`` for the workflow you wish to run and adjust parameters accordingly.  This file has reasonable default parameters for exome or targeted panels.  
- For WGS, please see notes below for ``isWGS``, ``wgsBinSize``, ``seg_param_undoSplits``, and ``seg_param_undoSD``
- Values starting with ``$__`` *must* be replaced with values for your run.  Please note that python [templates](https://docs.python.org/2/library/string.html#template-strings) can be useful for replacing these values


#### Explanation of fields in PoN workflow

Fields marked as Advanced almost never require a user to adjust the default value.
Reference used must be the same between PoN and case samples.

- ``pon_gatk_workflow.disable_sequence_dictionary_validation`` -- (Advanced) "true".  For end users, the value should not matter.  Developers leave this as true for debugging on smaller references.
- ``pon_gatk_workflow.isWGS`` -- Change to "true" if running WGS samples
- ``pon_gatk_workflow.noQC`` -- (Advanced) Change to "true" to disable quality control of samples and targets for PoN creation
- ``pon_gatk_workflow.max_open_files`` -- (Advanced) "100",  Maximum number of files to combine simultaneously while combining read counts
- ``pon_gatk_workflow.wgsBinSize`` --  Size of bins (in bp) for WGS read counts CLI.  *This must be the same values used for all case samples*.  Ignored if not running WGS.
- ``pon_gatk_workflow.enable_gc_correction`` --  "true" is recommended,  This is a true/false option change to false to bypass GC bias correction step.  All case samples must be run with the same value as used in PoN creation.
- ``pon_gatk_workflow.ref_fasta`` -- Path to reference fasta file. Broad internal: /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta
- ``pon_gatk_workflow.ref_fasta_fai`` --  Path to reference fasta fai file.  Broad internal: /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta.fai
- ``pon_gatk_workflow.pon_entity_id`` --  Name of the final PoN file
- ``pon_gatk_workflow.normal_bam_list`` --  tsv file consisting of corresponding bam and corresponding index files as described in pon_gatk_workflow.wdl
- ``pon_gatk_workflow.ref_fasta_dict`` --  Path to reference dict file.  Broad internal: /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dict
- ``pon_gatk_workflow.target_file`` --  Target file (NOT in bed format) that was used to describe the baits in capture (exome) samples.  Please run ``ConvertBedToTargetFile`` to convert a BED file to a target file.
- ``pon_gatk_workflow.gatk_jar`` -- Absolute path to gatk-protected.jar


#### Explanation of fields in case workflow

Fields marked as Advanced almost never require a user to adjust the default value.
Reference used must be the same between PoN and case samples.

- ``case_gatk_acnv_workflow.input_bam_list`` -- tsv file consisting of corresponding bam files and entity names as described in case_gatk_acnv_workflow.wdl
- ``case_gatk_acnv_workflow.CalculateTargetCoverage.grouping`` -- (Advanced) Always choose "SAMPLE" for this workflow.
- ``case_gatk_acnv_workflow.common_snp_list`` --  SNP list to use in ``GetBayesianHetCoverage``
- ``case_gatk_acnv_workflow.ref_fasta_fai`` --  Path to reference fasta fai file.  Broad internal: /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta.fai 
- ``case_gatk_acnv_workflow.target_file`` --  Same target file (NOT bed format) as was used for the PoN creation.  *This must be the same as was used to generate the PoN*.  Please run ``ConvertBedToTargetFile`` to convert a BED file to a target file.
- ``case_gatk_acnv_workflow.gatk_jar`` --  Absolute path to gatk-protected.jar
- ``case_gatk_acnv_workflow.PadTargets.mem`` -- (Advanced) Number of GB RAM to use when running PadTargets.  Always use 1, unless your target file is really large
- ``case_gatk_acnv_workflow.CalculateTargetCoverage.disable_all_read_filters`` -- (Advanced) Disables all read filters when making coverage profile.  This will probably have unexpected results.  Always use "false".
- ``case_gatk_acnv_workflow.CalculateTargetCoverage.mem`` -- Number of GB RAM to use when running CalculateTargetCoverage.  Most exome files do not need more than 4. 
- ``case_gatk_acnv_workflow.ref_fasta_dict`` --  Path to reference dict file.  Broad internal: /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dict 
- ``case_gatk_acnv_workflow.CalculateTargetCoverage.keepduplicatereads`` -- (Advanced) Whether coverage calculations should count duplicate reads.  Setting this to "true" has been shown (empirically) to improve performance.  "false" is also supported.
- ``case_gatk_acnv_workflow.PoN`` --  Path to GATK PoN file
- ``case_gatk_acnv_workflow.CalculateTargetCoverage.transform `` -- (Advanced) Always "PCOV"
- ``case_gatk_acnv_workflow.ref_fasta`` -- Path to reference fasta file.  Broad internal: ``/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta`` 
- ``case_gatk_acnv_workflow.NormalizeSomaticReadCounts.mem`` -- "4" Number of GB RAM to use when running 
- ``case_gatk_acnv_workflow.PadTargets.pd`` -- (Advanced) "250" has been shown to improve results.  This parameter is number of bases to pad the target file in ``PadTargets``.  *This must be the same as was used to generate the PoN*
- ``case_gatk_acnv_workflow.disable_sequence_dictionary_validation`` -- (Advanced) "true"   For end users, the value should not matter.  Developers leave this as true for debugging on smaller references.
- ``case_gatk_acnv_workflow.plots_dir`` -- (Advanced)  Path to store produced plots for CopyRatio and ACNV
- ``case_gatk_acnv_workflow.call_cnloh_dir`` -- (Advanced)  Path to directory to store CNLoH and splits results
- ``case_gatk_acnv_workflow.enable_gc_correction`` -- "true" is recommended. This is a true/false option change to false to bypass GC bias correction step.  *This must be the same as was used to generate the PoN*
- ``case_gatk_acnv_workflow.isWGS`` -- Change to true if running WGS samples.  You should also change ``case_gatk_acnv_workflow.seg_param_undoSplits`` and ``case_gatk_acnv_workflow.seg_param_undoSD`` to reduce (not eliminate) hyper-segmentation.
- ``case_gatk_acnv_workflow.wgsBinSize`` -- "10000" is recommended.  Size of bins (in bp) for WGS read counts CLI.  *This must be the same as was used to generate the PoN*.  Ignored if not running WGS.
- ``case_gatk_acnv_workflow.stringency`` -- (Advanced) ``GetHetBayesianCoverage`` parameter, default=30. End users should leave this value unchanged.
- ``case_gatk_acnv_workflow.read_depth_threshold`` -- (Advanced) ``GetHetBayesianCoverage`` parameter, default=15. End users should leave this value unchanged.

*For ``PerformSegmentation`` parameters, see pgs. 11-12 of https://bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf*
- ``case_gatk_acnv_workflow.seg_param_alpha`` -- (Advanced)  ``PerformSegmentation`` parameter, default=0.01
- ``case_gatk_acnv_workflow.seg_param_nperm`` -- (Advanced)  ``PerformSegmentation`` parameter, default=10000
- ``case_gatk_acnv_workflow.seg_param_pmethod`` -- (Advanced)  ``PerformSegmentation`` parameter, default=HYBRID
- ``case_gatk_acnv_workflow.seg_param_minWidth`` -- (Advanced)  ``PerformSegmentation`` parameter, default=2
- ``case_gatk_acnv_workflow.seg_param_kmax`` -- (Advanced)  ``PerformSegmentation`` parameter, default=25
- ``case_gatk_acnv_workflow.seg_param_nmin`` -- (Advanced)  ``PerformSegmentation`` parameter, default=200
- ``case_gatk_acnv_workflow.seg_param_eta`` -- (Advanced)  ``PerformSegmentation`` parameter, default=0.05
- ``case_gatk_acnv_workflow.seg_param_trim`` -- (Advanced)  ``PerformSegmentation`` parameter, default=0.025
- ``case_gatk_acnv_workflow.seg_param_undoSplits`` -- ``PerformSegmentation`` parameter, default=NONE   **For WGS: SDUNDO**
- ``case_gatk_acnv_workflow.seg_param_undoPrune`` -- (Advanced)  ``PerformSegmentation`` parameter, default=0.05
- ``case_gatk_acnv_workflow.seg_param_undoSD`` -- ``PerformSegmentation`` parameter, default=3 **For WGS: 2**

