## Running the M2 WDL

### Which WDL should you use for M2?
- Running one pair (or one sample in tumor-only mode): ``mutect2.wdl``
- Running several pairs (or multiple samples in tumor-only mode): ``mutect2_multi_sample.wdl``.
- Create a method (or method configuration ) in Workbench or FireCloud: ``mutect2.wdl`` and set it for entity type of pair (for pairs) or sample (for tumor-only).


### Setting up parameter json file for a run

To get started, *copy* the relevant ``*_template.json`` for the workflow you wish to run and adjust parameters accordingly.  DO NOT change the sample json file, nor should you commit your json file to this repo.
This file has reasonable default parameters.
- Values starting with ``$__`` *must* be replaced with values for your run.  Please note that python [templates](https://docs.python.org/2/library/string.html#template-strings) can be useful for replacing these values

*Please note that there are optional parameters that do not appear in the template files, since we do not want to specify, by default*

### Docker images
- "broadinstitute/gatk-protected:1.0.0.0-alpha1.2.4" (This is a private image!  Recommended use ``gatk_jar`` as ``/root/gatk-protected.jar``)
- "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282" (You must specify a ``gatk4_jar_override``)


### Parameter descriptions

#### mutect2_multi_sample

Recommended default values (where possible) are found in ``mutect2_multi_sample_template.json``

- ``Mutect2_Multi.gatk4_jar`` -- Location *within the docker file* of the GATK4 jar file.  If you wish you to use a different jar file, such as one on your local filesystem or a google bucket, specify that location with ``Mutect2_Multi.gatk4_jar_override``.  This parameter is ignored if ``Mutect2_Multi.gatk4_jar_override`` is specified.
- ``Mutect2_Multi.scatter_count`` -- Number of executions to split the Mutect2 task into.  The more you put here, the faster Mutect2 will return results, but at a higher cost of resources.
- ``Mutect2_Multi.intervals`` -- A file listing genomic intervals to search for somatic mutations.  This should be in the standard GATK4 format.
- ``Mutect2_Multi.ref_fasta`` -- reference fasta.  For Broad internal VM:  ``/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta``
- ``Mutect2_Multi.ref_fasta_index`` -- For Broad internal VM:  ``/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta.fai``
- ``Mutect2_Multi.ref_dict`` -- For Broad internal VM:  ``/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dict``
- ``Mutect2_Multi.pon`` -- (optional) Panel of normals VCF to use for false positive reduction.
- ``Mutect2_Multi.pon_index`` -- (optional, but required if ``Mutect2_Multi.pon`` is specified)  VCF index for the panel of normals.  Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2_Multi.dbsnp`` -- (optional)  dbSNP VCF file for pruning of false positives.
- ``Mutect2_Multi.dbsnp_index`` -- (optional, but required if ``Mutect2_Multi.dbsnp`` is specified)  VCF index for dbSNP.  Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2_Multi.cosmic`` -- (optional)  COSMIC VCF file to inform potential true positives.
- ``Mutect2_Multi.cosmic_index`` -- (optional, but required if ``Mutect2_Multi.cosmic`` is specified)  VCF index for COSMIC.  Please see GATK4 tool ``IndexFeatureFile`` for creation of an index.
- ``Mutect2_Multi.is_run_orientation_bias_filter`` -- ``true``/``false`` whether the orientation bias filter should be run.
- ``Mutect2_Multi.is_run_oncotator`` -- ``true``/``false`` whether the command-line version of oncotator should be run.  If ``false``, ``Mutect2_Multi.oncotator_docker`` parameter is ignored.
- ``Mutect2_Multi.m2_docker`` -- Docker image to use for Mutect2 tasks.  This is only used for backends configured to use docker.
- ``Mutect2_Multi.oncotator_docker`` -- Docker image to use for Oncotator tasks.  This is only used for backends configured to use docker.
- ``Mutect2_Multi.gatk4_jar_override`` -- (optional)  A GATK4 jar file to be used instead of the jar file in the docker image.  (See ``Mutect2_Multi.gatk4_jar``)  This can be very useful for developers.  Please note that you need to be careful that the docker image you use is compatible with the GATK4 jar file given here -- no automated checks are made.
- ``Mutect2_Multi.preemptible_attempts`` -- Number of times to attempt running a task on a preemptible VM.  This is only used for cloud backends in cromwell and is ignored for local and SGE backends.
- ``Mutect2_Multi.artifact_modes`` -- List of artifact modes to search for in the orientation bias filter.  For example to filter the OxoG artifact, you would specify ``["G/T"]``.  For both the FFPE artifact and the OxoG artifact, specify ``["G/T", "C/T"]``.  If you do not wish to search for any artifacts, please set ``Mutect2_Multi.is_run_orientation_bias_filter`` to ``false``.
- ``Mutect2_Multi.onco_ds_tar_gz`` -- (optional)  A tar.gz file of the oncotator datasources -- often quite large (>15GB).  This will be uncompressed as part of the oncotator task.  Depending on backend used, this can be specified as a path on the local filesystem of a cloud storage container (e.g. gs://...).  Typically the Oncotator default datasource can be downloaded at ``ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/``.  Do not put the FTP URL into the json file.
- ``Mutect2_Multi.onco_ds_local_db_dir`` -- (optional)  A direct path to the Oncotator datasource directory (uncompressed).  While this is the fastest approach, it cannot be used with docker unless your docker image already has the datasources in it.  For cromwell backends without docker, this can be a local filesystem path.  *This cannot be a cloud storage location*
- ``Mutect2_Multi.picard_jar`` -- A direct path to a picard jar for using ``CollectSequencingArtifactMetrics``.  This parameter requirement will be eliminated in the future.

 Note:  If neither ``Mutect2_Multi.onco_ds_tar_gz`` nor ``Mutect2_Multi.onco_ds_local_db_dir`` are specified, the Oncotator task will download and uncompress for each execution.

- ``Mutect2_Multi.pair_list`` -- a tab-separated table with no header in the following format:
 ```
 TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE</TAB>NORMAL_1_BAM</TAB>NORMAL_1_BAM_INDEX</TAB>NORMAL_1_SAMPLE</TAB>
 TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE</TAB>NORMAL_2_BAM</TAB>NORMAL_2_BAM_INDEX</TAB>NORMAL_2_SAMPLE</TAB>
 . . .
 ```
Temporarily, while waiting for a Cromwell bug to be resolved, tumor-only input looks like:
```
TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>
TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>
. . .
```
 That is, you actually write out "NO NORMAL" thrice per row

#### mutect2 (single pair/sample)

Recommended default values (where possible) are found in ``mutect2_template.json``

- ``Mutect2.gatk4_jar`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.intervals`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.ref_fasta`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.ref_fasta_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.ref_dict`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.tumor_bam`` -- File path or storage location (depending on backend) of the tumor bam file.
- ``Mutect2.tumor_bam_index`` --  File path or storage location (depending on backend) of the tumor bam file index.
- ``Mutect2.tumor_sample_name`` -- A name to identify the tumor sample being used.
- ``Mutect2.normal_bam`` -- (optional) File path or storage location (depending on backend) of the normal bam file.
- ``Mutect2.normal_bam_index`` --  (optional, but required if ``Mutect2.normal_bam`` is specified)  File path or storage location (depending on backend) of the normal bam file index.
- ``Mutect2.normal_sample_name`` --  (optional, but required if ``Mutect2.normal_bam`` is specified)  A name to identify the normal sample being used.
- ``Mutect2.pon`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.pon_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.scatter_count`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.dbsnp`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.dbsnp_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.cosmic`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.cosmic_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.is_run_orientation_bias_filter`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.is_run_oncotator`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.m2_docker`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.oncotator_docker`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.gatk4_jar_override`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.preemptible_attempts`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.onco_ds_tar_gz`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.onco_ds_local_db_dir`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.artifact_modes`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2.picard_jar`` -- Please see parameter description above in the mutect2_multi_sample.

#### mutect2-replicate-validation

*This script is usually used by developers and evaluators only.*

Recommended default values (where possible) are found in ``mutect2-replicate-validation_template.json``

- ``Mutect2ReplicateValidation.cosmic_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.dbsnp`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.gatk4_jar`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.is_run_orientation_bias_filter`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.Mutect2.onco_ds_tar_gz`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.scatter_count`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.intervals`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.is_run_oncotator`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.preemptible_attempts`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.oncotator_docker`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.pon`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.dbsnp_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.pon_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.m2_docker`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.cosmic`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.ref_fasta`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.artifact_modes`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.ref_fasta_index`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.Mutect2.onco_ds_local_db_dir`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.gatk4_jar_override`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.ref_dict`` -- Please see parameter description above in the mutect2_multi_sample.
- ``Mutect2ReplicateValidation.replicate_pair_list`` -- tab-separated values with six columns in the following format:
 ```
 REP_1_BAM</TAB>REP_1_BAM_INDEX</TAB>REP_1_SAMPLE</TAB>REP_2_BAM</TAB>REP_2_BAM_INDEX</TAB>REP_2_SAMPLE</TAB>
 REP_3_BAM</TAB>REP_3_BAM_INDEX</TAB>REP_3_SAMPLE</TAB>REP_4_BAM</TAB>REP_4_BAM_INDEX</TAB>REP_4_SAMPLE</TAB>
 . . .
 ```
For example:

```
gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V3.bam    gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V3.bai    SM-612V3        gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V4.bam gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V4.bai    SM-612V4
gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V3.bam    gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V3.bai    SM-612V3        gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V5.bam gs://broad-dsde-methods/takuto/na12878-crsp-ice/SM-612V5.bai    SM-612V5
```
### Example json

#### mutect2_multi_sample

- Local backend with docker
- Cromwell 0.25
- ``mutect2_multi_sample.wdl`` (though note this was actually only run on one pair)
- Uses a local copy of the tar.gz file for Oncotator.  This saves time downloading the .tar.gz file.
- Runs both the orientation bias filter and Oncotator.
- This puts in a dummy value for ``gatk4_jar``.  That makes sure that the correct jar is run, though this change is optional.


*You will need to change the values of the parameters.*  This json is just provided for illustration.

```
{
  "Mutect2_Multi.gatk4_jar_override": "/home/lichtens/test_onco_m2/gatk-protected/build/libs/gatk-protected.jar",
  "Mutect2_Multi.gatk4_jar": "DO NOT USE",
  "Mutect2_Multi.intervals": "/home/lichtens/test_onco_m2/intervals.interval_list",
  "Mutect2_Multi.ref_fasta": "/data/ref/Homo_sapiens_assembly19.fasta",
  "Mutect2_Multi.ref_fasta_index": "/data/ref/Homo_sapiens_assembly19.fasta.fai",
  "Mutect2_Multi.ref_dict": "/data/ref/Homo_sapiens_assembly19.dict",
  "Mutect2_Multi.pair_list": "/home/lichtens/test_onco_m2/pair_list",
  "Mutect2_Multi.pon": "/data/m1/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf",
  "Mutect2_Multi.pon_index": "/data/m1/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf.idx",
  "Mutect2_Multi.scatter_count": 2,
  "Mutect2_Multi.dbsnp": "/data/m1/dbsnp_134_b37.leftAligned.vcf",
  "Mutect2_Multi.dbsnp_index": "/data/m1/dbsnp_134_b37.leftAligned.vcf.idx",
  "Mutect2_Multi.cosmic": "/data/m1/hg19_cosmic_v54_120711.vcf",
  "Mutect2_Multi.cosmic_index": "/data/m1/hg19_cosmic_v54_120711.vcf.idx",
  "Mutect2_Multi.is_run_orientation_bias_filter": true,
  "Mutect2_Multi.is_run_oncotator": true,
  "Mutect2_Multi.m2_docker": "broadinstitute/gatk-protected:1.0.0.0-alpha1.2.4",
  "Mutect2_Multi.oncotator_docker": "broadinstitute/oncotator:1.9.2.0",
  "Mutect2_Multi.preemptible_attempts": 2,
  "Mutect2_Multi.onco_ds_tar_gz": "/data/onco_dir/oncotator_v1_ds_April052016.tar.gz"
}

```

Associated pair_list file (tab separated):

```
/home/lichtens/test_onco_m2/gatk-protected/src/test/resources/large/mutect/dream_synthetic_bams/tumor_1.bam	/home/lichtens/test_onco_m2/gatk-protected/src/test/resources/large/mutect/dream_synthetic_bams/tumor_1.bam.bai	 synthetic.challenge.set1.tumor  /home/lichtens/test_onco_m2/gatk-protected/src/test/resources/large/mutect/dream_synthetic_bams/normal_1.bam    /home/lichtens/test_onco_m2/gatk-protected/src/test/resources/large/mutect/dream_synthetic_bams/normal_1.bam.bai	synthetic.challenge.set1.normal
```
