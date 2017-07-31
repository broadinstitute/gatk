### Mutect2 autovalidation

## Introduction
The Mutect2 autovalidation comprises a sensitivity validation and a specificity validation.

In the sensitivity validation, we mix (in vitro, not in silico) several Hapmap samples in roughly equal proportions to simulate a tumor with varying allele fractions, sequence the resulting mixture, and run Mutect2 in tumor-only mode.  Sensitivity to "somatic" variations is then defined as sensitivity to the known germline variants of the constituent samples.  A validation consists of 5-plex, 10-plex, and 20-plex mixtures, each with several replicates.

In the specificity validation, we make several replicates of a non-tumor sample and for each pair of these replicates we run Mutect2 in tumor-normal mode, with one replicate arbitrarily assigned as the "tumor."  Since every call is by definition a false positive, this yields a measure of specificity.


## Requirements

The following files from a clone of the gatk git repository, copied into a single directory:
* scripts/mutect2_wdl/mutect2.wdl
* scripts/mutect2_wdl/unsupported/hapmap_sensitivity.wdl
* scripts/mutect2_wdl/unsupported/hapmap_sensitivity_all_plexes.wdl
* scripts/mutect2_wdl/unsupported/mutect2-replicate-validation.wdl
* scripts/mutect2_wdl/unsupported/calculate_sensitivity.py

The following resources:
* Three preprocessed Hapmap vcfs -- one each for the 5-plex, 10-plex and 20-plex mixtures.  These are produced by preprocess_hapmap.wdl but as long as the sample composition of the mixtures remains the same they do not need to be generated again.  That is, the proportions need not be the same, but the same 5, 10, and 20 Hapmap samples must be present.
* A reference .fasta file, along with accompanying .fasta.fai and .dict files.
* A gatk4 java .jar file.
* A Picard java .jar file.
* Three lists of .bam files -- one each for 5-plex, 10-plex and 20-plex replicates -- where each row has the format <bam_file.bam></TAB><bam_index.bai></TAB><sample_name>
* A list of .bam files of the specificity validation's replicates.
* An intervals file.
* A gnomAD vcf.
* A Mutect2 panel of normals .vcf corresponding to the intervals and sequencing protocol of the specificity replicates.

## Preparing the wdl input files
In the same directory as your wdl scripts, fill in a file called sensitivity.json as follows:

```
{

  "HapmapSensitivityAllPlexes.max_depth": "The maximum depth to consider for sensitivity.  1000 is a reasonable default.",
  "HapmapSensitivityAllPlexes.gatk": "[path to gatk .jar file]",
  "HapmapSensitivityAllPlexes.ref_fasta": "[path to reference .fasta file]",
  "HapmapSensitivityAllPlexes.ref_fasta_index": "[path to reference .fasta.fai file]",
  "HapmapSensitivityAllPlexes.ref_dict": "[path to reference .dict file]",
  "HapmapSensitivityAllPlexes.five_plex_bam_list": "[path to 5-plex bams list]",
  "HapmapSensitivityAllPlexes.ten_plex_bam_list": "[path to 10-plex bams list]",
  "HapmapSensitivityAllPlexes.twenty_plex_bam_list": "[path to 20-plex bams list]",
  "HapmapSensitivityAllPlexes.intervals": "[path to intervals file]",
  "HapmapSensitivityAllPlexes.scatter_count": "[How many ways to scatter runs on Mutect2 on each bam file]",
  "HapmapSensitivityAllPlexes.is_run_orientation_bias_filter": "true/false depending on whether you wish to run this filter",
  "HapmapSensitivityAllPlexes.artifact_modes": The artifact modes of the orientation bias filter eg: ["G/T", "C/T"],
  "HapmapSensitivityAllPlexes.picard_jar": "[path to Picard .jar file]",
  "HapmapSensitivityAllPlexes.five_plex_preprocessed": "[path to preprocessed 5-plex vcf]",
  "HapmapSensitivityAllPlexes.five_plex_preprocessed_idx": "[path to preprocessed 5-plex vcf index]",
  "HapmapSensitivityAllPlexes.ten_plex_preprocessed": "[path to preprocessed 10-plex vcf]",
  "HapmapSensitivityAllPlexes.ten_plex_preprocessed_idx": "[path to preprocessed 10-plex vcf index]",
  "HapmapSensitivityAllPlexes.twenty_plex_preprocessed": "[path to preprocessed 20-plex vcf]",
  "HapmapSensitivityAllPlexes.twenty_plex_preprocessed_idx": "[path to preprocessed 20-plex vcf index]",
  "HapmapSensitivityAllPlexes.python_script": "path to calculate_sensitivity.py",
  "HapmapSensitivityAllPlexes.m2_extra_args": "optionally, any additional Mutect2 command line arguments",
  "HapmapSensitivityAllPlexes.m2_extra_filtering_args": "--maxEventsInHaplotype 100 --max_germline_posterior 1.0"
}
```

Note the extra filtering arguments hard-coded into these inputs.  These are necessary to disable filtering of germline variants, because the "somatic" variants here are actually germline variants.

In the same directory as your wdl scripts, fill in a file called specificity.json as follows:

```
{
  "Mutect2ReplicateValidation.gatk4_jar": "[path to gatk .jar file in the docker image if running on the cloud eg /root/gatk.jar]",
  "Mutect2ReplicateValidation.gatk4_jar_override": "[path to local gatk .jar file when not running in the cloud]",
  "Mutect2ReplicateValidation.ref_fasta": "[path to reference .fasta file]",
  "Mutect2ReplicateValidation.ref_fasta_index": "[path to reference .fasta.fai file]",
  "Mutect2ReplicateValidation.ref_dict": "[path to reference .dict file]",
  "Mutect2ReplicateValidation.replicate_pair_list": "[path to replicate bams list]",
  "Mutect2ReplicateValidation.intervals": "[path to intervals file]",
  "Mutect2ReplicateValidation.pon": "[path to panel of normals vcf]",
  "Mutect2ReplicateValidation.pon_index": "[path to panel of normals vcf index]",
  "Mutect2ReplicateValidation.gnomad": "[path to panel of gnomAD vcf]",
  "Mutect2ReplicateValidation.gnomad_index": "[path to panel of gnomAD vcf index]",
  "Mutect2ReplicateValidation.scatter_count": "[How many ways to scatter runs on Mutect2 on each bam file]",
  "Mutect2ReplicateValidation.is_run_orientation_bias_filter": "true/false depending on whether you wish to run this filter",
  "Mutect2ReplicateValidation.artifact_modes": The artifact modes of the orientation bias filter eg: ["G/T", "C/T"],
  "Mutect2ReplicateValidation.preemptible_attempts": "2",
  "Mutect2ReplicateValidation.m2_docker": "[gatk docker image eg broadinstitute/gatk:4.beta.3]",
  "Mutect2ReplicateValidation.picard_jar": "[path to Picard .jar file]",
  "Mutect2ReplicateValidation.m2_extra_args": "optionally, any additional Mutect2 command line arguments",
  "Mutect2ReplicateValidation.m2_extra_filtering_args": "optionally, any additional Mutect2 command line arguments"  
}
```

Note that the docker image path is not used when the validations are run locally.

## Running in Cromwell
* Run hapmap_sensitivity_all_plexes.wdl with the parameters in sensitivity.json
* Run mutect2-replicate-validation.wdl with the parameters in specificity.json

## Outputs
The sensitivity validation outputs include many vcfs of true positives and false negatives used for debugging and improving Mutect.  The relevant outputs for validation results are:
* {snp, indel}\_{table, plot}\_{5, 10, 20, all}\_plex: tables in tsv format and graphs in png format of sensitivity versus depth and allele fraction for snvs and indels at each plex and aggregated over all plexes.
* {snp, indel}\_jaccard\_{5, 10, 20}\_plex: matrices in tsv format of the snv and indel jaccard index between each pair of replicates for each plex.  The jaccard index is the overlap of callsets divided by the union.
* summary\_{5, 10, 20}\_plex: the overall sensitivity (not binned by depth and allele fraction) for snvs and indels for each replicate of each plex.

The specificty validation's primary output, summary.txt, is a tsv file containing the rates of snv and indel false positives for each replicate pair.