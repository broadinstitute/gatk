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

Additionally, the gatk git repository has a script called gatk (in the root directory of the repo) that is used to invoke the gatk.  If running on the cloud this is in the gatk docker image and you don't have to do anything.  If running on SGE, you must copy this script to a directory that is in your $PATH.

The following resources:
* Three preprocessed Hapmap vcfs -- one each for the 5-plex, 10-plex and 20-plex mixtures.  These are produced by preprocess_hapmap.wdl but as long as the sample composition of the mixtures remains the same they do not need to be generated again.  That is, the proportions need not be the same, but the same 5, 10, and 20 Hapmap samples must be present.
* A reference .fasta file, along with accompanying .fasta.fai and .dict files.
* A gatk4 java .jar file.
* Three lists of .bam files -- one each for 5-plex, 10-plex and 20-plex replicates -- where each row has the format <bam_file.bam></TAB><bam_index.bai>
* A list of .bam files of the specificity validation's replicates, where each row has the format <replicate_i.bam></TAB><replicate_i.bai></TAB><replicate_j.bam></TAB><replicate_j.bai>, with one row for each *ordered* pair i, j eg (1,2), (1,3), (2,1), (2,3), (3,1), (3,2) if there are three replicates.
* An intervals file.
* A gnomAD vcf.
* A Mutect2 panel of normals .vcf corresponding to the intervals and sequencing protocol of the specificity replicates.

## Preparing the wdl input files
In the same directory as your wdl scripts, fill in a file called sensitivity.json as follows:

```
{
  "HapmapSensitivityAllPlexes.gatk_override": "[Path to a gatk jar file.  Omitting this line uses the gatk jar in the docker image.]",
  "HapmapSensitivityAllPlexes.gatk_docker": "[gatk docker image eg broadinstitute/gatk:4.beta.3 -- this is not used in SGE but you still have to fill it in.]",
  "HapmapSensitivityAllPlexes.intervals": "[path to intervals file]",
  "HapmapSensitivityAllPlexes.ref_fasta": "[path to reference .fasta file]",
  "HapmapSensitivityAllPlexes.ref_fai": "[path to reference .fasta.fai file]",
  "HapmapSensitivityAllPlexes.ref_dict": "[path to reference .dict file]",
  "HapmapSensitivityAllPlexes.five_plex_bam_list": "[path to 5-plex bams list]",
  "HapmapSensitivityAllPlexes.ten_plex_bam_list": "[path to 10-plex bams list]",
  "HapmapSensitivityAllPlexes.twenty_plex_bam_list": "[path to 20-plex bams list]",
  "HapmapSensitivityAllPlexes.max_depth": "The maximum depth to consider for sensitivity.  1000 is a reasonable default.",
  "HapmapSensitivityAllPlexes.depth_bins": "Discrete depths at which to bin statistics.  [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800] is reasonable for many exomes",
  "HapmapSensitivityAllPlexes.depth_bin_width": "The width of depth bins.  Half the spacing betweens depths is reasonable.",
  "HapmapSensitivityAllPlexes.scatter_count": "[How many ways to scatter runs on Mutect2 on each bam file]",
  "HapmapSensitivityAllPlexes.run_orientation_bias_filter": "true/false depending on whether you wish to run this filter",
  "HapmapSensitivityAllPlexes.artifact_modes": The artifact modes of the orientation bias filter eg: ["G/T", "C/T"],
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
  "Mutect2ReplicateValidation.gatk_override": "[Path to a gatk jar file.  Omitting this line uses the gatk jar in the docker image.]",
  "Mutect2ReplicateValidation.gatk_docker": "[gatk docker image eg broadinstitute/gatk:4.beta.3 -- this is not used in SGE but you still have to fill it in.]",
  "Mutect2ReplicateValidation.ref_fasta": "[path to reference .fasta file]",
  "Mutect2ReplicateValidation.ref_fai": "[path to reference .fasta.fai file]",
  "Mutect2ReplicateValidation.ref_dict": "[path to reference .dict file]",
  "Mutect2ReplicateValidation.replicate_pair_list": "[path to replicate bams list]",
  "Mutect2ReplicateValidation.intervals": "[path to intervals file]",
  "Mutect2ReplicateValidation.pon": "[path to panel of normals vcf]",
  "Mutect2ReplicateValidation.pon_index": "[path to panel of normals vcf index]",
  "Mutect2ReplicateValidation.gnomad": "[path to panel of gnomAD vcf]",
  "Mutect2ReplicateValidation.gnomad_index": "[path to panel of gnomAD vcf index]",
  "Mutect2ReplicateValidation.scatter_count": "[How many ways to scatter runs on Mutect2 on each bam file]",
  "Mutect2ReplicateValidation.run_orientation_bias_filter": "true/false depending on whether you wish to run this filter",
  "Mutect2ReplicateValidation.artifact_modes": The artifact modes of the orientation bias filter eg: ["G/T", "C/T"],
  "Mutect2ReplicateValidation.preemptible_attempts": "2",
  "Mutect2ReplicateValidation.m2_extra_args": "optionally, any additional Mutect2 command line arguments",
  "Mutect2ReplicateValidation.m2_extra_filtering_args": "optionally, any additional Mutect2 command line arguments"  
}
```

Note that the docker image path is not used when the validations are run on an SGE cluster.  When running on SGE, a valid docker path must still be given or else cromwell will fail.

To summarize the differences between running in the cloud and on SGE:
* Your jsons must include a valid gatk_docker in both cases, however, when running on SGE this docker image is not actually used.
* When running in SGE you must put a gatk_override jar file in your jsons.  When running in the cloud you may include one but if you omit this line from your jsons the gatk jar in the docker image will be used.
* When running in SGE you must make sure to copy the gatk script in the root directory of the gatk git repo into a folder that is in your bash $PATH variable.

## Running in Cromwell
* Run hapmap_sensitivity_all_plexes.wdl with the parameters in sensitivity.json
* Run mutect2-replicate-validation.wdl with the parameters in specificity.json

## Outputs
The sensitivity validation outputs include many vcfs of true positives and false negatives used for debugging and improving Mutect.  The relevant outputs for validation results are:
* {snp, indel}\_{table, plot}\_{5, 10, 20, all}\_plex: tables in tsv format and graphs in png format of sensitivity versus depth and allele fraction for snvs and indels at each plex and aggregated over all plexes.
* {snp, indel}\_jaccard\_{5, 10, 20}\_plex: matrices in tsv format of the snv and indel jaccard index between each pair of replicates for each plex.  The jaccard index is the overlap of callsets divided by the union.
* summary\_{5, 10, 20}\_plex: the overall sensitivity (not binned by depth and allele fraction) for snvs and indels for each replicate of each plex.

The specificty validation's primary output, summary.txt, is a tsv file containing the rates of snv and indel false positives for each replicate pair.
