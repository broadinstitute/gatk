export output_dir="/Users/tsato/workspace/gatk/src/test/java/org/broadinstitute/hellbender/tools/walkers/consensus"
export consensusDir="/Users/tsato/workspace/gatk/src/test/java/org/broadinstitute/hellbender/tools/walkers/consensus";
export fgbioGroupByUMIScript=$consensusDir"/fgbio_group_reads_by_umi.sh";
export fgbioConsensusScript=$consensusDir"/fgbio_consensus.sh";

export longDeletionDir="/dsde/working/tsato/consensus/indel/NA12878_truth/NA12878-INT/feb-2020-long-indels";

export bam="$longDeletionDir"/SM-EYXL1.RG_merged.chr22_9bp_deletion.bam
export output_bam="$longDeletionDir"/SM-EYXL1.RG_merged.chr22_9bp_deletion.grouped.bam
"$fgbioGroupByUMIScript" "$bam" "$output_bam" "$longDeletionDir"