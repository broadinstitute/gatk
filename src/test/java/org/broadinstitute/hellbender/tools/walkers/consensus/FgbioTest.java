package org.broadinstitute.hellbender.tools.walkers.consensus;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.testng.annotations.Test;

public class FgbioTest extends CommandLineProgramTest {
    private final String hg19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
    private final String tp53TestDir = "/dsde/working/tsato/consensus/tp53/test/";
    private final String tp53_AGA_TGA = tp53TestDir + "Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.AGG-TGA.grouped.bam";
    private final String tp53_full_grouped = tp53TestDir + "Jonna_Grimsby_A04_denovo_bloodbiopsy_1pct_rep1.tp53.umi_grouped.bam";

    private final String consensusDir = "/Users/tsato/workspace/gatk/src/test/java/org/broadinstitute/hellbender/tools/walkers/consensus/";
    private final String fgbioGroupByUMIScript = consensusDir + "fgbio_group_reads_by_umi.sh";
    private final String fgbioConsensusScript = consensusDir + "fgbio_consensus.sh";
    // /Users/tsato/workspace/gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/sort.sh
    private final String sortScript = consensusDir + "sort.sh";
    private final String realignScript = consensusDir + "realign.sh";
    private final String mergeBamAlignmentScript = consensusDir + "merge_bam_alignment.sh";

    private final String homeDir = "/dsde/working/tsato/consensus/bams/synthetic-test/";
    private final String testDir = "/dsde/working/tsato/consensus/tp53/test/";

    // 2/26/20 - Wed: why do long reads disappear?
    @Test
    public void testNA12878LongDeletion(){
        final String longDeletionDir = "/dsde/working/tsato/consensus/indel/NA12878_truth/NA12878-INT/feb-2020-long-indels/";
        // final String longDeletionDir = "/Users/tsato/Desktop/now/liquid_biopsy_files/feb-2020-long-indels/";
        final String basename = "SM-EYXL1.RG_merged.chr22_9bp_deletion";
        final String bam = longDeletionDir + basename + ".bam";
        final String groupedBam = longDeletionDir + basename + ".grouped.bam";
        final String debug = "debug";
        runProcess(new ProcessController(), new String[]{ fgbioGroupByUMIScript, bam, groupedBam, longDeletionDir});
        runProcess(new ProcessController(), new String[]{ fgbioConsensusScript, groupedBam, basename, longDeletionDir });
        // This must be a mistake--right?
        // runProcess(new ProcessController(), new String[]{ realignScript, groupedBam, basename, longDeletionDir });
        final String consensusBam = longDeletionDir + basename + ".consensus.bam";
        runProcess(new ProcessController(), new String[]{ realignScript, consensusBam, "SM-EYXL1.RG_merged.chr22.deletion.consensus", longDeletionDir });
        // TODO: input to the file should be the output name
        int d = 3;
    }
}
