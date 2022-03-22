package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class ClipReadsForRSEMIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final String home = "/Volumes/dsde_working/tsato/hydro.gen/Analysis/874_twist_RNA/rsem/";
        //final String bam = home + "SM-LLM7I_Kapa_SM-LQZYG_test_1.bam";
        //final String bam = home + "SM-KYN26_Kapa_SM-LQZYH_test.bam";
        //final String output = home + "SM-LLM7I_Kapa_SM-LQZYG_test_1_out.bam";

        // Case 3
        final String bam = home + "test/SM-HESUP_SSIV_SM-LQZYY_genome_test.bam";
        final String output = home + "test/SM-HESUP_SSIV_SM-LQZYY_genome_test_out.bam";
        // 2101:10221:8484 has 2S104M40S and 43S103M, which would be clipped as is.

        // final String bam = home + "test/SM-HESUP_SSIV_SM-LQZYY_transcriptome_test.bam";
        // final String output = home + "test/SM-HESUP_SSIV_SM-LQZYY_transcriptome_test_out.bam";

        // final String bam = home + "test/SM-LVFDN_15_Min_Low_Mid_transcriptome.duplicate_marked.bam";
        // final String output = home + "test/SM-LVFDN_15_Min_Low_Mid_transcriptome_out.bam";

        args.addInput(bam);
        args.addOutput(output);
        runCommandLine(args, ClipReadsForRSEM.class.getSimpleName());
    }
}