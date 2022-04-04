package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

public class PostProcessReadsForRSEMIntegrationTest extends CommandLineProgramTest {
    @Test
    public void test(){
        final ArgumentsBuilder args = new ArgumentsBuilder();
        final String home = "/Volumes/dsde_working/tsato/hydro.gen/Analysis/874_twist_RNA/rsem/";

        // TEst
        // Case 4
        final String bam = home + "test/SM-LVFDN_15_Min_Low_Mid_transcriptome.duplicate_marked.bam";
        final String output = home + "test/SM-LVFDN_15_Min_Low_Mid_transcriptome_out.bam";

        // Case 5
        args.addInput(bam);
        args.addOutput(output);
        runCommandLine(args, PostProcessReadsForRSEM.class.getSimpleName());
    }
}