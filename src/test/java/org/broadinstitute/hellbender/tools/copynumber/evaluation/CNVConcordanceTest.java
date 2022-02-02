package org.broadinstitute.hellbender.tools.copynumber.evaluation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class CNVConcordanceTest extends CommandLineProgramTest {
    @Test
    public void test(){
        final String hg19 = "/Users/tsato/workspace/cfCNV/resources/Homo_sapiens_assembly19.fasta";
        final String home = "/Users/tsato/workspace/cfCNV/evaluation/tool/";
        // final String eval = home + "1227495080.called.seg";
        final String eval = home +  "MP65_5753_in_KP_61683_3pct_duplicate_marked.called.seg";
        // final String truth = home + "SM-74P4M.called.seg";
        final String truth = home + "MP65_truth_373505753_duplicate_marked.called.seg";

        final String output1 = home + "output1.seg";
        final String output2 = home + "output2.seg";
        final String output3 = home + "output3.tsv";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("eval", eval)
                .add("truth", truth)
                .add("output1", output1)
                .add("output2", output2)
                .add("output3", output3);
        runCommandLine(args, CNVConcordance.class.getSimpleName());
    }

}