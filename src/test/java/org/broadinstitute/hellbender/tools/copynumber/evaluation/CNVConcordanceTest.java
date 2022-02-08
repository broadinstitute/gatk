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
        final String dir = "/Users/tsato/workspace/cfCNV/analysis/ms_call_segments/";
        final String evalSeg = dir + "MP65_5753_in_KP_61683_3pct_duplicate_marked_called_by_me.seg";
        final String truthSeg = dir + "373505753_duplicate_marked_called_by_me.seg";

        final String te = dir + "truth_annotated_with_eval.seg";

        // To be implemented
        final String et = dir + "eval_annotated_with_truth.seg";
        final String summary = dir + "summary.txt";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("eval", evalSeg)
                .add("truth", truthSeg)
                .add(CNVConcordance.TRUTH_ANNOTATED_WITH_EVAL_SHORT_NAME, te)
                .add(CNVConcordance.SUMMARY_ARG_NAME, summary);
        runCommandLine(args, CNVConcordance.class.getSimpleName());
    }

}