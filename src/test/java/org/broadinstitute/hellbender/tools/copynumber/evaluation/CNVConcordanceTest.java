package org.broadinstitute.hellbender.tools.copynumber.evaluation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class CNVConcordanceTest extends CommandLineProgramTest {
    final String hg19 = "/Users/tsato/workspace/cfCNV/resources/Homo_sapiens_assembly19.fasta";

    @Test
    public void test(){
        final String dir = "/Users/tsato/workspace/cfCNV/analysis/ichor_output/";
        final String evalSeg = dir + "RP-1544_Pt107_M7_with_SM-LJGJ7_KP_61893.seg.txt";
        final String truthSeg = dir + "RP-1544_Pt107_M7_v1_WGS_OnPrem.seg.txt";

        final String te = dir + "truth_annotated_with_eval.seg";
        final String summary = dir + "summary.txt";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("eval", evalSeg)
                .add("truth", truthSeg)
                .add(CNVConcordance.FILE_TYPE_NAME, "ICHORCNA")
                .add(CNVConcordance.TRUTH_ANNOTATED_WITH_EVAL_SHORT_NAME, te)
                .add(CNVConcordance.SUMMARY_ARG_NAME, summary);
        runCommandLine(args, CNVConcordance.class.getSimpleName());
    }

    @Test
    public void testMS(){
        final String dir = "/Users/tsato/workspace/cfCNV/analysis/ms_call_segments/";
        final String evalSeg = dir + "MP65_5753_in_KP_61683_3pct_duplicate_marked_called_by_me.seg";
        final String truthSeg = dir + "373505753_duplicate_marked_called_by_me.seg";

        final String te = dir + "truth_annotated_with_eval.seg";
        final String summary = dir + "summary.txt";

        // To be implemented
        final String et = dir + "eval_annotated_with_truth.seg";


        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("R", hg19)
                .add("eval", evalSeg)
                .add("truth", truthSeg)
                .add(CNVConcordance.FILE_TYPE_NAME, "MODELSEGMENTS")
                .add(CNVConcordance.TRUTH_ANNOTATED_WITH_EVAL_SHORT_NAME, te)
                .add(CNVConcordance.SUMMARY_ARG_NAME, summary);
        runCommandLine(args, CNVConcordance.class.getSimpleName());
    }

}