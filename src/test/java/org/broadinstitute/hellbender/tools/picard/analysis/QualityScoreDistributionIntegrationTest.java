package org.broadinstitute.hellbender.tools.picard.analysis;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class QualityScoreDistributionIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/QualityScoreDistribution");

    //Note: the 'expected' results in this test come from running picard 1.130

    //Note: we don't test the contents of the chart pdf

    @DataProvider(name = "QualityScoreDistribution")
    private Iterator<Object[]> makeQualityScoreDistributionData(){
        final List<Object[]> list= new ArrayList<>();
        list.add(new Object[]{"first5000a.bam", "qualscoredist.txt", null, true, false, false});
        list.add(new Object[]{"first5000a.cram", "qualscoredist.txt", b37_reference_20_21, true, false, false});
        list.add(new Object[]{"originalQuals.chr1.1-1K.bam", "originalQuals.chr1.1-1K.QualityScoreDistribution.txt", null, true, false, false});

        list.add(new Object[]{"example_pfFail_reads.bam", "pfFailBam.pf.txt", null, true, false, false});
        list.add(new Object[]{"example_pfFail_reads.bam", "pfFailBam.pfOnly.txt", null, true, true, false});

        list.add(new Object[]{"unmapped.bam", "unmappedBam.ALIGNED_READS_ONLY_false.txt", null, true, false, false});
        list.add(new Object[]{"unmapped.bam", "unmappedBam.ALIGNED_READS_ONLY_true.txt", null, false, false, true});

        return list.iterator();
    }

    @Test(dataProvider = "QualityScoreDistribution", groups = {"R"})
    public void test(final String unsortedBamName, final String expectedFileName, final String referenceName,
                     final boolean makePdf, final boolean pfReadsOnly, final boolean alignedReadsOnly) throws IOException {
        final File unsortedBam = new File(TEST_DATA_DIR, unsortedBamName);
        final File expectedFile = new File(TEST_DATA_DIR, expectedFileName);

        //Note we compare to non-spark outputs
        final File outfile = BaseTest.createTempFile("test", ".metrics");
        final File pdf = BaseTest.createTempFile("test", ".pdf");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--" + "input");
        args.add(unsortedBam.getCanonicalPath());
        args.add("--" + "output");
        args.add(outfile.getCanonicalPath());
        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.add("-R");
            args.add(REF.getAbsolutePath());
        }
        args.add("--" + "CHART_OUTPUT");
        args.add(pdf.getCanonicalPath());
        args.add("--" + "PRODUCE_PLOT");
        args.add(makePdf);

        args.add("--" + "PF_READS_ONLY");
        args.add(pfReadsOnly);
        args.add("--" + "ALIGNED_READS_ONLY");
        args.add(alignedReadsOnly);

        this.runCommandLine(args.getArgsArray());

        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
