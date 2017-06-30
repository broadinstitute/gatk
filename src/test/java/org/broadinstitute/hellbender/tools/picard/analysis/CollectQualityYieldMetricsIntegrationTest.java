package org.broadinstitute.hellbender.tools.picard.analysis;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.TestResources;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public final class CollectQualityYieldMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectQualityYieldMetrics");

    //Note: the 'expected' results in this test come from running picard 1.130
    //Note: we don't test the contents of the chart pdf

    @DataProvider(name = "CollectQualityYieldMetrics")
    private Iterator<Object[]> makeCollectQualityYieldMetricsData(){
        final List<Object[]> list= new ArrayList<>();

        list.add(new Object[]{"valid.bam", "valid.CollectQualityYieldMetrics.txt", null, true});
        list.add(new Object[]{"valid.cram", "valid.CollectQualityYieldMetrics.txt", new File(TEST_DATA_DIR, "valid.fasta").getAbsolutePath(), true});

        list.add(new Object[]{"collect_quality_yield_metrics.sam", "collect_quality_yield_metrics.originalquals.txt", null, true});
        list.add(new Object[]{"collect_quality_yield_metrics.sam", "collect_quality_yield_metrics.quals.txt", null, false});
        list.add(new Object[]{"collect_quality_yield_metrics.bam", "collect_quality_yield_metrics.originalquals.txt", null, true});
        list.add(new Object[]{"collect_quality_yield_metrics.bam", "collect_quality_yield_metrics.quals.txt", null, false});
        list.add(new Object[]{"collect_quality_yield_metrics.cram", "collect_quality_yield_metrics.originalquals.txt", TestResources.hg19_chr1_1M_Reference, true});
        list.add(new Object[]{"collect_quality_yield_metrics.cram", "collect_quality_yield_metrics.quals.txt", TestResources.hg19_chr1_1M_Reference, false});
        return list.iterator();
    }
    @Test(dataProvider = "CollectQualityYieldMetrics")
    public void test(final String inName, final String outName, final String referenceName, final boolean useOQ) throws IOException {
        final File input = new File(TEST_DATA_DIR, inName);
        final File expectedFile = new File(TEST_DATA_DIR, outName);   //file created using picard 1.130
        final File outfile = BaseTest.createTempFile("testCollectQualityYield", ".metrics");

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--I");
        args.add(input.getCanonicalPath());
        args.add("--O");
        args.add(outfile.getCanonicalPath());
        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.add("--R");
            args.add(REF.getAbsolutePath());
        }
        args.add("-OQ");
        args.add(String.valueOf(useOQ));

        runCommandLine(args.getArgsList());
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
