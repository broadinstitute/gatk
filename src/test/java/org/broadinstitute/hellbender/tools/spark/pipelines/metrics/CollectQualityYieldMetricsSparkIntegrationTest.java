package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public final class CollectQualityYieldMetricsSparkIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectQualityYieldMetrics");

    //Note: the 'expected' results in this test come from running picard 1.130
    //Note: we don't test the contents of the chart pdf
    //NOTE: these tests use the same data and results as the non-spark ones, by design

    @DataProvider(name = "CollectQualityYieldMetrics")
    private Iterator<Object[]> makeCollectQualityYieldMetricsData(){
        final List<Object[]> list= new ArrayList<>();
        list.add(new Object[]{"collect_quality_yield_metrics.sam", "collect_quality_yield_metrics.originalquals.txt", true});
        list.add(new Object[]{"collect_quality_yield_metrics.sam", "collect_quality_yield_metrics.quals.txt", false});
        return list.iterator();
    }
    @Test(dataProvider = "CollectQualityYieldMetrics")
    public void test(final String inName, final String outName, final boolean useOQ) throws IOException {
        final File input = new File(TEST_DATA_DIR, inName);
        final File expectedFile = new File(TEST_DATA_DIR, outName);   //file created using picard 1.130
        final File outfile = BaseTest.createTempFile("testCollectQualityYield", ".metrics");
        final String[] args = new String[]{
                "--I", input.getAbsolutePath(),
                "--O", outfile.getAbsolutePath(),
                "-OQ", String.valueOf(useOQ),
        };
        runCommandLine(args);
        IntegrationTestSpec.assertEqualTextFiles(outfile, expectedFile, "#");
    }
}
