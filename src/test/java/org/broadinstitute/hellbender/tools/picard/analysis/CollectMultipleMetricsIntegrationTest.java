package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.TestSparkProgramGroup;
import org.broadinstitute.hellbender.metrics.InsertSizeMetrics;
import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public final class CollectMultipleMetricsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/analysis/CollectInsertSizeMetrics");

    @Override
    public String getTestedClassName() {
        return CollectMultipleMetrics.class.getSimpleName();
    }

    @DataProvider(name="metricsTestFiles")
    public Object[][] insertSizeMetricsFiles() {
        return new Object[][] {
                // single level collection
                {"insert_size_metrics_test.sam", null, "expectedInsertSizeMetricsL1.txt"},
                {"insert_size_metrics_test.bam", null, "expectedInsertSizeMetricsL1.txt"},
                {"insert_size_metrics_test.cram", hg19_chr1_1M_Reference, "expectedInsertSizeMetricsL1.txt"},
        };
    }

    @Test(dataProvider="metricsTestFiles")
    public void test(
            final String fileName,
            final String referenceName,
            final String expectedInsertSizeResults) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        String outBase = setupMultipleCollector(args, fileName, referenceName);

        // for now, run the only two conforming collectors that we have
        args.add("--PROGRAM");
        args.add("CollectInsertSizeMetrics");

        this.runCommandLine(args.getArgsArray());

        validateInsertSizeMetrics(outBase, expectedInsertSizeResults);
    }

    private String setupMultipleCollector(
            final ArgumentsBuilder args,
            final String fileName,
            final String referenceName) throws IOException
    {
        // create a directory to contain the results since there will be multiple collectors
        // and each may create multiple files
        final File outDir = BaseTest.createTempDir("collectMultiMetricsTest" );
        String outBase = outDir.getAbsolutePath() + "/collectMultiMetrics";

        // IO arguments
        args.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(new File(TEST_DATA_DIR, fileName).getAbsolutePath());

        // The output arg specifies only the basename from which output file(s)
        // will derived
        args.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(outBase);

        if (null != referenceName) {
            final File REF = new File(referenceName);
            args.add("-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME);
            args.add(REF.getAbsolutePath());
        }

        return outBase;
    }

    private void validateInsertSizeMetrics(final String outBase, final String expectedResults) throws IOException {
        String localOut = outBase + "." + InsertSizeMetrics.getUniqueNameSuffix() + ".txt";

        IntegrationTestSpec.assertEqualTextFiles(
                new File(localOut),
                new File(TEST_DATA_DIR, expectedResults),
                "#");
    }

    // Test implementation of MetricsCollectorSpark used for testing CollectMultipleMetricsSpark
    // with a custom collector added programmatically
    @CommandLineProgramProperties(programGroup= TestSparkProgramGroup.class,
            summary="test custom collector", oneLineSummary = "test custom collector")
    public static class TestCustomCollector extends SinglePassSamProgram {
        long count = 0;
        public void initialize(
                MetricsArgumentCollection inputArgs, SAMFileHeader samHeader, List<Header> defaultHeaders) {}
        @Override
        protected void setup(final SAMFileHeader header, final File samFile) {}
        @Override
        protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
            count++;
        };
        @Override
        public void finish() {}
    }

    @Test(dataProvider="metricsTestFiles")
    public void testCustomCollectorAPI(
            final String fileName,
            final String referenceName,
            final String expectedInsertSizeResults) throws IOException {
        {
            // Test CollectMultipleMetricsSpark with a custom collector
            final TestCustomCollector testCollector = new TestCustomCollector();

            ArgumentsBuilder args = new ArgumentsBuilder();
            setupMultipleCollector(args, fileName, referenceName);

            // CollectMultipleMetrics provider that creates an initializes a custom collector
            CollectMultipleMetrics.ProgramInterface customProvider = new CollectMultipleMetrics.ProgramInterface() {
                @Override
                public SinglePassSamProgram makeInstance(final String outbase) {
                    return testCollector;
                }
            };

            // Manually create a tool and programmatically set the custome collector as the one to run
            CollectMultipleMetrics multipleCollectorTool = new CollectMultipleMetrics();
            multipleCollectorTool.setProgramsToRun(Collections.singletonList(customProvider));
            multipleCollectorTool.instanceMain(args.getArgsArray());

            Assert.assertEquals(testCollector.count, 52L);
        }
    }

}
