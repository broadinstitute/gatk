package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.tools.picard.sam.markduplicates.MarkDuplicatesIntegrationTest;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesSparkTester;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.testers.AbstractMarkDuplicatesCommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.testers.AbstractMarkDuplicatesTester;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


public class MarkDuplicatesSparkIntegrationTest extends AbstractMarkDuplicatesCommandLineProgramTest {

    @Override
    protected AbstractMarkDuplicatesTester getTester() {
        return new MarkDuplicatesSparkTester();
    }

    @Override
    protected CommandLineProgram getCommandLineProgramInstance() {
        return new MarkDuplicatesSpark();
    }

    @Test(dataProvider = "testMDdata")
    @Override
    public void testMDOrder(final File input, final File expectedOutput) throws Exception {
        // Override this test case to provide a --shardedOutput false argument, so that we write a single, sorted
        // bam (since sharded output is not sorted, and this test case is sensitive to order).
        testMDOrderImpl(input, expectedOutput, "--shardedOutput false");
    }

    @DataProvider(name = "md")
    public Object[][] md(){
        return new Object[][]{
            // The first two values are total reads and duplicate reads. The list is an encoding of the metrics
            // file output by this bam file. These metrics files all match the outputs of picard mark duplicates.

             //Note: in each of those cases, we'd really want to pass null as the last parameter (not 0L) but IntelliJ
             // does not like it and skips the test (rendering issue) - so we pass 0L and account for it at test time
             // (see comment in testMarkDuplicatesSparkIntegrationTestLocal)
            {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam"), 20, 0,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16406", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L))},
            {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.bam"), 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR,"example.chr1.1-1K.markedDups.bam"), 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR, "optical_dupes.bam"), 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
            {new File(MarkDuplicatesIntegrationTest.TEST_DATA_DIR, "optical_dupes_casava.bam"), 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
        };
    }

    @Test(groups = "spark", dataProvider = "md")
    public void testMarkDuplicatesSparkIntegrationTestLocal(
        final File input, final long totalExpected, final long dupsExpected,
        Map<String, List<String>> metricsExpected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+ StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(input.getPath());
        args.add("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.add(outputFile.getAbsolutePath());

        args.add("--METRICS_FILE");
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.add(metricsFile.getAbsolutePath());

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        int totalReads = 0;
        int duplicateReads = 0;
        try ( final ReadsDataSource outputReads = new ReadsDataSource(outputFile) ) {
            for ( GATKRead read : outputReads ) {
                ++totalReads;

                if ( read.isDuplicate() ) {
                    ++duplicateReads;
                }
            }
        }

        Assert.assertEquals(totalReads, totalExpected, "Wrong number of reads in output BAM");
        Assert.assertEquals(duplicateReads, dupsExpected, "Wrong number of duplicate reads in output BAM");

        final MetricsFile<DuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<DuplicationMetrics> nonEmptyMetrics = metricsOutput.getMetrics().stream().filter(
                metric ->
                    metric.UNPAIRED_READS_EXAMINED != 0L ||
                    metric.READ_PAIRS_EXAMINED != 0L ||
                    metric.UNMAPPED_READS != 0L ||
                    metric.UNPAIRED_READ_DUPLICATES != 0L ||
                    metric.READ_PAIR_DUPLICATES != 0L ||
                    metric.READ_PAIR_OPTICAL_DUPLICATES != 0L ||
                    (metric.PERCENT_DUPLICATION != null && metric.PERCENT_DUPLICATION != 0.0 && !Double.isNaN(metric.PERCENT_DUPLICATION)) ||
                    (metric.ESTIMATED_LIBRARY_SIZE != null && metric.ESTIMATED_LIBRARY_SIZE != 0L)
        ).collect(Collectors.toList());

        Assert.assertEquals(nonEmptyMetrics.size(), metricsExpected.size(),
                            "Wrong number of metrics with non-zero fields.");
        for (int i = 0; i < nonEmptyMetrics.size(); i++ ){
            final DuplicationMetrics observedMetrics = nonEmptyMetrics.get(i);
            List<?> expectedList = metricsExpected.get(observedMetrics.LIBRARY);
            Assert.assertNotNull(expectedList, "Unexpected library found: " + observedMetrics.LIBRARY);
            Assert.assertEquals(observedMetrics.UNPAIRED_READS_EXAMINED, expectedList.get(0));
            Assert.assertEquals(observedMetrics.READ_PAIRS_EXAMINED, expectedList.get(1));
            Assert.assertEquals(observedMetrics.UNMAPPED_READS, expectedList.get(2));
            Assert.assertEquals(observedMetrics.UNPAIRED_READ_DUPLICATES, expectedList.get(3));
            Assert.assertEquals(observedMetrics.READ_PAIR_DUPLICATES, expectedList.get(4));
            Assert.assertEquals(observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, expectedList.get(5));
            Assert.assertEquals(observedMetrics.PERCENT_DUPLICATION, expectedList.get(6));

            //Note: IntelliJ does not like it when a parameter for a test is null (can't print it and skips the test)
            //so we work around it by passing in an 'expected 0L' and only comparing to it if the actual value is non-null
            if (observedMetrics.ESTIMATED_LIBRARY_SIZE != null && (Long)expectedList.get(7) != 0L)  {
                Assert.assertEquals(observedMetrics.ESTIMATED_LIBRARY_SIZE, expectedList.get(7));
            }
        }
    }
}
