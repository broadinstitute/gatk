package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.tools.walkers.markduplicates.AbstractMarkDuplicatesCommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.testers.MarkDuplicatesSparkTester;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.GATKDuplicationMetrics;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.MarkDuplicates;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

@Test(groups = "spark")
public class MarkDuplicatesSparkIntegrationTest extends AbstractMarkDuplicatesCommandLineProgramTest {

    @Override
    protected MarkDuplicatesSparkTester getTester() {
        MarkDuplicatesSparkTester markDuplicatesSparkTester = new MarkDuplicatesSparkTester();
        markDuplicatesSparkTester.addArg("--"+ MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME);
        return markDuplicatesSparkTester;
    }

    @Override
    protected CommandLineProgram getCommandLineProgramInstance() {
        return new MarkDuplicatesSpark();
    }

    @Override
    protected boolean markSecondaryAndSupplementaryRecordsLikeTheCanonical() { return true; }

    @Test(dataProvider = "testMDdata", groups = "spark")
    @Override
    public void testMDOrder(final File input, final File expectedOutput) throws Exception {
        // Override this test case to provide a --sharded-output false argument, so that we write a single, sorted
        // bam (since sharded output is not sorted, and this test case is sensitive to order).
        testMDOrderImpl(input, expectedOutput, "--" + GATKSparkTool.SHARDED_OUTPUT_LONG_NAME +" false");
    }

    @DataProvider(name = "md")
    public Object[][] md(){
        return new Object[][]{
            // The first two values are total reads and duplicate reads. The list is an encoding of the metrics
            // file output by this bam file. These metrics files all match the outputs of picard mark duplicates.

             //Note: in each of those cases, we'd really want to pass null as the last parameter (not 0L) but IntelliJ
             // does not like it and skips the test (rendering issue) - so we pass 0L and account for it at test time
             // (see comment in testMarkDuplicatesSparkIntegrationTestLocal)
            {new File[]{new File(TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam")}, 20, 0,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16406", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L))},
            {new File[]{new File(TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.bam")}, 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File[]{new File(TEST_DATA_DIR,"example.chr1.1-1K.markedDups.bam")}, 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File[]{new File(TEST_DATA_DIR, "optical_dupes.bam")}, 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
            {new File[]{new File(TEST_DATA_DIR, "optical_dupes_casava.bam")}, 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
            {new File[]{new File(TEST_DATA_DIR, "optical_dupes.queryname.bam"), new File(TEST_DATA_DIR, "example.chr1.1-1K.markedDups.queryname.bam")}, 94, 8,
             ImmutableMap.builder().put("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))
                                    .put("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L))
                                    .put("Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L))
                                    .put("Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L))
                                    .put("Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L))
                                    .put("Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L)).build()},
            // Testing that that multiple input file behavior still functions when both files are mixed between queryname and querygroup sorting
            {new File[]{new File(TEST_DATA_DIR, "optical_dupes.queryname.bam"), new File(TEST_DATA_DIR, "example.chr1.1-1K.markedDups.querygrouped.bam")}, 94, 8,
                ImmutableMap.builder().put("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))
                        .put("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L))
                        .put("Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L))
                        .put("Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L))
                        .put("Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L))
                        .put("Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L)).build()},
        };
    }

    // Tests asserting that without --do-not-mark-unmapped-mates argument that unmapped mates are still duplicate marked with their partner
    @Test
    public void testMappedPairAndMappedFragmentAndMatePairSecondUnmapped() {
        final MarkDuplicatesSparkTester tester = new MarkDuplicatesSparkTester(true);
        tester.addMatePair(1, 10040, 10040, false, true, true, true, "76M", null, false, false, false, false, false, DEFAULT_BASE_QUALITY); // first a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, ELIGIBLE_BASE_QUALITY); // mapped OK
        tester.addMappedFragment(1, 10040, true, DEFAULT_BASE_QUALITY); // duplicate
        tester.runTest();
    }

    @Test( dataProvider = "md")
    public void testMarkDuplicatesSparkIntegrationTestLocal(
        final File[] inputFiles, final long totalExpected, final long dupsExpected,
        Map<String, List<String>> metricsExpected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        for (File input : inputFiles) {
            args.add(StandardArgumentDefinitions.INPUT_LONG_NAME,input.getPath());
        }
        args.addFlag(StandardArgumentDefinitions.OUTPUT_LONG_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.addRaw(outputFile.getAbsolutePath());

        args.addRaw("--"+StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.addRaw(metricsFile.getAbsolutePath());

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        int totalReads = 0;
        int duplicateReads = 0;
        try ( final ReadsDataSource outputReads = new ReadsPathDataSource(outputFile.toPath()) ) {
            for ( GATKRead read : outputReads ) {
                ++totalReads;

                if ( read.isDuplicate() ) {
                    ++duplicateReads;
                }
            }
        }

        Assert.assertEquals(totalReads, totalExpected, "Wrong number of reads in output BAM");
        Assert.assertEquals(duplicateReads, dupsExpected, "Wrong number of duplicate reads in output BAM");

        final MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<GATKDuplicationMetrics> nonEmptyMetrics = getGatkDuplicationMetrics(metricsOutput);

        Assert.assertEquals(nonEmptyMetrics.size(), metricsExpected.size(),
                            "Wrong number of metrics with non-zero fields.");
        compareMetricsToExpected(metricsExpected, nonEmptyMetrics);
    }

    protected void compareMetricsToExpected(Map<String, List<String>> metricsExpected, List<GATKDuplicationMetrics> nonEmptyMetrics) {
        for (int i = 0; i < nonEmptyMetrics.size(); i++) {
            final GATKDuplicationMetrics observedMetrics = nonEmptyMetrics.get(i);
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
            if (observedMetrics.ESTIMATED_LIBRARY_SIZE != null && (Long) expectedList.get(7) != 0L) {
                Assert.assertEquals(observedMetrics.ESTIMATED_LIBRARY_SIZE, expectedList.get(7));
            }
        }
    }


    protected List<GATKDuplicationMetrics> getGatkDuplicationMetrics(MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput) {
        return metricsOutput.getMetrics().stream().filter(
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
    }

    @Test( dataProvider = "md")
    public void testMarkDuplicatesSparkMarkingOpticalDuplicatesWithTagging(
            final File[] inputFiles, final long totalExpected, final long dupsExpected,
            Map<String, List<String>> metricsExpected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        for (File input : inputFiles) {
            args.add(StandardArgumentDefinitions.INPUT_LONG_NAME,input.getPath());
        }
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.addRaw(outputFile.getAbsolutePath());

        args.addRaw("--" + StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.addRaw(metricsFile.getAbsolutePath());

        args.addRaw("--" + MarkDuplicatesSparkArgumentCollection.DUPLICATE_TAGGING_POLICY_LONG_NAME);
        args.addRaw(MarkDuplicates.DuplicateTaggingPolicy.OpticalOnly);

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        Set<String> opticalDuplicateNamesObserved = new HashSet<>();

        final MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<GATKDuplicationMetrics> nonEmptyMetrics = getGatkDuplicationMetrics(metricsOutput);


        int totalReads = 0;
        int duplicateReads = 0;
        try (final ReadsDataSource outputReads = new ReadsPathDataSource(outputFile.toPath())) {
            for (GATKRead read : outputReads) {
                ++totalReads;

                if (read.isDuplicate()) {
                    ++duplicateReads;
                }

                if (MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(read.getAttributeAsString(MarkDuplicates.DUPLICATE_TYPE_TAG))) {
                    opticalDuplicateNamesObserved.add(read.getName());
                }
            }
        }

        // Count the number of optical duplicates in the metrics and then compare against duplicates observed in the file
        int expectedOpticalDuplicatesGroups = 0;
        for (int i = 0; i < nonEmptyMetrics.size(); i++ ){
            final GATKDuplicationMetrics observedMetrics = nonEmptyMetrics.get(i);
            List<?> expectedList = metricsExpected.get(observedMetrics.LIBRARY);
            expectedOpticalDuplicatesGroups += (Long) (expectedList.get(5));
            Assert.assertEquals(observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, expectedList.get(5));
        }
        Assert.assertEquals(opticalDuplicateNamesObserved.size(), expectedOpticalDuplicatesGroups, "Wrong number of optical duplicate marked reads in output BAM");
    }

    @Test( dataProvider = "md")
    // Testing the DUPLICATE_TAGGING_POLICY_LONG_NAME = ALL option.
    public void testMarkDuplicatesSparkMarkingAllDuplicatesWithTagging(
            final File[] inputFiles, final long totalExpected, final long dupsExpected,
            Map<String, List<String>> metricsExpected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        for (File input : inputFiles) {
            args.add(StandardArgumentDefinitions.INPUT_LONG_NAME,input.getPath());
        }
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.addRaw(outputFile.getAbsolutePath());

        args.addRaw("--" + StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.addRaw(metricsFile.getAbsolutePath());

        args.addRaw("--" + MarkDuplicatesSparkArgumentCollection.DUPLICATE_TAGGING_POLICY_LONG_NAME);
        args.addRaw(MarkDuplicates.DuplicateTaggingPolicy.All);

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        final Set<String> opticalDuplicateNamesObserved = new HashSet<>();
        final Set<String> libraryDuplicateNamesObserved = new HashSet<>();

        final MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<GATKDuplicationMetrics> nonEmptyMetrics = getGatkDuplicationMetrics(metricsOutput);


        int totalReads = 0;
        int duplicateReads = 0;
        try (final ReadsDataSource outputReads = new ReadsPathDataSource(outputFile.toPath())) {
            for (GATKRead read : outputReads) {
                ++totalReads;

                if (read.isDuplicate()) {
                    ++duplicateReads;
                }

                if (MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(read.getAttributeAsString(MarkDuplicates.DUPLICATE_TYPE_TAG))) {
                    opticalDuplicateNamesObserved.add(read.getName());
                }
                if (MarkDuplicates.DUPLICATE_TYPE_LIBRARY.equals(read.getAttributeAsString(MarkDuplicates.DUPLICATE_TYPE_TAG))) {
                    libraryDuplicateNamesObserved.add(read.getName());
                }
            }
        }

        // Count the number of optical and library duplicates in the metrics and then compare against duplicates observed in the file
        int expectedOpticalDuplicatesGroups = 0;
        int expectedLibraryDuplicatesGroups = 0;
        for (int i = 0; i < nonEmptyMetrics.size(); i++ ){
            final GATKDuplicationMetrics observedMetrics = nonEmptyMetrics.get(i);
            List<?> expectedList = metricsExpected.get(observedMetrics.LIBRARY);
            expectedLibraryDuplicatesGroups += (Long) (expectedList.get(4));
            expectedOpticalDuplicatesGroups += (Long) (expectedList.get(5));
            Assert.assertEquals(observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, expectedList.get(5));
        }
        Assert.assertEquals(opticalDuplicateNamesObserved.size(), expectedOpticalDuplicatesGroups, "Wrong number of optical duplicate marked reads in output BAM");
        Assert.assertEquals(libraryDuplicateNamesObserved.size(), expectedLibraryDuplicatesGroups - expectedOpticalDuplicatesGroups, "Wrong number of duplicate marked reads in output BAM");
    }


    @Test( dataProvider = "md")
    public void testMarkDuplicatesSparkDeletingDuplicateReads(
            final File[] inputFiles, final long totalExpected, final long dupsExpected,
            Map<String, List<String>> metricsExpected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        for (File input : inputFiles) {
            args.add(StandardArgumentDefinitions.INPUT_LONG_NAME,input.getPath());
        }
        args.addRaw("--"+StandardArgumentDefinitions.OUTPUT_LONG_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.addRaw(outputFile.getAbsolutePath());

        args.addRaw("--"+StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.addRaw(metricsFile.getAbsolutePath());

        args.addRaw("--"+MarkDuplicatesSparkArgumentCollection.REMOVE_ALL_DUPLICATE_READS);

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        int totalReads = 0;
        int duplicateReads = 0;
        try ( final ReadsDataSource outputReads = new ReadsPathDataSource(outputFile.toPath()) ) {
            for ( GATKRead read : outputReads ) {
                ++totalReads;

                if ( read.isDuplicate() ) {
                    ++duplicateReads;
                }
            }
        }

        Assert.assertEquals(totalReads, totalExpected - dupsExpected, "Wrong number of reads in output BAM");
        Assert.assertEquals(duplicateReads, 0, "Wrong number of duplicate reads in output BAM");

        final MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<GATKDuplicationMetrics> nonEmptyMetrics = getGatkDuplicationMetrics(metricsOutput);

        // Assert that the metrics haven't changed at all
        Assert.assertEquals(nonEmptyMetrics.size(), metricsExpected.size(),
                "Wrong number of metrics with non-zero fields.");
        compareMetricsToExpected(metricsExpected, nonEmptyMetrics);
    }


    @Test( dataProvider = "md")
    public void testMarkDuplicatesSparkDeletingOpticalDuplicateReads(
            final File[] inputFiles, final long totalExpected, final long dupsExpected,
            Map<String, List<String>> metricsExpected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        for (File input : inputFiles) {
            args.add(StandardArgumentDefinitions.INPUT_LONG_NAME,input.getPath());
        }
        args.addRaw("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.addRaw(outputFile.getAbsolutePath());

        args.addRaw("--" + StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.addRaw(metricsFile.getAbsolutePath());

        args.addRaw("--" + MarkDuplicatesSparkArgumentCollection.REMOVE_SEQUENCING_DUPLICATE_READS);

        runCommandLine(args.getArgsArray());

        Assert.assertTrue(outputFile.exists(), "Can't find expected MarkDuplicates output file at " + outputFile.getAbsolutePath());

        int totalReads = 0;
        int duplicateReads = 0;
        try (final ReadsDataSource outputReads = new ReadsPathDataSource(outputFile.toPath())) {
            for (GATKRead read : outputReads) {
                ++totalReads;

                if (read.isDuplicate()) {
                    ++duplicateReads;
                }
            }
        }

        final MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<GATKDuplicationMetrics> nonEmptyMetrics = getGatkDuplicationMetrics(metricsOutput);

        // Asserting
        int expectedOpticalDuplicatesGroups = 0;
        for (int i = 0; i < nonEmptyMetrics.size(); i++ ){
            final GATKDuplicationMetrics observedMetrics = nonEmptyMetrics.get(i);
            List<?> expectedList = metricsExpected.get(observedMetrics.LIBRARY);
            expectedOpticalDuplicatesGroups += (Long) expectedList.get(5);
        }
        // NOTE: this test will fail if we add a more comprehensive example set with optical duplicates containing secondary/supplementary reads because of how the test is counting optical duplicate reads.
        Assert.assertEquals(totalReads, totalExpected - expectedOpticalDuplicatesGroups*2, "Wrong number of reads in output BAM");
        Assert.assertEquals(duplicateReads, dupsExpected - expectedOpticalDuplicatesGroups*2, "Wrong number of duplicate reads in output BAM");
    }


    @Test
    public void testHashCollisionHandling() {
        // This test asserts that the handling of two read pairs with the same start positions but on different in such a way
        // that they might cause hash collisions are handled properly.
        final File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(getTestFile("hashCollisionedReads.bam"));
        runCommandLine(args);

        try ( final ReadsDataSource outputReadsSource = new ReadsPathDataSource(output.toPath()) ) {
            final List<GATKRead> actualReads = new ArrayList<>();
            for ( final GATKRead read : outputReadsSource ) {
                Assert.assertFalse(read.isDuplicate());
                actualReads.add(read);
            }

            Assert.assertEquals(actualReads.size(), 4, "Wrong number of reads output");
        }
    }

    @Test (expectedExceptions = UserException.class)
    public void testAssertCorrectSortOrderMultipleBams() {
        // This test asserts that two bams with different sort orders will throw an exception for MarkDuplicatesSpark if both
        // are supplied as inputs to the tool (currently we require all bams in multi-inputs to be querygroup/queryname sorted).
        final File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(new File(TEST_DATA_DIR,"optical_dupes.bam"));
        args.addInput(new File(TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam"));
        runCommandLine(args);
    }

    @Test
    public void testAssertCorrectSortOrderMultipleBamsOverriding() {
        final File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(new File(TEST_DATA_DIR,"optical_dupes.bam"));
        args.addInput(new File(TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam"));
        args.addFlag(MarkDuplicatesSpark.ALLOW_MULTIPLE_SORT_ORDERS_IN_INPUT_ARG);
        runCommandLine(args);
    }

    @Test
    public void testNullOpticalDuplicates() {
        final File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(new File(TEST_DATA_DIR,"optical_dupes.bam"));
        args.addFlag(MarkDuplicatesSpark.ALLOW_MULTIPLE_SORT_ORDERS_IN_INPUT_ARG);
        args.add(OpticalDuplicatesArgumentCollection.READ_NAME_REGEX_LONG_NAME, "null");
        args.addRaw("--"+StandardArgumentDefinitions.METRICS_FILE_LONG_NAME);
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.addRaw(metricsFile.getAbsolutePath());
        runCommandLine(args);

        final MetricsFile<GATKDuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            System.err.println("Metrics file not found: " + ex);
        }
        final List<GATKDuplicationMetrics> nonEmptyMetrics = getGatkDuplicationMetrics(metricsOutput);

        // Assert that optical duplicates were not counted
        Assert.assertEquals(nonEmptyMetrics.get(0).READ_PAIR_DUPLICATES, 1);
        Assert.assertEquals(nonEmptyMetrics.get(0).READ_PAIR_OPTICAL_DUPLICATES, 0);

    }

    @Test
    public void testAssertAssumeUnsortedFilesAreQueryGroupedFiles() {
        final File output = createTempFile("supplementaryReadUnmappedMate", "bam");
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addOutput(output);
        args.addInput(new File(TEST_DATA_DIR,"optical_dupes.queryname.bam"));
        args.addInput(new File(TEST_DATA_DIR,"optical_dupes.unsorted.querygrouped.sam"));
        args.addFlag(MarkDuplicatesSpark.TREAT_UNSORTED_AS_ORDERED);
        runCommandLine(args);
    }
}
