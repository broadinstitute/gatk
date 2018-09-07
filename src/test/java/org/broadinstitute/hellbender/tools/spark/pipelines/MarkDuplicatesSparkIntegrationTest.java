package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSpark;
import org.broadinstitute.hellbender.tools.walkers.markduplicates.AbstractMarkDuplicatesCommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.testers.MarkDuplicatesSparkTester;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

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
            {new File(TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.noDups.bam"), 20, 0,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(0L, 3L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16406", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(0L, 1L, 0L, 0L, 0L, 0L, 0.0, 0L))},
            {new File(TEST_DATA_DIR,"example.chr1.1-1K.unmarkedDups.bam"), 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File(TEST_DATA_DIR,"example.chr1.1-1K.markedDups.bam"), 90, 6,
             ImmutableMap.of("Solexa-16419", ImmutableList.of(4L, 4L, 4L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16416", ImmutableList.of(2L, 2L, 2L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16404", ImmutableList.of(3L, 9L, 3L, 0L, 2L, 0L, 0.190476, 17L),
                             "Solexa-16406", ImmutableList.of(1L, 10L, 1L, 0L, 0L, 0L, 0.0, 0L),
                             "Solexa-16412", ImmutableList.of(3L, 6L, 3L, 0L, 1L, 0L, 0.133333, 15L))},
            {new File(TEST_DATA_DIR, "optical_dupes.bam"), 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
            {new File(TEST_DATA_DIR, "optical_dupes_casava.bam"), 4, 2,
             ImmutableMap.of("mylib", ImmutableList.of(0L, 2L, 0L, 0L, 1L, 1L, 0.5, 0L))},
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
}
