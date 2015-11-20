package org.broadinstitute.hellbender.tools.picard.sam.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.TestUtil;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesTester;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.testers.AbstractMarkDuplicatesCommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.testers.AbstractMarkDuplicatesTester;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by MarkDuplicatesTester (see getTester).
 */
public final class MarkDuplicatesIntegrationTest extends AbstractMarkDuplicatesCommandLineProgramTest {
    protected static String TEST_BASE_NAME = null;

    @BeforeClass
    public void setUp() {
        TEST_BASE_NAME = "MarkDuplicates";
    }

    protected AbstractMarkDuplicatesTester getTester() {
        return new MarkDuplicatesTester();
    }

    @DataProvider(name="strictlyBadBams")
    public Object[][] strictlyBadBams() {
        return new Object[][]{
                {new File(getTestDataDir(), "non_strict.bam")},
        };
    }

    @Override
    protected CommandLineProgram getCommandLineProgramInstance() {
        return new MarkDuplicates();
    }

    // NB: this test should return different results than MarkDuplicatesWithMateCigar

    /**
     * Test that PG header records are created & chained appropriately (or not created), and that the PG record chains
     * are as expected.  MarkDuplicates is used both to merge and to mark dupes in this case.
     * @param suppressPg If true, do not create PG header record.
     * @param expectedPnVnByReadName For each read, info about the expect chain of PG records.
     */
    @Test(dataProvider = "pgRecordChainingTest")
    public void pgRecordChainingTest(final boolean suppressPg,
                                     final Map<String, List<ExpectedPnAndVn>> expectedPnVnByReadName) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        try {
            // Run MarkDuplicates, merging the 3 input files, and either enabling or suppressing PG header
            // record creation according to suppressPg.
            final MarkDuplicates markDuplicates = new MarkDuplicates();
            final ArgumentsBuilder args = new ArgumentsBuilder();
            for (int i = 1; i <= 3; ++i) {
                args.add("--input");
                args.add(new File(TEST_DATA_DIR, "merge" + i + ".sam").getAbsolutePath());
            }
            final File outputSam = new File(outputDir, TEST_BASE_NAME + ".sam");
            args.add("--output");
            args.add(outputSam.getAbsolutePath());
            args.add("--METRICS_FILE");
            args.add(new File(outputDir, TEST_BASE_NAME + ".duplicate_metrics").getAbsolutePath());
            if (suppressPg) {
                args.add("--PROGRAM_RECORD_ID");
                args.add("null");
            }

            // I generally prefer to call doWork rather than invoking the argument parser, but it is necessary
            // in this case to initialize the command line.
            // Note that for the unit test, version won't come through because it is obtained through jar
            // manifest, and unit test doesn't run code from a jar.
            markDuplicates.instanceMain(args.getArgsArray());

            // Read the MarkDuplicates output file, and get the PG ID for each read.  In this particular test,
            // the PG ID should be the same for both ends of a pair.
            final SamReader reader = SamReaderFactory.makeDefault().open(outputSam);

            final Map<String, String> pgIdForReadName = new HashMap<>();
            for (final SAMRecord rec : reader) {
                final String existingPgId = pgIdForReadName.get(rec.getReadName());
                final String thisPgId = rec.getStringAttribute(SAMTag.PG.name());
                if (existingPgId != null) {
                    Assert.assertEquals(thisPgId, existingPgId);
                } else {
                    pgIdForReadName.put(rec.getReadName(), thisPgId);
                }
            }
            final SAMFileHeader header = reader.getFileHeader();
            CloserUtil.close(reader);

            // Confirm that for each read name, the chain of PG records contains exactly the number that is expected,
            // and that values in the PG chain are as expected.
            for (final Map.Entry<String, List<ExpectedPnAndVn>> entry : expectedPnVnByReadName.entrySet()) {
                final String readName = entry.getKey();
                final List<ExpectedPnAndVn> expectedList = entry.getValue();
                String pgId = pgIdForReadName.get(readName);
                for (final ExpectedPnAndVn expected : expectedList) {
                    final SAMProgramRecord programRecord = header.getProgramRecord(pgId);
                    if (expected.expectedPn != null) Assert.assertEquals(programRecord.getProgramName(), expected.expectedPn);
                    if (expected.expectedVn != null) Assert.assertEquals(programRecord.getProgramVersion(), expected.expectedVn);
                    pgId = programRecord.getPreviousProgramGroupId();
                }
                Assert.assertNull(pgId);
            }

        } finally {
            TestUtil.recursiveDelete(outputDir);
        }
    }

    /**
     * Represents an expected PN value and VN value for a PG record.  If one of thexe is null, any value is allowed
     * in the PG record being tested.
     */
    private static class ExpectedPnAndVn {
        final String expectedPn;
        final String expectedVn;

        private ExpectedPnAndVn(final String expectedPn, final String expectedVn) {
            this.expectedPn = expectedPn;
            this.expectedVn = expectedVn;
        }
    }

    @DataProvider(name = "pgRecordChainingTest")
    public Object[][] pgRecordChainingTestDataProvider() {
        // Two test cases: One in which PG record generation is enabled, the other in which it is turned off.
        final Map<String, List<ExpectedPnAndVn>> withPgMap = new HashMap<>();
        withPgMap.put("1AAXX.1.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, null), new ExpectedPnAndVn(TEST_BASE_NAME, "1"), new ExpectedPnAndVn("bwa", "1")));
        withPgMap.put("1AAXX.2.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, null), new ExpectedPnAndVn("bwa", "2")));
        withPgMap.put("1AAXX.3.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, null)));

        final Map<String, List<ExpectedPnAndVn>> suppressPgMap = new HashMap<>();
        suppressPgMap .put("1AAXX.1.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, "1"), new ExpectedPnAndVn("bwa", "1")));
        suppressPgMap .put("1AAXX.2.1", Arrays.asList(new ExpectedPnAndVn("bwa", "2")));
        suppressPgMap .put("1AAXX.3.1", new ArrayList<>(0));
        return new Object[][] {
                { false, withPgMap},
                { true, suppressPgMap}
        };
    }

    @Test(dataProvider = "testOpticalDuplicateDetectionDataProvider")
    public void testOpticalDuplicateDetection(final File sam, final File reference, final String outputExtension, final long expectedNumOpticalDuplicates) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final File outputSam = new File(outputDir, TEST_BASE_NAME + outputExtension);
        outputSam.deleteOnExit();
        final File metricsFile = new File(outputDir, TEST_BASE_NAME + ".duplicate_metrics");
        metricsFile.deleteOnExit();
        // Run MarkDuplicates, merging the 3 input files, and either enabling or suppressing PG header
        // record creation according to suppressPg.
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.setupOpticalDuplicateFinder();
        markDuplicates.INPUT = CollectionUtil.makeList(sam);
        if (null != reference) {
            markDuplicates.REFERENCE_SEQUENCE = reference;
        }
        markDuplicates.OUTPUT = outputSam;
        markDuplicates.METRICS_FILE = metricsFile;
        markDuplicates.TMP_DIR = CollectionUtil.makeList(outputDir);
        // Needed to suppress calling CommandLineProgram.getVersion(), which doesn't work for code not in a jar
        markDuplicates.PROGRAM_RECORD_ID = null;
        Assert.assertEquals(markDuplicates.doWork(), null);
        Assert.assertEquals(markDuplicates.numOpticalDuplicates(), expectedNumOpticalDuplicates);
    }

    @Test(dataProvider = "strictlyBadBams")
    public void testLenientStringency(File input) {
        // TODO: This test should really be a PicardCommandLineProgramTest, not here (see #1170).
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(input.getPath());
        args.add("--"+"VALIDATION_STRINGENCY");
        args.add(ValidationStringency.LENIENT);
        args.add("--"+StandardArgumentDefinitions.OUTPUT_SHORT_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.add(outputFile.getAbsolutePath());

        args.add("--METRICS_FILE");
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.add(metricsFile.getAbsolutePath());

        runCommandLine(args.getArgsArray());
        // We don't care about the results, only that the program doesn't crash because of stringency.
    }

    @Test(dataProvider = "strictlyBadBams", expectedExceptions = SAMFormatException.class)
    public void testStrictStringency(File input) {
        // TODO: This test should really be a PicardCommandLineProgramTest, not here (see #1170).
        // For these inputs, we expect an expection to be thrown because the bams have some issues.
        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(input.getPath());
        args.add("--"+"VALIDATION_STRINGENCY");
        args.add(ValidationStringency.STRICT);
        args.add("--"+StandardArgumentDefinitions.OUTPUT_SHORT_NAME);

        File outputFile = createTempFile("markdups", ".bam");
        outputFile.delete();
        args.add(outputFile.getAbsolutePath());

        args.add("--METRICS_FILE");
        File metricsFile = createTempFile("markdups_metrics", ".txt");
        args.add(metricsFile.getAbsolutePath());

        runCommandLine(args.getArgsArray());
        // We don't care about the results, only that the program it throws an exception because of stringency.
    }

    @DataProvider(name="testOpticalDuplicateDetectionDataProvider")
    public Object[][] testOpticalDuplicateDetectionDataProvider() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "optical_dupes.sam"), null, ".sam", 1L},
                {new File(TEST_DATA_DIR, "optical_dupes.cram"), new File(TEST_DATA_DIR, "optical_dupes.fasta"), ".cram", 1L},
                {new File(TEST_DATA_DIR, "optical_dupes_casava.sam"), null, ".sam", 1L}
        };
    }

    @Test(dataProvider = "testDuplicateDetectionDataProvider")
    public void testDuplicateDetection(final File sam, final long expectedNumOpticalDuplicates) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final File outputSam = new File(outputDir, TEST_BASE_NAME + ".sam");
        outputSam.deleteOnExit();
        final File metricsFile = new File(outputDir, TEST_BASE_NAME + ".duplicate_metrics");
        metricsFile.deleteOnExit();
        // Run MarkDuplicates, merging the 3 input files, and either enabling or suppressing PG header
        // record creation according to suppressPg.
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.setupOpticalDuplicateFinder();
        markDuplicates.INPUT = CollectionUtil.makeList(sam);
        markDuplicates.OUTPUT = outputSam;
        markDuplicates.METRICS_FILE = metricsFile;
        markDuplicates.TMP_DIR = CollectionUtil.makeList(outputDir);
        // Needed to suppress calling CommandLineProgram.getVersion(), which doesn't work for code not in a jar
        markDuplicates.PROGRAM_RECORD_ID = null;
        Assert.assertEquals(markDuplicates.doWork(), null);
        Assert.assertEquals(markDuplicates.numOpticalDuplicates(), expectedNumOpticalDuplicates);
    }

    @DataProvider(name="testDuplicateDetectionDataProvider")
    public Object[][] testDuplicateDetectionDataProvider() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "markDups.test2reads.bam"), 0L},
                {new File(TEST_DATA_DIR, "example.chr1.1-1K.unmarkedDups.noDups.bam"), 0L},
                {new File(TEST_DATA_DIR, "example.chr1.1-1K.markedDups.bam"), 0L},
                {new File(TEST_DATA_DIR, "example.chr1.1-1K.unmarkedDups.bam"), 0L},
        };
    }


}
