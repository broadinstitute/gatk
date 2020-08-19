package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.*;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.function.Function;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.testutils.XorWrapper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class ReadsPathDataSourceUnitTest extends GATKBaseTest {

    private static final String FIRST_TEST_BAM_NAME = "reads_data_source_test1.bam";
    private static final String SECOND_TEST_BAM_NAME = "reads_data_source_test2.bam";
    private static final String THIRD_TEST_BAM_NAME = "reads_data_source_test3.bam";

    // This bam has mapped reads from various contigs, plus a few unmapped reads with no mapped mate
    private static final String UNMAPPED_TEST_BAM_NAME = "reads_data_source_test1_with_unmapped.bam";

    // This is a snippet of the CEUTrio.HiSeq.WGS.b37.NA12878 bam from large, with mapped reads
    // from chromosome 20 (with one mapped read having an unmapped mate), plus several unmapped
    // reads with no mapped mate.
    private static final String CEU_SNIPPET_NAME = "CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam";

    private static final String READS_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private final GATKPath FIRST_TEST_BAM = new GATKPath(READS_DATA_SOURCE_TEST_DIRECTORY + FIRST_TEST_BAM_NAME);
    private final GATKPath SECOND_TEST_BAM = new GATKPath(READS_DATA_SOURCE_TEST_DIRECTORY + SECOND_TEST_BAM_NAME);
    private final GATKPath THIRD_TEST_BAM = new GATKPath(READS_DATA_SOURCE_TEST_DIRECTORY + THIRD_TEST_BAM_NAME);
    private final GATKPath FIRST_TEST_SAM = new GATKPath(READS_DATA_SOURCE_TEST_DIRECTORY + "invalid_coord_sort_order.sam");

    private final GATKPath UNMAPPED_TEST_BAM = new GATKPath(READS_DATA_SOURCE_TEST_DIRECTORY + UNMAPPED_TEST_BAM_NAME);
    private final GATKPath CEU_SNIPPET_BAM = new GATKPath(READS_DATA_SOURCE_TEST_DIRECTORY + CEU_SNIPPET_NAME);

    private static final String ENDPOINT = "htsget://127.0.0.1:3000/reads/";
    private static final String LOCAL_PREFIX = "gatktest.";

    private final GATKPath FIRST_TEST_BAM_HTSGET = new GATKPath(ENDPOINT + LOCAL_PREFIX + FIRST_TEST_BAM_NAME);
    private final GATKPath SECOND_TEST_BAM_HTSGET = new GATKPath(ENDPOINT + LOCAL_PREFIX + SECOND_TEST_BAM_NAME);
    private final GATKPath THIRD_TEST_BAM_HTSGET = new GATKPath(ENDPOINT + LOCAL_PREFIX + THIRD_TEST_BAM_NAME);
    private final GATKPath UNMAPPED_TEST_BAM_HTSGET = new GATKPath(ENDPOINT + LOCAL_PREFIX + UNMAPPED_TEST_BAM_NAME);
    private final GATKPath CEU_SNIPPET_BAM_HTSGET = new GATKPath(ENDPOINT + LOCAL_PREFIX + CEU_SNIPPET_NAME);


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFile() {
        new ReadsPathDataSource((GATKPath) null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFileList() {
        new ReadsPathDataSource((List<GATKPath>) null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleEmptyFileList() {
        new ReadsPathDataSource(Collections.emptyList());
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testHandleNonExistentFile() {
        new ReadsPathDataSource(new GATKPath(GATKBaseTest.getSafeNonExistentPath("nonexistent.bam")));
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleUnindexedFileWithIntervals() {
        // Cannot initialize a reads source with intervals unless all files are indexed
        final Path unindexed = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "unindexed.bam");
        Assert.assertNull(SamFiles.findIndex(unindexed), "Expected file to have no index, but found an index file. " + unindexed.toAbsolutePath());
        final ReadsDataSource readsSource = new ReadsPathDataSource(new GATKPath(unindexed));
        readsSource.setTraversalBounds(Collections.singletonList(new SimpleInterval("1", 1, 5)));
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleUnindexedFileQuery() {
        // Construction should succeed, since we don't pass in any intervals, but the query should throw.
        final Path unindexed = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "unindexed.bam");
        Assert.assertNull(SamFiles.findIndex(unindexed), "Expected file to have no index, but found an index file" + unindexed.toAbsolutePath());
        final ReadsDataSource readsSource = new ReadsPathDataSource(new GATKPath(unindexed));
        readsSource.query(new SimpleInterval("1", 1, 5));
    }

    @Test
    public void testDefaultSamReaderValidationStringency() {
        // Default validation stringency = SILENT results in no validation errors on invalid coordinate sort
        final ReadsDataSource readsSource = new ReadsPathDataSource(FIRST_TEST_SAM);
        //noinspection StatementWithEmptyBody
        for (@SuppressWarnings("unused") final GATKRead read : readsSource) {
        }
    }

    @Test(expectedExceptions = SAMFormatException.class)
    public void testCustomSamReaderFactory() {
        // Custom SamReaderFactory with validation stringency = STRICT fails on invalid coordinate sort
        final ReadsDataSource readsSource = new ReadsPathDataSource(
            FIRST_TEST_SAM,
            SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT));
        //noinspection StatementWithEmptyBody
        for (@SuppressWarnings("unused") final GATKRead read : readsSource) {
        }
    }

    @DataProvider(name = "supportsSerialIteration")
    public Object[][] supportsSerialIteration() {
        // Files, expected to support serial iteration (false if any input is a .sam)
        return new Object[][]{
            {Collections.singletonList(FIRST_TEST_BAM), true},
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM), true},
            {Collections.singletonList(FIRST_TEST_SAM), false},
            {Arrays.asList(FIRST_TEST_BAM, FIRST_TEST_SAM), false},
            {Arrays.asList(FIRST_TEST_BAM, FIRST_TEST_BAM_HTSGET), true}
        };
    }

    @Test(dataProvider = "supportsSerialIteration")
    public void testSupportsSerialIteration(final List<GATKPath> inputs, final boolean expectedSupportsSerialIteration) {
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(inputs)) {
            Assert.assertEquals(readsSource.supportsSerialIteration(), expectedSupportsSerialIteration);
        }
    }

    @DataProvider(name = "singleFileCompleteTraversalData")
    public Object[][] singleFileCompleteTraversalData() {
        // Files, with expected read names in the expected order
        return new Object[][]{
            {FIRST_TEST_BAM, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")},
            {SECOND_TEST_BAM, Arrays.asList("l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v")},
            {THIRD_TEST_BAM, Arrays.asList("w", "x", "y", "z")},
            {FIRST_TEST_BAM_HTSGET, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")},
            {SECOND_TEST_BAM_HTSGET, Arrays.asList("l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v")},
            {THIRD_TEST_BAM_HTSGET, Arrays.asList("w", "x", "y", "z")},
        };
    }

    @Test(dataProvider = "singleFileCompleteTraversalData")
    public void testSingleFileCompleteTraversal(final GATKPath samFile, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(samFile)) {
            traverseOnce(readsSource, samFile, expectedReadNames);
        }
    }

    @Test(dataProvider = "singleFileCompleteTraversalData")
    public void testSingleFileSerialTraversal(final GATKPath samFile, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(samFile)) {
            Assert.assertTrue(readsSource.supportsSerialIteration());

            traverseOnce(readsSource, samFile, expectedReadNames);
            traverseOnce(readsSource, samFile, expectedReadNames);
            traverseOnce(readsSource, samFile, expectedReadNames);
        }
    }

    private static void traverseOnce(final ReadsDataSource readsSource, final GATKPath samFile, final List<String> expectedReadNames) {
        final List<GATKRead> reads = new ArrayList<>();
        for (final GATKRead read : readsSource) {
            reads.add(read);
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in complete traversal of " + samFile);

        // Make sure we got the reads we expected in the right order
        for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
            Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in complete traversal of " + samFile.toString() + " not equal to expected read");
        }
    }

    @DataProvider(name = "singleFileTraversalWithIntervalsData")
    public Object[][] singleFileTraversalWithIntervalsData() {
        // Files, with intervals, and expected read names in the expected order
        return new Object[][]{
            {FIRST_TEST_BAM,
                Arrays.asList(new SimpleInterval("1", 200, 210), new SimpleInterval("2", 550, 700), new SimpleInterval("4", 700, 701)),
                Arrays.asList("a", "b", "c", "f", "g", "h", "k")
            },
            {FIRST_TEST_BAM,
                Arrays.asList(new SimpleInterval("1", 205, 209), new SimpleInterval("3", 400, 410)),
                Arrays.asList("a", "b", "j")
            },
            {FIRST_TEST_BAM,
                Arrays.asList(new SimpleInterval("1", 999, 1200), new SimpleInterval("2", 530, 625), new SimpleInterval("4", 1000, 1200)),
                Arrays.asList("d", "e", "f", "g")
            },
            {FIRST_TEST_BAM,
                Arrays.asList(new SimpleInterval("1", 900, 1100), new SimpleInterval("1", 1000, 1200)),
                Arrays.asList("d", "e")
            },
            {FIRST_TEST_BAM,
                Collections.singletonList(new SimpleInterval("1", 1000, 1099)),
                Collections.singletonList("d")
            },
            {FIRST_TEST_BAM,
                Arrays.asList(new SimpleInterval("1", 2000, 3000), new SimpleInterval("1", 4000, 5000), new SimpleInterval("2", 1000, 2000)),
                Collections.emptyList()
            },
            {FIRST_TEST_BAM_HTSGET,
                Arrays.asList(new SimpleInterval("1", 200, 210), new SimpleInterval("2", 550, 700), new SimpleInterval("4", 700, 701)),
                Arrays.asList("a", "b", "c", "f", "g", "h", "k")
            },
            {FIRST_TEST_BAM_HTSGET,
                Arrays.asList(new SimpleInterval("1", 205, 209), new SimpleInterval("3", 400, 410)),
                Arrays.asList("a", "b", "j")
            },
            {FIRST_TEST_BAM_HTSGET,
                Arrays.asList(new SimpleInterval("1", 999, 1200), new SimpleInterval("2", 530, 625), new SimpleInterval("4", 1000, 1200)),
                Arrays.asList("d", "e", "f", "g")
            },
            {FIRST_TEST_BAM_HTSGET,
                Arrays.asList(new SimpleInterval("1", 900, 1100), new SimpleInterval("1", 1000, 1200)),
                Arrays.asList("d", "e")
            },
            {FIRST_TEST_BAM_HTSGET,
                Collections.singletonList(new SimpleInterval("1", 1000, 1099)),
                Collections.singletonList("d")
            },
            {FIRST_TEST_BAM_HTSGET,
                Arrays.asList(new SimpleInterval("1", 2000, 3000), new SimpleInterval("1", 4000, 5000), new SimpleInterval("2", 1000, 2000)),
                Collections.emptyList()
            },
        };
    }

    @Test(dataProvider = "singleFileTraversalWithIntervalsData")
    public void testSingleFileTraversalWithIntervals(final GATKPath samFile, final List<SimpleInterval> intervals, final List<String> expectedReadNames) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(samFile)) {
            Assert.assertTrue(readsSource.isQueryableByInterval(), "This reads source should be queryable by interval");

            readsSource.setTraversalBounds(intervals);

            final List<GATKRead> reads = new ArrayList<>();
            for (final GATKRead read : readsSource) {
                reads.add(read);
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFile);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFile + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "singleFileQueryByIntervalData")
    public Object[][] singleFileQueryByIntervalData() {
        // Files, with a single query interval, and expected read names in the expected order
        return new Object[][]{
            {FIRST_TEST_BAM,
                new SimpleInterval("1", 200, 209),
                Arrays.asList("a", "b")
            },
            {FIRST_TEST_BAM,
                new SimpleInterval("1", 285, 1100),
                Arrays.asList("c", "d", "e")
            },
            {FIRST_TEST_BAM,
                new SimpleInterval("2", 550, 649),
                Arrays.asList("f", "g")
            },
            {FIRST_TEST_BAM,
                new SimpleInterval("3", 399, 400),
                Collections.singletonList("j")
            },
            {FIRST_TEST_BAM,
                new SimpleInterval("4", 100, 200),
                Collections.emptyList()
            },
            {FIRST_TEST_BAM,
                new SimpleInterval("4", 600, 699),
                Collections.emptyList()
            },
            {FIRST_TEST_BAM_HTSGET,
                new SimpleInterval("1", 200, 209),
                Arrays.asList("a", "b")
            },
            {FIRST_TEST_BAM_HTSGET,
                new SimpleInterval("1", 285, 1100),
                Arrays.asList("c", "d", "e")
            },
            {FIRST_TEST_BAM_HTSGET,
                new SimpleInterval("2", 550, 649),
                Arrays.asList("f", "g")
            },
            {FIRST_TEST_BAM_HTSGET,
                new SimpleInterval("3", 399, 400),
                Collections.singletonList("j")
            },
            {FIRST_TEST_BAM_HTSGET,
                new SimpleInterval("4", 100, 200),
                Collections.emptyList()
            },
            {FIRST_TEST_BAM_HTSGET,
                new SimpleInterval("4", 600, 699),
                Collections.emptyList()
            },
        };
    }

    @Test(dataProvider = "singleFileQueryByIntervalData")
    public void testSingleFileQueryByInterval(final GATKPath samFile, final SimpleInterval interval, final List<String> expectedReadNames) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(samFile)) {
            Assert.assertTrue(readsSource.isQueryableByInterval(), "This reads source should be queryable by interval");
            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
        }
    }

    @Test(dataProvider = "singleFileQueryByIntervalData")
    public void testSingleFileQueryByIntervalSerialIteration(final GATKPath samFile, final SimpleInterval interval, final List<String> expectedReadNames) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(samFile)) {
            Assert.assertTrue(readsSource.isQueryableByInterval(), "This reads source should be queryable by interval");
            Assert.assertTrue(readsSource.supportsSerialIteration());

            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
        }
    }

    private static void traverseOnceByInterval(final ReadsDataSource readsSource, final GATKPath samFile, final SimpleInterval interval, final List<String> expectedReadNames) {
        final List<GATKRead> reads = new ArrayList<>();
        final Iterator<GATKRead> queryIterator = readsSource.query(interval);
        while (queryIterator.hasNext()) {
            reads.add(queryIterator.next());
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in query by interval of " + samFile);

        // Make sure we got the reads we expected in the right order
        for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
            Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in query by interval of " + samFile + " not equal to expected read");
        }
    }

    @DataProvider(name = "multipleFilesCompleteTraversalData")
    public Object[][] multipleFilesCompleteTraversalData() {
        // Files, with expected read names in the expected order
        return new Object[][]{
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM), Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "t", "i", "j", "u", "v", "k")},
            {Arrays.asList(SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.asList("l", "m", "n", "o", "p", "q", "r", "s", "w", "t", "x", "u", "v", "y", "z")},
            {Arrays.asList(FIRST_TEST_BAM, THIRD_TEST_BAM), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "w", "x", "i", "j", "y", "k", "z")},
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")},
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET), Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "t", "i", "j", "u", "v", "k")},
            {Arrays.asList(SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET), Arrays.asList("l", "m", "n", "o", "p", "q", "r", "s", "w", "t", "x", "u", "v", "y", "z")},
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "w", "x", "i", "j", "y", "k", "z")},
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET), Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")},
        };
    }

    @Test(dataProvider = "multipleFilesCompleteTraversalData")
    public void testMultipleFilesCompleteTraversal(final List<GATKPath> samFiles, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(samFiles)) {
            final List<GATKRead> reads = new ArrayList<>();

            for (final GATKRead read : readsSource) {
                reads.add(read);
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in complete traversal of " + samFiles);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in complete traversal of " + samFiles + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "MultipleFilesTraversalWithIntervalsData")
    public Object[][] getMultipleFilesTraversalWithIntervalsData() {
        // Files, with intervals, and expected read names in the expected order
        return new Object[][]{
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                Arrays.asList(new SimpleInterval("1", 205, 207), new SimpleInterval("1", 400, 1000), new SimpleInterval("4", 500, 704)),
                Arrays.asList("a", "b", "l", "n", "d", "y", "k")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                Arrays.asList(new SimpleInterval("4", 500, 704), new SimpleInterval("1", 400, 1000), new SimpleInterval("1", 205, 207)),
                Arrays.asList("a", "b", "l", "n", "d", "y", "k")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                Arrays.asList(new SimpleInterval("2", 500, 600), new SimpleInterval("2", 2099, 2200), new SimpleInterval("3", 50, 100), new SimpleInterval("3", 300, 500)),
                Arrays.asList("f", "p", "g", "s", "w", "i", "j", "u")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                Arrays.asList(new SimpleInterval("1", 1, 300), new SimpleInterval("1", 100, 500), new SimpleInterval("1", 200, 600)),
                Arrays.asList("a", "b", "l", "c", "m", "n")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                Arrays.asList(new SimpleInterval("1", 11000, 12000), new SimpleInterval("3", 1000, 2000)),
                Collections.emptyList()
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                Arrays.asList(new SimpleInterval("1", 1, 16000), new SimpleInterval("2", 1, 16000), new SimpleInterval("3", 1, 16000), new SimpleInterval("4", 1, 16000)),
                Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                Arrays.asList(new SimpleInterval("1", 205, 207), new SimpleInterval("1", 400, 1000), new SimpleInterval("4", 500, 704)),
                Arrays.asList("a", "b", "l", "n", "d", "y", "k")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                Arrays.asList(new SimpleInterval("4", 500, 704), new SimpleInterval("1", 400, 1000), new SimpleInterval("1", 205, 207)),
                Arrays.asList("a", "b", "l", "n", "d", "y", "k")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                Arrays.asList(new SimpleInterval("2", 500, 600), new SimpleInterval("2", 2099, 2200), new SimpleInterval("3", 50, 100), new SimpleInterval("3", 300, 500)),
                Arrays.asList("f", "p", "g", "s", "w", "i", "j", "u")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                Arrays.asList(new SimpleInterval("1", 1, 300), new SimpleInterval("1", 100, 500), new SimpleInterval("1", 200, 600)),
                Arrays.asList("a", "b", "l", "c", "m", "n")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                Arrays.asList(new SimpleInterval("1", 11000, 12000), new SimpleInterval("3", 1000, 2000)),
                Collections.emptyList()
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                Arrays.asList(new SimpleInterval("1", 1, 16000), new SimpleInterval("2", 1, 16000), new SimpleInterval("3", 1, 16000), new SimpleInterval("4", 1, 16000)),
                Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")
            },
        };
    }

    @Test(dataProvider = "MultipleFilesTraversalWithIntervalsData")
    public void testMultipleFilesTraversalWithIntervals(final List<GATKPath> samFiles, final List<SimpleInterval> intervals, final List<String> expectedReadNames) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(samFiles)) {
            Assert.assertTrue(readsSource.isQueryableByInterval(), "This reads source should be queryable by interval");

            readsSource.setTraversalBounds(intervals);

            final List<GATKRead> reads = new ArrayList<>();
            for (final GATKRead read : readsSource) {
                reads.add(read);
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFiles);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFiles + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "multipleFilesQueryByIntervalData")
    public Object[][] multipleFilesQueryByIntervalData() {
        // Files, with a single query interval, and expected read names in the expected order
        return new Object[][]{
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                new SimpleInterval("1", 285, 1000),
                Arrays.asList("c", "m", "n", "d")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                new SimpleInterval("3", 200, 300),
                Arrays.asList("t", "x", "i")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                new SimpleInterval("1", 9000, 11000),
                Collections.singletonList("o")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                new SimpleInterval("3", 1, 16000),
                Arrays.asList("w", "t", "x", "i", "j", "u", "v")
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                new SimpleInterval("2", 10000, 12000),
                Collections.emptyList()
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                new SimpleInterval("1", 285, 1000),
                Arrays.asList("c", "m", "n", "d")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                new SimpleInterval("3", 200, 300),
                Arrays.asList("t", "x", "i")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                new SimpleInterval("1", 9000, 11000),
                Collections.singletonList("o")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                new SimpleInterval("3", 1, 16000),
                Arrays.asList("w", "t", "x", "i", "j", "u", "v")
            },
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, SECOND_TEST_BAM_HTSGET, THIRD_TEST_BAM_HTSGET),
                new SimpleInterval("2", 10000, 12000),
                Collections.emptyList()
            },
        };
    }

    @Test(dataProvider = "multipleFilesQueryByIntervalData")
    public void testMultipleFilesQueryByInterval(final List<GATKPath> samFiles, final SimpleInterval interval, final List<String> expectedReadNames) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(samFiles)) {
            Assert.assertTrue(readsSource.isQueryableByInterval(), "This reads source should be queryable by interval");

            final List<GATKRead> reads = new ArrayList<>();
            final Iterator<GATKRead> queryIterator = readsSource.query(interval);
            while (queryIterator.hasNext()) {
                reads.add(queryIterator.next());
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in query by interval of " + samFiles);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in query by interval of " + samFiles + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "traversalWithUnmappedReadsTestData")
    public Object[][] traversalWithUnmappedReadsTestData() {
        // This bam has only mapped reads
        final GATKPath mappedBam = FIRST_TEST_BAM;
        return new Object[][]{
            // One interval, no unmapped
            {UNMAPPED_TEST_BAM, Collections.singletonList(new SimpleInterval("1", 200, 1000)), false, Arrays.asList("a", "b", "c", "d")},
            // One interval, with unmapped
            {UNMAPPED_TEST_BAM, Collections.singletonList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d", "u1", "u2", "u3", "u4", "u5")},
            // Multiple intervals, no unmapped
            {UNMAPPED_TEST_BAM, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), false, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k")},
            // Multiple intervals, with unmapped
            {UNMAPPED_TEST_BAM, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), true, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k", "u1", "u2", "u3", "u4", "u5")},
            // Interval with no overlapping reads, no unmapped
            {UNMAPPED_TEST_BAM, Collections.singletonList(new SimpleInterval("1", 3000, 4000)), false, Collections.emptyList()},
            // Interval with no overlapping reads, with unmapped
            {UNMAPPED_TEST_BAM, Collections.singletonList(new SimpleInterval("1", 3000, 4000)), true, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            // Interval with no overlapping reads, with unmapped, but no unmapped reads in bam
            {mappedBam, Collections.singletonList(new SimpleInterval("1", 3000, 4000)), true, Collections.emptyList()},
            // Interval with overlapping reads, with unmapped, but no unmapped reads in bam
            {mappedBam, Collections.singletonList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d")},
            // Null intervals, with unmapped
            {UNMAPPED_TEST_BAM, null, true, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            // Empty intervals, with unmapped
            {UNMAPPED_TEST_BAM, Collections.emptyList(), true, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            // Null intervals, with unmapped, but no unmapped reads in bam
            {mappedBam, null, true, Collections.emptyList()},
            // Empty intervals, with unmapped, but no unmapped reads in bam
            {mappedBam, Collections.emptyList(), true, Collections.emptyList()},
            // Null intervals, no unmapped (an unbounded traversal, so we expect all the reads)
            {UNMAPPED_TEST_BAM, null, false, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "u1", "u2", "u3", "u4", "u5")},
            // Empty intervals, no unmapped (an unbounded traversal, so we expect all the reads)
            {UNMAPPED_TEST_BAM, Collections.emptyList(), false, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "u1", "u2", "u3", "u4", "u5")},
            // Interval containing mapped read with unmapped mate, no unmapped
            {CEU_SNIPPET_BAM, Collections.singletonList(new SimpleInterval("20", 10000011, 10000013)), false, Arrays.asList("a", "b", "c", "d", "e", "f", "f")},
            // Interval containing mapped read with unmapped mate, with unmapped
            {CEU_SNIPPET_BAM, Collections.singletonList(new SimpleInterval("20", 10000011, 10000013)), true, Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i")},
            // Interval not containing mapped read with unmapped mate, no unmapped
            {CEU_SNIPPET_BAM, Collections.singletonList(new SimpleInterval("20", 10000009, 10000011)), false, Arrays.asList("a", "b", "c", "d", "e")},
            // Interval not containing mapped read with unmapped mate, with unmapped
            {CEU_SNIPPET_BAM, Collections.singletonList(new SimpleInterval("20", 10000009, 10000011)), true, Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i")}
        };
    }

    @Test(dataProvider = "traversalWithUnmappedReadsTestData")
    public void testTraversalWithUnmappedReads(final GATKPath samFile, final List<SimpleInterval> queryIntervals, final boolean queryUnmapped, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(samFile)) {
            readsSource.setTraversalBounds(queryIntervals, queryUnmapped);

            if ((queryIntervals != null && !queryIntervals.isEmpty()) || queryUnmapped) {
                Assert.assertTrue(readsSource.traversalIsBounded());
            } else {
                Assert.assertFalse(readsSource.traversalIsBounded());
            }

            final List<GATKRead> reads = new ArrayList<>();
            for (final GATKRead read : readsSource) {
                reads.add(read);
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFile);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFile + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "queryUnmappedTestData")
    public Object[][] queryUnmappedTestData() {
        return new Object[][]{
            {UNMAPPED_TEST_BAM, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            {CEU_SNIPPET_BAM, Arrays.asList("g", "h", "h", "i", "i")},
            {UNMAPPED_TEST_BAM_HTSGET, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            {CEU_SNIPPET_BAM_HTSGET, Arrays.asList("g", "h", "h", "i", "i")},
        };
    }

    @Test(dataProvider = "queryUnmappedTestData")
    public void testQueryUnmapped(final GATKPath samFile, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(samFile)) {
            final List<GATKRead> reads = new ArrayList<>();
            final Iterator<GATKRead> queryIterator = readsSource.queryUnmapped();
            while (queryIterator.hasNext()) {
                reads.add(queryIterator.next());
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in queryUnmapped on " + samFile);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in queryUnmapped on " + samFile + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "mergedHeaderIntervalQueries")
    public Object[][] mergedHeaderQueries() {
        return new Object[][]{

            // Single interval that only overlaps one of the input sources
            {new SimpleInterval[]{new SimpleInterval("EXTRA_CONTIG_1", 1, 10)},
                new String[]{"EXTRA_CONTIG_1_READ"}, new int[]{4}},
            {new SimpleInterval[]{new SimpleInterval("EXTRA_CONTIG_2", 1, 10)},
                new String[]{"EXTRA_CONTIG_2_READ"}, new int[]{5}},

            // Multiple intervals, each of which only overlaps one of the input sources
            {new SimpleInterval[]{new SimpleInterval("EXTRA_CONTIG_1", 1, 10), new SimpleInterval("EXTRA_CONTIG_2", 1, 10)},
                new String[]{"EXTRA_CONTIG_1_READ", "EXTRA_CONTIG_2_READ"},
                new int[]{4, 5}},

            // Single interval doesn't overlap ANY input source
            {new SimpleInterval[]{new SimpleInterval("TOTALLY_FAKE_CONTIG", 1, 10)},
                new String[]{}, new int[]{}}
        };
    }

    @Test(dataProvider = "mergedHeaderIntervalQueries")
    public void testMergedQueryWithFileSpecificContigs(
        final SimpleInterval[] intervals,
        final String[] expectedReadNames,
        final int[] expectedSequenceIndex) {
        // create two files, each with a read referencing a sequence that is not present in the other
        final GATKPath testFile1 = getFileWithAddedContig(FIRST_TEST_BAM, "EXTRA_CONTIG_1", "test1", ".bam");
        final GATKPath testFile2 = getFileWithAddedContig(SECOND_TEST_BAM, "EXTRA_CONTIG_2", "test2", ".bam");

        try (final ReadsDataSource readsSource = new ReadsPathDataSource(Arrays.asList(testFile1, testFile2))) {
            final SAMFileHeader samHeader = readsSource.getHeader();
            final SAMSequenceDictionary sequenceDictionary = readsSource.getSequenceDictionary();

            // 4 sequences in the original dictionaries, plus 1 added to each of the 2 inputs == 6
            Assert.assertEquals(sequenceDictionary.getSequences().size(), 6);
            Assert.assertEquals(samHeader.getSequenceDictionary(), sequenceDictionary);

            readsSource.setTraversalBounds(Arrays.asList(intervals));
            int count = 0;
            for (final GATKRead read : readsSource) {
                Assert.assertEquals(expectedReadNames[count], read.getName());
                Assert.assertEquals(
                    ReadUtils.getReferenceIndex(read, samHeader),
                    expectedSequenceIndex[count]
                );
                count++;
            }
            Assert.assertEquals(count, expectedReadNames.length);
        }
    }

    // Copy the reads in the inputFile to the output file, adding an extra contig to the sequence dictionary
    // and a read referencing that contig
    private static GATKPath getFileWithAddedContig(
        final GATKPath inputPath,
        final String extraContig,
        final String outputName,
        @SuppressWarnings("SameParameterValue") final String extension) {
        final File outputFile = GATKBaseTest.createTempFile(outputName, extension);
        try (final ReadsDataSource readsSource = new ReadsPathDataSource(inputPath)) {
            final SAMFileHeader header = readsSource.getHeader();
            final SAMSequenceRecord fakeSequenceRec = new SAMSequenceRecord(extraContig, 100);
            header.addSequence(fakeSequenceRec);

            try (final SAMFileGATKReadWriter gatkReadWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(
                outputFile, null, header, true, true, false))) {
                for (final GATKRead read : readsSource) {
                    gatkReadWriter.addRead(read);
                }
                // use the contig name in the read name to make it easy to see where this read came from
                final SAMRecord samRec = new SAMRecord(header);
                samRec.setReadName(extraContig + "_READ");
                samRec.setReferenceName(extraContig);
                samRec.setAlignmentStart(5);
                samRec.setReadBases(new byte[]{'a', 'c', 'g', 't'});
                gatkReadWriter.addRead(new SAMRecordToGATKReadAdapter(samRec));
            }
        }
        return new GATKPath(outputFile);
    }

    @DataProvider(name = "manuallySpecifiedIndexTestData")
    public Object[][] manuallySpecifiedIndexTestData() {
        final String BAM_DIR = READS_DATA_SOURCE_TEST_DIRECTORY + "readIndexTest/";
        final String INDEX_DIR = BAM_DIR + "indices/";

        final List<Path> bams = Arrays.asList(IOUtils.getPath(BAM_DIR + "reads_data_source_test1.bam"),
                IOUtils.getPath(BAM_DIR + "reads_data_source_test2.bam"));

        final List<Path> indices = Arrays.asList(IOUtils.getPath(INDEX_DIR + "reads_data_source_test1.bam.bai"),
                IOUtils.getPath(INDEX_DIR + "reads_data_source_test2.bam.bai"));

        return new Object[][]{
                {bams, indices}
        };
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData")
    public void testManuallySpecifiedIndices(final List<GATKPath> bams, final List<GATKPath> indices) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(bams, indices)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while (queryReads.hasNext()) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 5, "Wrong number of reads returned in query");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData")
    public void testManuallySpecifiedIndicesWithCustomReaderFactory(final List<GATKPath> bams, final List<GATKPath> indices) {
        final SamReaderFactory customFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);

        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(bams, indices, customFactory)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while (queryReads.hasNext()) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 5, "Wrong number of reads returned in query");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData")
    public void testManuallySpecifiedIndicesWithCustomReaderFactoryAndNullWrappers(final List<GATKPath> bams, final List<GATKPath> indices) {
        final SamReaderFactory customFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
        // ReadsDataSource should not be using the wrapper since the files are not on the Google cloud.
        // So we pass this invalid wrapper: if the code tries to use it, it'll blow up.
        final Function<SeekableByteChannel, SeekableByteChannel> nullWrapper = (SeekableByteChannel) -> null;

        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(bams, indices, customFactory, nullWrapper, nullWrapper)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while (queryReads.hasNext()) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 5, "Wrong number of reads returned in query");
        }
    }

    @DataProvider(name = "manuallySpecifiedMixedIndexTestData")
    public Object[][] manuallySpecifiedMixedIndexTestData() {
        final String BAM_DIR = READS_DATA_SOURCE_TEST_DIRECTORY + "readIndexTest/";
        final String INDEX_DIR = BAM_DIR + "indices/";

        final GATKPath BAM_1 = new GATKPath(BAM_DIR + "reads_data_source_test1.bam");
        final GATKPath BAM_2 = new GATKPath(BAM_DIR + "reads_data_source_test2.bam");

        final GATKPath INDEX_1 = new GATKPath(INDEX_DIR + "reads_data_source_test1.bam.bai");
        final GATKPath INDEX_2 = new GATKPath(INDEX_DIR + "reads_data_source_test2.bam.bai");

        return new Object[][]{
            {Arrays.asList(BAM_1, SECOND_TEST_BAM_HTSGET), Collections.singletonList(INDEX_1)},
            {Arrays.asList(FIRST_TEST_BAM_HTSGET, BAM_2), Collections.singletonList(INDEX_2)},
        };
    }

    @Test(dataProvider = "manuallySpecifiedMixedIndexTestData")
    public void testMixedHtsgetAndFilePaths(final List<GATKPath> bams, final List<GATKPath> indices) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(bams, indices)) {
            Assert.assertFalse(readsSource.indicesAvailable(), "Expected source to not have indices but found indices");
            Assert.assertTrue(readsSource.isQueryableByInterval(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while (queryReads.hasNext()) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 5, "Wrong number of reads returned in query");
        }
    }

    @DataProvider(name = "cloudXorTestData")
    public Object[][] cloudXorTestData() {
        final String BAM_DIR = getGCPTestInputPath() + "org/broadinstitute/hellbender/engine/";

        final List<Path> bams = Collections.singletonList(IOUtils.getPath(BAM_DIR + "reads_data_source_test4.xor.bam"));

        final List<Path> indices = Collections.singletonList(IOUtils.getPath(BAM_DIR + "reads_data_source_test4.xor.bam.bai"));

        return new Object[][]{
            {bams, indices}
        };
    }

    @Test(dataProvider = "cloudXorTestData", groups = {"bucket"})
    public void testCloudBamWithCustomReaderFactoryAndWrappers(final List<GATKPath> bams, final List<GATKPath> indices) {
        final SamReaderFactory customFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
        // The input files are XOR'd with a constant. We use a wrapper to XOR it back.
        // If the code uses the wrong wrapper, or omits one, then the test will fail.
        final Function<SeekableByteChannel, SeekableByteChannel> xorData = XorWrapper.forKey((byte) 74);
        final Function<SeekableByteChannel, SeekableByteChannel> xorIndex = XorWrapper.forKey((byte) 80);

        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(bams, indices, customFactory, xorData, xorIndex)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while (queryReads.hasNext()) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 2, "Wrong number of reads returned in query");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData")
    public void testManuallySpecifiedIndicesNullIndexListOK(final List<GATKPath> bams, @SuppressWarnings("unused") final List<GATKPath> indices) {
        try (final ReadsPathDataSource readsSource = new ReadsPathDataSource(bams, (List<GATKPath>) null)) {
            Assert.assertFalse(readsSource.indicesAvailable(), "Bams not indexed and explicit indices not provided, but indicesAvailable() returns true");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData", expectedExceptions = UserException.class)
    public void testManuallySpecifiedIndicesEmptyIndexList(final List<GATKPath> bams, @SuppressWarnings("unused") final List<GATKPath> indices) {
        new ReadsPathDataSource(bams, Collections.emptyList());
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData", expectedExceptions = UserException.class)
    public void testManuallySpecifiedIndicesWrongNumberOfIndices(final List<GATKPath> bams, final List<GATKPath> indices) {
        final List<GATKPath> wrongIndices = new ArrayList<>();
        wrongIndices.add(indices.get(0)); // Add one index, but not the other

        new ReadsPathDataSource(bams, wrongIndices);
    }


    @DataProvider(name = "readHeaders")
    public Object[][] getHeadersForDetectOrder() {
        final SAMFileHeader unknown = new SAMFileHeader();
        unknown.setSortOrder(SAMFileHeader.SortOrder.unknown);
        final SAMFileHeader coordinate = new SAMFileHeader();
        coordinate.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMFileHeader queryname = new SAMFileHeader();
        queryname.setSortOrder(SAMFileHeader.SortOrder.queryname);

        return new Object[][]{
            // empty list
            {Collections.emptyList(), SAMFileHeader.SortOrder.unsorted},
            // single header
            {Collections.singletonList(unknown), SAMFileHeader.SortOrder.unknown},
            {Collections.singletonList(coordinate), SAMFileHeader.SortOrder.coordinate},
            {Collections.singletonList(queryname), SAMFileHeader.SortOrder.queryname},
            // equal header list
            {Arrays.asList(unknown, unknown), SAMFileHeader.SortOrder.unknown},
            {Arrays.asList(coordinate, coordinate), SAMFileHeader.SortOrder.coordinate},
            {Arrays.asList(queryname, queryname), SAMFileHeader.SortOrder.queryname},
            // different header list
            {Arrays.asList(unknown, coordinate, queryname), SAMFileHeader.SortOrder.unsorted},
            {Arrays.asList(coordinate, coordinate, queryname), SAMFileHeader.SortOrder.unsorted}
        };
    }

    @Test(dataProvider = "readHeaders")
    public void testIdentifySortOrder(final List<SAMFileHeader> headers, final SAMFileHeader.SortOrder expected) {
        Assert.assertEquals(ReadsPathDataSource.identifySortOrder(headers), expected);
    }
}
