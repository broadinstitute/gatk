package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.*;
import java.nio.channels.SeekableByteChannel;
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
import java.nio.file.Path;
import java.util.*;

public final class ReadsDataSourceUnitTest extends GATKBaseTest {
    private static final String READS_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private final Path FIRST_TEST_BAM = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "reads_data_source_test1.bam");
    private final Path SECOND_TEST_BAM = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "reads_data_source_test2.bam");
    private final Path THIRD_TEST_BAM = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "reads_data_source_test3.bam");
    private final Path FIRST_TEST_SAM = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "invalid_coord_sort_order.sam");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFile() {
        Path nullFile = null;
        ReadsDataSource readsSource = new ReadsDataSource(nullFile);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFileList() {
        List<Path> nullList = null;
        ReadsDataSource readsSource = new ReadsDataSource(nullList);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleEmptyFileList() {
        ReadsDataSource readsSource = new ReadsDataSource(new ArrayList<>());
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testHandleNonExistentFile() {
        new ReadsDataSource(GATKBaseTest.getSafeNonExistentPath("nonexistent.bam"));
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleUnindexedFileWithIntervals() {
        // Cannot initialize a reads source with intervals unless all files are indexed
        final Path unindexed = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "unindexed.bam");
        Assert.assertNull(SamFiles.findIndex(unindexed), "Expected file to have no index, but found an index file. " + unindexed.toAbsolutePath());
        ReadsDataSource readsSource = new ReadsDataSource(unindexed);
        readsSource.setTraversalBounds(Arrays.asList(new SimpleInterval("1", 1, 5)));
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleUnindexedFileQuery() {
        // Construction should succeed, since we don't pass in any intervals, but the query should throw.
        final Path unindexed = IOUtils.getPath(READS_DATA_SOURCE_TEST_DIRECTORY + "unindexed.bam");
        Assert.assertNull(SamFiles.findIndex(unindexed), "Expected file to have no index, but found an index file" + unindexed.toAbsolutePath());
        ReadsDataSource readsSource = new ReadsDataSource(unindexed);
        readsSource.query(new SimpleInterval("1", 1, 5));
    }

    @Test
    public void testDefaultSamReaderValidationStringency() {
        // Default validation stringency = SILENT results in no validation errors on invalid coordinate sort
        final ReadsDataSource readsSource = new ReadsDataSource(FIRST_TEST_SAM);
        //noinspection StatementWithEmptyBody
        for ( @SuppressWarnings("unused") final GATKRead read : readsSource ) {
        }
    }

    @Test(expectedExceptions = SAMFormatException.class)
    public void testCustomSamReaderFactory() {
        // Custom SamReaderFactory with validation stringency = STRICT fails on invalid coordinate sort
        final ReadsDataSource readsSource = new ReadsDataSource(
                FIRST_TEST_SAM,
                SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT));
        //noinspection StatementWithEmptyBody
        for ( @SuppressWarnings("unused") final GATKRead read : readsSource ) {
        }
    }

    @DataProvider(name = "SingleFileCompleteTraversalData")
    public Object[][] getSingleFileCompleteTraversalData() {
        // Files, with expected read names in the expected order
        return new Object[][] {
                { FIRST_TEST_BAM, Arrays.<String>asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k") },
                { SECOND_TEST_BAM, Arrays.<String>asList("l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v") },
                { THIRD_TEST_BAM, Arrays.<String>asList("w", "x", "y", "z") }
        };
    }

    @Test(dataProvider = "SingleFileCompleteTraversalData")
    public void testSingleFileCompleteTraversal( final Path samFile, final List<String> expectedReadNames ) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFile)) {
            List<GATKRead> reads = new ArrayList<>();
            for ( GATKRead read : readsSource ) {
                reads.add(read);
            }
    
            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in complete traversal of " + samFile.toAbsolutePath());
    
            // Make sure we got the reads we expected in the right order
            for ( int readIndex = 0; readIndex < reads.size(); ++readIndex ) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in complete traversal of " + samFile.toAbsolutePath() + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "SingleFileTraversalWithIntervalsData")
    public Object[][] getSingleFileTraversalWithIntervalsData() {
        // Files, with intervals, and expected read names in the expected order
        return new Object[][] {
                { FIRST_TEST_BAM,
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 200, 210), new SimpleInterval("2", 550, 700), new SimpleInterval("4", 700, 701)),
                  Arrays.<String>asList("a", "b", "c", "f", "g", "h", "k")
                },
                { FIRST_TEST_BAM,
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 205, 209), new SimpleInterval("3", 400, 410)),
                  Arrays.<String>asList("a", "b", "j")
                },
                { FIRST_TEST_BAM,
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 999, 1200), new SimpleInterval("2", 530, 625), new SimpleInterval("4", 1000, 1200)),
                  Arrays.<String>asList("d", "e", "f", "g")
                },
                { FIRST_TEST_BAM,
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 900, 1100), new SimpleInterval("1", 1000, 1200)),
                  Arrays.<String>asList("d", "e")
                },
                { FIRST_TEST_BAM,
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 1000, 1099)),
                  Arrays.<String>asList("d")
                },
                { FIRST_TEST_BAM,
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 2000, 3000), new SimpleInterval("1", 4000, 5000), new SimpleInterval("2", 1000, 2000)),
                  Arrays.<String>asList()
                }
        };
    }

    @Test(dataProvider = "SingleFileTraversalWithIntervalsData")
    public void testSingleFileTraversalWithIntervals( final Path samFile, final List<SimpleInterval> intervals, final List<String> expectedReadNames ) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFile)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Indices should be reported as available for this reads source");

            readsSource.setTraversalBounds(intervals);

            List<GATKRead> reads = new ArrayList<>();
            for ( GATKRead read : readsSource ) {
                reads.add(read);
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFile.toAbsolutePath());

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFile.toAbsolutePath() + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "SingleFileQueryByIntervalData")
    public Object[][] getSingleFileQueryByIntervalData() {
        // Files, with a single query interval, and expected read names in the expected order
        return new Object[][]{
                { FIRST_TEST_BAM,
                  new SimpleInterval("1", 200, 209),
                  Arrays.<String>asList("a", "b")
                },
                { FIRST_TEST_BAM,
                  new SimpleInterval("1", 285, 1100),
                  Arrays.<String>asList("c", "d", "e")
                },
                { FIRST_TEST_BAM,
                  new SimpleInterval("2", 550, 649),
                  Arrays.<String>asList("f", "g")
                },
                { FIRST_TEST_BAM,
                  new SimpleInterval("3", 399, 400),
                  Arrays.<String>asList("j")
                },
                { FIRST_TEST_BAM,
                  new SimpleInterval("4", 100, 200),
                  Arrays.<String>asList()
                },
                { FIRST_TEST_BAM,
                  new SimpleInterval("4", 600, 699),
                  Arrays.<String>asList()
                }
        };
    }

    @Test(dataProvider = "SingleFileQueryByIntervalData")
    public void testSingleFileQueryByInterval( final Path samFile, final SimpleInterval interval, final List<String> expectedReadNames ) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFile)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Indices should be reported as available for this reads source");

            List<GATKRead> reads = new ArrayList<>();
            Iterator<GATKRead> queryIterator = readsSource.query(interval);
            while (queryIterator.hasNext()) {
                reads.add(queryIterator.next());
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in query by interval of " + samFile.toAbsolutePath());

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in query by interval of " + samFile.toAbsolutePath() + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "MultipleFilesCompleteTraversalData")
    public Object[][] getMultipleFilesCompleteTraversalData() {
        // Files, with expected read names in the expected order
        return new Object[][] {
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM), Arrays.<String>asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "t", "i", "j", "u", "v", "k") },
                { Arrays.<Path>asList(SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.<String>asList("l", "m", "n", "o", "p", "q", "r", "s", "w", "t", "x", "u", "v", "y", "z") },
                { Arrays.<Path>asList(FIRST_TEST_BAM, THIRD_TEST_BAM), Arrays.<String>asList("a", "b", "c", "d", "e", "f", "g", "h", "w", "x", "i", "j", "y", "k", "z") },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.<String>asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z") }
        };
    }

    @Test(dataProvider = "MultipleFilesCompleteTraversalData")
    public void testMultipleFilesCompleteTraversal(final List<Path> samFiles, final List<String> expectedReadNames) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFiles)) {
            List<GATKRead> reads = new ArrayList<>();

            for (GATKRead read : readsSource) {
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
        return new Object[][] {
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 205, 207), new SimpleInterval("1", 400, 1000), new SimpleInterval("4", 500, 704)),
                  Arrays.<String>asList("a", "b", "l", "n", "d", "y", "k")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("4", 500, 704), new SimpleInterval("1", 400, 1000), new SimpleInterval("1", 205, 207)),
                  Arrays.<String>asList("a", "b", "l", "n", "d", "y", "k")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("2", 500, 600), new SimpleInterval("2", 2099, 2200), new SimpleInterval("3", 50, 100), new SimpleInterval("3", 300, 500)),
                  Arrays.<String>asList("f", "p", "g", "s", "w", "i", "j", "u")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 1, 300), new SimpleInterval("1", 100, 500), new SimpleInterval("1", 200, 600)),
                  Arrays.<String>asList("a", "b", "l", "c", "m", "n")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 11000, 12000), new SimpleInterval("3", 1000, 2000)),
                  Arrays.<String>asList()
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 1, 16000), new SimpleInterval("2", 1, 16000), new SimpleInterval("3", 1, 16000), new SimpleInterval("4", 1, 16000)),
                  Arrays.<String>asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")
                }
        };
    }

    @Test(dataProvider = "MultipleFilesTraversalWithIntervalsData")
    public void testMultipleFilesTraversalWithIntervals( final List<Path> samFiles, final List<SimpleInterval> intervals, final List<String> expectedReadNames ) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFiles)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Indices should be reported as available for this reads source");

            readsSource.setTraversalBounds(intervals);

            List<GATKRead> reads = new ArrayList<>();
            for (GATKRead read : readsSource) {
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

    @DataProvider(name = "MultipleFilesQueryByIntervalData")
    public Object[][] getMultipleFilesQueryByIntervalData() {
        // Files, with a single query interval, and expected read names in the expected order
        return new Object[][] {
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("1", 285, 1000),
                  Arrays.<String>asList("c", "m", "n", "d")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("3", 200, 300),
                  Arrays.<String>asList("t", "x", "i")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("1", 9000, 11000),
                  Arrays.<String>asList("o")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("3", 1, 16000),
                  Arrays.<String>asList("w", "t", "x", "i", "j", "u", "v")
                },
                { Arrays.<Path>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("2", 10000, 12000),
                  Arrays.<String>asList()
                }
        };
    }

    @Test(dataProvider = "MultipleFilesQueryByIntervalData")
    public void testMultipleFilesQueryByInterval( final List<Path> samFiles, final SimpleInterval interval, final List<String> expectedReadNames ) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFiles)) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Indices should be reported as available for this reads source");

            List<GATKRead> reads = new ArrayList<>();
            Iterator<GATKRead> queryIterator = readsSource.query(interval);
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

    @DataProvider(name = "TraversalWithUnmappedReadsTestData")
    public Object[][] traversalWithUnmappedReadsTestData() {
        // This bam has only mapped reads
        final Path mappedBam = IOUtils.getPath(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1.bam");

        // This bam has mapped reads from various contigs, plus a few unmapped reads with no mapped mate
        final Path unmappedBam = IOUtils.getPath(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_with_unmapped.bam");

        // This is a snippet of the CEUTrio.HiSeq.WGS.b37.NA12878 bam from large, with mapped reads
        // from chromosome 20 (with one mapped read having an unmapped mate), plus several unmapped
        // reads with no mapped mate.
        final Path ceuSnippet = IOUtils.getPath(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam");

        return new Object[][] {
                // One interval, no unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000)), false, Arrays.asList("a", "b", "c", "d") },
                // One interval, with unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d", "u1", "u2", "u3", "u4", "u5") },
                // Multiple intervals, no unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), false, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k") },
                // Multiple intervals, with unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), true, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k", "u1", "u2", "u3", "u4", "u5") },
                // Interval with no overlapping reads, no unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 3000, 4000)), false, Collections.<String>emptyList() },
                // Interval with no overlapping reads, with unmapped
                { unmappedBam, Arrays.asList(new SimpleInterval("1", 3000, 4000)), true, Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // Interval with no overlapping reads, with unmapped, but no unmapped reads in bam
                { mappedBam, Arrays.asList(new SimpleInterval("1", 3000, 4000)), true, Collections.<String>emptyList() },
                // Interval with overlapping reads, with unmapped, but no unmapped reads in bam
                { mappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d") },
                // Null intervals, with unmapped
                { unmappedBam, null, true, Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // Empty intervals, with unmapped
                { unmappedBam, Collections.<SimpleInterval>emptyList(), true, Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                // Null intervals, with unmapped, but no unmapped reads in bam
                { mappedBam, null, true, Collections.<String>emptyList() },
                // Empty intervals, with unmapped, but no unmapped reads in bam
                { mappedBam, Collections.<SimpleInterval>emptyList(), true, Collections.<String>emptyList() },
                // Null intervals, no unmapped (an unbounded traversal, so we expect all the reads)
                { unmappedBam, null, false, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "u1", "u2", "u3", "u4", "u5") },
                // Empty intervals, no unmapped (an unbounded traversal, so we expect all the reads)
                { unmappedBam, Collections.<SimpleInterval>emptyList(), false, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "u1", "u2", "u3", "u4", "u5") },
                // Interval containing mapped read with unmapped mate, no unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000011, 10000013)), false, Arrays.asList("a", "b", "c", "d", "e", "f", "f")},
                // Interval containing mapped read with unmapped mate, with unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000011, 10000013)), true, Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i")},
                // Interval not containing mapped read with unmapped mate, no unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000009, 10000011)), false, Arrays.asList("a", "b", "c", "d", "e")},
                // Interval not containing mapped read with unmapped mate, with unmapped
                { ceuSnippet, Arrays.asList(new SimpleInterval("20", 10000009, 10000011)), true, Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i")}
        };
    }

    @Test(dataProvider = "TraversalWithUnmappedReadsTestData")
    public void testTraversalWithUnmappedReads( final Path samFile, final List<SimpleInterval> queryIntervals, final boolean queryUnmapped, final List<String> expectedReadNames ) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFile)) {
            readsSource.setTraversalBounds(queryIntervals, queryUnmapped);

            if ( (queryIntervals != null && ! queryIntervals.isEmpty()) || queryUnmapped ) {
                Assert.assertTrue(readsSource.traversalIsBounded());
            }
            else {
                Assert.assertFalse(readsSource.traversalIsBounded());
            }

            List<GATKRead> reads = new ArrayList<>();
            for ( GATKRead read : readsSource ) {
                reads.add(read);
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFile.toAbsolutePath());

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFile.toAbsolutePath() + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "QueryUnmappedTestData")
    public Object[][] queryUnmappedTestData() {
        return new Object[][] {
                { IOUtils.getPath(publicTestDir + "org/broadinstitute/hellbender/engine/reads_data_source_test1_with_unmapped.bam"), Arrays.asList("u1", "u2", "u3", "u4", "u5") },
                { IOUtils.getPath(publicTestDir + "org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam"), Arrays.asList("g", "h", "h", "i", "i") }
        };
    }

    @Test(dataProvider = "QueryUnmappedTestData")
    public void testQueryUnmapped( final Path samFile, final List<String> expectedReadNames ) {
        try (ReadsDataSource readsSource = new ReadsDataSource(samFile)) {
            List<GATKRead> reads = new ArrayList<>();
            Iterator<GATKRead> queryIterator = readsSource.queryUnmapped();
            while (queryIterator.hasNext()) {
                reads.add(queryIterator.next());
            }

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in queryUnmapped on " + samFile.toAbsolutePath());

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in queryUnmapped on " + samFile.toAbsolutePath() + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "MergedHeaderIntervalQueries")
    public Object[][] mergedHeaderQueries() {
        return new Object[][] {

                // Single interval that only overlpas one of the input sources
                { new SimpleInterval[] {new SimpleInterval("EXTRA_CONTIG_1", 1,10)},
                         new String[] {"EXTRA_CONTIG_1_READ"}, new int[] {4} },
                { new SimpleInterval[] {new SimpleInterval("EXTRA_CONTIG_2", 1,10)},
                         new String[] {"EXTRA_CONTIG_2_READ"}, new int[] {5} },

                // Multiple intervals, each of which only overlaps one of the input sources
                { new SimpleInterval[] {new SimpleInterval("EXTRA_CONTIG_1", 1,10), new SimpleInterval("EXTRA_CONTIG_2", 1,10)},
                         new String[] {"EXTRA_CONTIG_1_READ", "EXTRA_CONTIG_2_READ"},
                         new int[] { 4, 5} },

                // Single interval doesn't overlap ANY input source
                { new SimpleInterval[] {new SimpleInterval("TOTALLY_FAKE_CONTIG", 1,10)},
                         new String[] {}, new int[] {} }
        };
    }

    @Test(dataProvider = "MergedHeaderIntervalQueries")
    public void testMergedQueryWithFileSpecificContigs(
            final SimpleInterval intervals[],
            final String expectedReadNames[],
            final int expectedSequenceIndex[]) {
        // create two files, each with a read referencing a sequence that is not present in the other
        final File testFile1 = getFileWithAddedContig(FIRST_TEST_BAM, "EXTRA_CONTIG_1", "test1", ".bam");
        final File testFile2 = getFileWithAddedContig(SECOND_TEST_BAM, "EXTRA_CONTIG_2", "test2", ".bam");

        try (final ReadsDataSource readsSource = new ReadsDataSource(Arrays.asList(testFile1.toPath(), testFile2.toPath()))) {
            SAMFileHeader samHeader = readsSource.getHeader();
            SAMSequenceDictionary sequenceDictionary = readsSource.getSequenceDictionary();

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
    private File getFileWithAddedContig(
            final Path inputPath,
            final String extraContig,
            final String outputName,
            final String extension) {
        final File outputFile = GATKBaseTest.createTempFile(outputName, extension);
        try (ReadsDataSource readsSource = new ReadsDataSource(inputPath)) {
            SAMFileHeader header = readsSource.getHeader();
            SAMSequenceRecord fakeSequenceRec = new SAMSequenceRecord(extraContig, 100);
            header.addSequence(fakeSequenceRec);

            try (final SAMFileGATKReadWriter gatkReadWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(
                    outputFile, null, header, true, true, false))) {
                for (final GATKRead read : readsSource) {
                    gatkReadWriter.addRead(read);
                }
                // use the contig name in the read name to make it easy to see where this read came from
                SAMRecord samRec = new SAMRecord(header);
                samRec.setReadName(extraContig + "_READ");
                samRec.setReferenceName(extraContig );
                samRec.setAlignmentStart(5);
                samRec.setReadBases(new byte[] {'a', 'c', 'g', 't'});
                gatkReadWriter.addRead(new SAMRecordToGATKReadAdapter(samRec));
            }
        }
        return outputFile;
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
    public void testManuallySpecifiedIndices( final List<Path> bams, final List<Path> indices ) {
        try ( final ReadsDataSource readsSource = new ReadsDataSource(bams, indices) ) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while ( queryReads.hasNext() ) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 5, "Wrong number of reads returned in query");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData")
    public void testManuallySpecifiedIndicesWithCustomReaderFactory( final List<Path> bams, final List<Path> indices ) {
        final SamReaderFactory customFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);

        try ( final ReadsDataSource readsSource = new ReadsDataSource(bams, indices, customFactory) ) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while ( queryReads.hasNext() ) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 5, "Wrong number of reads returned in query");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData")
    public void testManuallySpecifiedIndicesWithCustomReaderFactoryAndNullWrappers( final List<Path> bams, final List<Path> indices ) {
        final SamReaderFactory customFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
        // ReadsDataSource should not be using the wrapper since the files are not on the Google cloud.
        // So we pass this invalid wrapper: if the code tries to use it, it'll blow up.
        Function<SeekableByteChannel, SeekableByteChannel> nullWrapper = (SeekableByteChannel) -> null;

        try ( final ReadsDataSource readsSource = new ReadsDataSource(bams, indices, customFactory, nullWrapper, nullWrapper) ) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while ( queryReads.hasNext() ) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 5, "Wrong number of reads returned in query");
        }
    }

    @DataProvider(name="cloudXorTestData")
    public Object[][] cloudXorTestData() {
        final String BAM_DIR = getGCPTestInputPath() + "org/broadinstitute/hellbender/engine/";
        final String INDEX_DIR = BAM_DIR;

        final List<Path> bams = Arrays.asList(IOUtils.getPath(BAM_DIR + "reads_data_source_test4.xor.bam"));

        final List<Path> indices = Arrays.asList(IOUtils.getPath(INDEX_DIR + "reads_data_source_test4.xor.bam.bai"));

        return new Object[][]{
            {bams, indices}
        };
    }

    @Test(dataProvider = "cloudXorTestData", groups={"bucket"})
    public void testCloudBamWithCustomReaderFactoryAndWrappers( final List<Path> bams, final List<Path> indices ) {
        final SamReaderFactory customFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
        // The input files are XOR'd with a constant. We use a wrapper to XOR it back.
        // If the code uses the wrong wrapper, or omits one, then the test will fail.
        Function<SeekableByteChannel, SeekableByteChannel> xorData = XorWrapper.forKey((byte)74);
        Function<SeekableByteChannel, SeekableByteChannel> xorIndex = XorWrapper.forKey((byte)80);

        try ( final ReadsDataSource readsSource = new ReadsDataSource(bams, indices, customFactory, xorData, xorIndex) ) {
            Assert.assertTrue(readsSource.indicesAvailable(), "Explicitly-provided indices not detected for bams: " + bams);

            final Iterator<GATKRead> queryReads = readsSource.query(new SimpleInterval("1", 1, 300));
            int queryCount = 0;
            while ( queryReads.hasNext() ) {
                ++queryCount;
                queryReads.next();
            }
            Assert.assertEquals(queryCount, 2, "Wrong number of reads returned in query");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData")
    public void testManuallySpecifiedIndicesNullIndexListOK( final List<Path> bams, final List<Path> indices ) {
        try ( final ReadsDataSource readsSource = new ReadsDataSource(bams, (List<Path>)null) ) {
            Assert.assertFalse(readsSource.indicesAvailable(), "Bams not indexed and explicit indices not provided, but indicesAvailable() returns true");
        }
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData", expectedExceptions = UserException.class)
    public void testManuallySpecifiedIndicesEmptyIndexList( final List<Path> bams, final List<Path> indices ) {
        final ReadsDataSource readsSource = new ReadsDataSource(bams, Collections.emptyList());
    }

    @Test(dataProvider = "manuallySpecifiedIndexTestData", expectedExceptions = UserException.class)
    public void testManuallySpecifiedIndicesWrongNumberOfIndices( final List<Path> bams, final List<Path> indices ) {
        final List<Path> wrongIndices = new ArrayList<>();
        wrongIndices.add(indices.get(0)); // Add one index, but not the other

        final ReadsDataSource readsSource = new ReadsDataSource(bams, wrongIndices);
    }


    @DataProvider(name = "readHeaders")
    public Object[][] getHeadersForDetectOrder() {
        final SAMFileHeader unknown = new SAMFileHeader();
        unknown.setSortOrder(SAMFileHeader.SortOrder.unknown);
        final SAMFileHeader coordinate = new SAMFileHeader();
        coordinate.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMFileHeader queryname = new SAMFileHeader();
        queryname.setSortOrder(SAMFileHeader.SortOrder.queryname);

        return new Object[][] {
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
        Assert.assertEquals(ReadsDataSource.identifySortOrder(headers), expected);
    }
}
