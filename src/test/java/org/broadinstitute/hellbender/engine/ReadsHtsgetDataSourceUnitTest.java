package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.util.*;

public class ReadsHtsgetDataSourceUnitTest extends GATKBaseTest {
    private static final String READS_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";

    private static final String ENDPOINT = "htsget://0.0.0.0:3000/reads/";
    private static final String LOCAL_PREFIX = "gatktest.";

    private final GATKPath FIRST_TEST_BAM = new GATKPath(ENDPOINT + LOCAL_PREFIX + "reads_data_source_test1.bam");
    private final GATKPath SECOND_TEST_BAM = new GATKPath(ENDPOINT + LOCAL_PREFIX + "reads_data_source_test2.bam");
    private final GATKPath THIRD_TEST_BAM = new GATKPath(ENDPOINT + LOCAL_PREFIX + "reads_data_source_test3.bam");

    private final GATKPath FIRST_TEST_UNMAPPED_BAM = new GATKPath(ENDPOINT + LOCAL_PREFIX + "reads_data_source_test1_with_unmapped.bam");
    private final GATKPath CEU_SNIPPET = new GATKPath(ENDPOINT + LOCAL_PREFIX + "CEUTrio.HiSeq.WGS.b37.NA12878.snippet_with_unmapped.bam");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFile() {
        new ReadsHtsgetDataSource((GATKPath) null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFileList() {
        new ReadsHtsgetDataSource((List<GATKPath>) null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleEmptyFileList() {
        new ReadsHtsgetDataSource(Collections.emptyList());
    }

    @Test(expectedExceptions = UserException.class)
    public void testQueryAfterAlreadyStopped() {
        final ReadsHtsgetDataSource source = new ReadsHtsgetDataSource(FIRST_TEST_BAM);
        source.close();
        source.iterator();
    }

    @DataProvider(name = "SingleFileCompleteTraversalData")
    public Object[][] getSingleFileCompleteTraversalData() {
        // Files, with expected read names in the expected order
        return new Object[][]{
            {FIRST_TEST_BAM, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")},
            {SECOND_TEST_BAM, Arrays.asList("l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v")},
            {THIRD_TEST_BAM, Arrays.asList("w", "x", "y", "z")}
        };
    }

    @Test(dataProvider = "SingleFileCompleteTraversalData")
    public void testSingleFileCompleteTraversal(final GATKPath samFile, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsHtsgetDataSource(samFile)) {
            traverseOnce(readsSource, samFile, expectedReadNames);
        }
    }

    @Test(dataProvider = "SingleFileCompleteTraversalData")
    public void testSingleFileSerialTraversal(final GATKPath samFile, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsHtsgetDataSource(samFile)) {
            Assert.assertTrue(readsSource.supportsSerialIteration());

            traverseOnce(readsSource, samFile, expectedReadNames);
            traverseOnce(readsSource, samFile, expectedReadNames);
            traverseOnce(readsSource, samFile, expectedReadNames);
        }
    }

    private static void traverseOnce(final ReadsDataSource readsSource, final GATKPath samFile, final List<String> expectedReadNames) {
        final List<GATKRead> reads = new ArrayList<>();
        readsSource.forEach(reads::add);

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in complete traversal of " + samFile);

        // Make sure we got the reads we expected in the right order
        for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
            Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in complete traversal of " + samFile + " not equal to expected read");
        }
    }

    @DataProvider(name = "SingleFileTraversalWithIntervalsData")
    public Object[][] getSingleFileTraversalWithIntervalsData() {
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
                Collections.<String>emptyList()
            },
            {FIRST_TEST_BAM,
                Arrays.asList(new SimpleInterval("1", 200, 205), new SimpleInterval("1", 207, 210)),
                Arrays.asList("a", "b", "c")
            },
            {FIRST_TEST_BAM,
                Arrays.asList(new SimpleInterval("2", 500, 502), new SimpleInterval("2", 552, 700)),
                Arrays.asList("f", "g", "h")
            }
        };
    }

    @Test(dataProvider = "SingleFileTraversalWithIntervalsData")
    public void testSingleFileTraversalWithIntervals(final GATKPath samFile, final List<SimpleInterval> intervals, final List<String> expectedReadNames) {
        try (final ReadsHtsgetDataSource readsSource = new ReadsHtsgetDataSource(samFile)) {
            readsSource.setTraversalBounds(intervals);

            final List<GATKRead> reads = new ArrayList<>();
            readsSource.forEach(reads::add);

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFile);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFile + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "SingleFileQueryByIntervalData")
    public Object[][] getSingleFileQueryByIntervalData() {
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
                Collections.<String>emptyList()
            },
            {FIRST_TEST_BAM,
                new SimpleInterval("4", 600, 699),
                Collections.<String>emptyList()
            }
        };
    }

    @Test(dataProvider = "SingleFileQueryByIntervalData")
    public void testSingleFileQueryByInterval(final GATKPath samFile, final SimpleInterval interval, final List<String> expectedReadNames) {
        try (final ReadsHtsgetDataSource readsSource = new ReadsHtsgetDataSource(samFile)) {
            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
        }
    }

    @Test(dataProvider = "SingleFileQueryByIntervalData")
    public void testSingleFileQueryByIntervalSerialIteration(final GATKPath samFile, final SimpleInterval interval, final List<String> expectedReadNames) {
        try (final ReadsHtsgetDataSource readsSource = new ReadsHtsgetDataSource(samFile)) {
            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
            traverseOnceByInterval(readsSource, samFile, interval, expectedReadNames);
        }
    }

    private static void traverseOnceByInterval(final ReadsDataSource readsSource, final GATKPath samFile, final SimpleInterval interval, final List<String> expectedReadNames) {
        final List<GATKRead> reads = new ArrayList<>();
        final Iterator<GATKRead> queryIterator = readsSource.query(interval);
        queryIterator.forEachRemaining(reads::add);

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in query by interval of " + samFile);

        // Make sure we got the reads we expected in the right order
        for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
            Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in query by interval of " + samFile + " not equal to expected read");
        }
    }

    @DataProvider(name = "MultipleFilesCompleteTraversalData")
    public Object[][] getMultipleFilesCompleteTraversalData() {
        // Files, with expected read names in the expected order
        return new Object[][]{
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM), Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "t", "i", "j", "u", "v", "k")},
            {Arrays.asList(SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.asList("l", "m", "n", "o", "p", "q", "r", "s", "w", "t", "x", "u", "v", "y", "z")},
            {Arrays.asList(FIRST_TEST_BAM, THIRD_TEST_BAM), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "w", "x", "i", "j", "y", "k", "z")},
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")}
        };
    }

    @Test(dataProvider = "MultipleFilesCompleteTraversalData")
    public void testMultipleFilesCompleteTraversal(final List<GATKPath> samFiles, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsHtsgetDataSource(samFiles)) {
            final List<GATKRead> reads = new ArrayList<>();
            readsSource.forEach(reads::add);

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
                Collections.<String>emptyList()
            },
            {Arrays.asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                Arrays.asList(new SimpleInterval("1", 1, 16000), new SimpleInterval("2", 1, 16000), new SimpleInterval("3", 1, 16000), new SimpleInterval("4", 1, 16000)),
                Arrays.asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")
            }
        };
    }

    @Test(dataProvider = "MultipleFilesTraversalWithIntervalsData")
    public void testMultipleFilesTraversalWithIntervals(final List<GATKPath> samFiles, final List<SimpleInterval> intervals, final List<String> expectedReadNames) {
        try (final ReadsHtsgetDataSource readsSource = new ReadsHtsgetDataSource(samFiles)) {
            readsSource.setTraversalBounds(intervals);
            final List<GATKRead> reads = new ArrayList<>();
            readsSource.forEach(reads::add);

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
                Collections.<String>emptyList()
            }
        };
    }

    @Test(dataProvider = "MultipleFilesQueryByIntervalData")
    public void testMultipleFilesQueryByInterval(final List<GATKPath> samFiles, final SimpleInterval interval, final List<String> expectedReadNames) {
        try (final ReadsHtsgetDataSource readsSource = new ReadsHtsgetDataSource(samFiles)) {
            final List<GATKRead> reads = new ArrayList<>();
            final Iterator<GATKRead> queryIterator = readsSource.query(interval);
            queryIterator.forEachRemaining(reads::add);

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
        final GATKPath mappedBam = FIRST_TEST_BAM;
        // This bam has mapped reads from various contigs, plus a few unmapped reads with no mapped mate
        final GATKPath unmappedBam = FIRST_TEST_UNMAPPED_BAM;

        // This is a snippet of the CEUTrio.HiSeq.WGS.b37.NA12878 bam from large, with mapped reads
        // from chromosome 20 (with one mapped read having an unmapped mate), plus several unmapped
        // reads with no mapped mate.
        final GATKPath ceuSnippet = CEU_SNIPPET;

        return new Object[][]{
            // One interval, no unmapped
            {unmappedBam, Collections.singletonList(new SimpleInterval("1", 200, 1000)), false, Arrays.asList("a", "b", "c", "d")},
            // One interval, with unmapped
            {unmappedBam, Collections.singletonList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d", "u1", "u2", "u3", "u4", "u5")},
            // Multiple intervals, no unmapped
            {unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), false, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k")},
            // Multiple intervals, with unmapped
            {unmappedBam, Arrays.asList(new SimpleInterval("1", 200, 1000), new SimpleInterval("2", 500, 700), new SimpleInterval("4", 700, 701)), true, Arrays.asList("a", "b", "c", "d", "f", "g", "h", "k", "u1", "u2", "u3", "u4", "u5")},
            // Interval with no overlapping reads, no unmapped
            {unmappedBam, Collections.singletonList(new SimpleInterval("1", 3000, 4000)), false, Collections.<String>emptyList()},
            // Interval with no overlapping reads, with unmapped
            {unmappedBam, Collections.singletonList(new SimpleInterval("1", 3000, 4000)), true, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            // Interval with no overlapping reads, with unmapped, but no unmapped reads in bam
            {mappedBam, Collections.singletonList(new SimpleInterval("1", 3000, 4000)), true, Collections.<String>emptyList()},
            // Interval with overlapping reads, with unmapped, but no unmapped reads in bam
            {mappedBam, Collections.singletonList(new SimpleInterval("1", 200, 1000)), true, Arrays.asList("a", "b", "c", "d")},
            // Null intervals, with unmapped
            {unmappedBam, null, true, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            // Empty intervals, with unmapped
            {unmappedBam, Collections.<SimpleInterval>emptyList(), true, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            // Null intervals, with unmapped, but no unmapped reads in bam
            {mappedBam, null, true, Collections.<String>emptyList()},
            // Empty intervals, with unmapped, but no unmapped reads in bam
            {mappedBam, Collections.<SimpleInterval>emptyList(), true, Collections.<String>emptyList()},
            // Null intervals, no unmapped (an unbounded traversal, so we expect all the reads)
            {unmappedBam, null, false, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "u1", "u2", "u3", "u4", "u5")},
            // Empty intervals, no unmapped (an unbounded traversal, so we expect all the reads)
            {unmappedBam, Collections.<SimpleInterval>emptyList(), false, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "u1", "u2", "u3", "u4", "u5")},
            // Interval containing mapped read with unmapped mate, no unmapped
            {ceuSnippet, Collections.singletonList(new SimpleInterval("20", 10000011, 10000013)), false, Arrays.asList("a", "b", "c", "d", "e", "f", "f")},
            // Interval containing mapped read with unmapped mate, with unmapped
            {ceuSnippet, Collections.singletonList(new SimpleInterval("20", 10000011, 10000013)), true, Arrays.asList("a", "b", "c", "d", "e", "f", "f", "g", "h", "h", "i", "i")},
            // Interval not containing mapped read with unmapped mate, no unmapped
            {ceuSnippet, Collections.singletonList(new SimpleInterval("20", 10000009, 10000011)), false, Arrays.asList("a", "b", "c", "d", "e")},
            // Interval not containing mapped read with unmapped mate, with unmapped
            {ceuSnippet, Collections.singletonList(new SimpleInterval("20", 10000009, 10000011)), true, Arrays.asList("a", "b", "c", "d", "e", "g", "h", "h", "i", "i")}
        };
    }

    @Test(dataProvider = "TraversalWithUnmappedReadsTestData")
    public void testTraversalWithUnmappedReads(final GATKPath samFile, final List<SimpleInterval> queryIntervals, final boolean queryUnmapped, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsHtsgetDataSource(samFile)) {
            readsSource.setTraversalBounds(queryIntervals, queryUnmapped);

            final boolean traversalShouldBeBounded = (queryIntervals != null && !queryIntervals.isEmpty()) || queryUnmapped;
            Assert.assertEquals(traversalShouldBeBounded, readsSource.traversalIsBounded());

            final List<GATKRead> reads = new ArrayList<>();
            readsSource.forEach(reads::add);

            reads.forEach(System.out::println);

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFile);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFile + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "QueryUnmappedTestData")
    public Object[][] queryUnmappedTestData() {
        return new Object[][]{
            {FIRST_TEST_UNMAPPED_BAM, Arrays.asList("u1", "u2", "u3", "u4", "u5")},
            {CEU_SNIPPET, Arrays.asList("g", "h", "h", "i", "i")}
        };
    }

    @Test(dataProvider = "QueryUnmappedTestData")
    public void testQueryUnmapped(final GATKPath samFile, final List<String> expectedReadNames) {
        try (final ReadsDataSource readsSource = new ReadsHtsgetDataSource(samFile)) {
            final List<GATKRead> reads = new ArrayList<>();
            final Iterator<GATKRead> queryIterator = readsSource.queryUnmapped();
            queryIterator.forEachRemaining(reads::add);

            // Make sure we got the right number of reads
            Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in queryUnmapped on " + samFile);

            // Make sure we got the reads we expected in the right order
            for (int readIndex = 0; readIndex < reads.size(); ++readIndex) {
                Assert.assertEquals(reads.get(readIndex).getName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in queryUnmapped on " + samFile + " not equal to expected read");
            }
        }
    }

    @DataProvider(name = "MergedHeaderIntervalQueries")
    public Object[][] mergedHeaderQueries() {
        return new Object[][]{
            // Single interval that only overlaps one of the input sources
            {
                Collections.singletonList(new SimpleInterval("EXTRA_CONTIG_1", 1, 10)),
                Collections.singletonList("EXTRA_CONTIG_1_READ"),
                Collections.singletonList(4)
            },
            {
                Collections.singletonList(new SimpleInterval("EXTRA_CONTIG_2", 1, 10)),
                Collections.singletonList("EXTRA_CONTIG_2_READ"),
                Collections.singletonList(5)
            },
            // Multiple intervals, each of which only overlaps one of the input sources
            {
                Arrays.asList(new SimpleInterval("EXTRA_CONTIG_1", 1, 10), new SimpleInterval("EXTRA_CONTIG_2", 1, 10)),
                Arrays.asList("EXTRA_CONTIG_1_READ", "EXTRA_CONTIG_2_READ"),
                Arrays.asList(4, 5)
            },
            // Single interval doesn't overlap ANY input source
            {
                Collections.singletonList(new SimpleInterval("TOTALLY_FAKE_CONTIG", 1, 10)),
                Collections.emptyList(),
                Collections.emptyList()
            }
        };
    }

    @Test(dataProvider = "MergedHeaderIntervalQueries")
    public void testMergedQueryWithFileSpecificContigs(
        final List<SimpleInterval> intervals,
        final List<String> expectedReadNames,
        final List<Integer> expectedSequenceIndex) {
        // create two files, each with a read referencing a sequence that is not present in the other
        final File testFile1 = getFileWithAddedContig(FIRST_TEST_BAM, "EXTRA_CONTIG_1", "test1", ".bam");
        final File testFile2 = getFileWithAddedContig(SECOND_TEST_BAM, "EXTRA_CONTIG_2", "test2", ".bam");

        final GATKPath testPath1 = new GATKPath(ENDPOINT + LOCAL_PREFIX + testFile1.getName());
        final GATKPath testPath2 = new GATKPath(ENDPOINT + LOCAL_PREFIX + testFile2.getName());

        try (final ReadsDataSource readsSource = new ReadsHtsgetDataSource(Arrays.asList(testPath1, testPath2))) {
            final SAMFileHeader samHeader = readsSource.getHeader();
            final SAMSequenceDictionary sequenceDictionary = readsSource.getSequenceDictionary();

            // 4 sequences in the original dictionaries, plus 1 added to each of the 2 inputs == 6
            Assert.assertEquals(sequenceDictionary.getSequences().size(), 6);
            Assert.assertEquals(samHeader.getSequenceDictionary(), sequenceDictionary);

            readsSource.setTraversalBounds(intervals);
            int count = 0;
            for (final GATKRead read : readsSource) {
                Assert.assertEquals(expectedReadNames.get(count), read.getName());
                Assert.assertEquals(
                    ReadUtils.getReferenceIndex(read, samHeader),
                    expectedSequenceIndex.get(count).intValue()
                );
                count++;
            }
            Assert.assertEquals(count, expectedReadNames.size());
        }
    }

    // Copy the reads in the inputFile to the output file, adding an extra contig to the sequence dictionary
    // and a read referencing that contig
    private static File getFileWithAddedContig(
        final GATKPath inputPath,
        final String extraContig,
        final String outputName,
        final String extension) {
        final File outputFile = IOUtils.createTempFileInDirectory(outputName, extension, new File(READS_DATA_SOURCE_TEST_DIRECTORY));
        try (final ReadsDataSource readsSource = new ReadsHtsgetDataSource(inputPath)) {
            final SAMFileHeader header = readsSource.getHeader();
            final SAMSequenceRecord fakeSequenceRec = new SAMSequenceRecord(extraContig, 100);
            header.addSequence(fakeSequenceRec);

            try (final SAMFileGATKReadWriter gatkReadWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(
                outputFile, null, header, true, true, false))) {

                readsSource.forEach(gatkReadWriter::addRead);

                // use the contig name in the read name to make it easy to see where this read came from
                final SAMRecord samRec = new SAMRecord(header);
                samRec.setReadName(extraContig + "_READ");
                samRec.setReferenceName(extraContig);
                samRec.setAlignmentStart(5);
                samRec.setReadBases(new byte[]{'a', 'c', 'g', 't'});
                gatkReadWriter.addRead(new SAMRecordToGATKReadAdapter(samRec));
            }
        }
        return outputFile;
    }
}

// dataproctestutils
