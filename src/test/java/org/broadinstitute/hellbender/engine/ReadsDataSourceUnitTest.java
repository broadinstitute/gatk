package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public final class ReadsDataSourceUnitTest extends BaseTest {
    private static final String READS_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final File FIRST_TEST_BAM = new File(READS_DATA_SOURCE_TEST_DIRECTORY + "reads_data_source_test1.bam");
    private static final File SECOND_TEST_BAM = new File(READS_DATA_SOURCE_TEST_DIRECTORY + "reads_data_source_test2.bam");
    private static final File THIRD_TEST_BAM = new File(READS_DATA_SOURCE_TEST_DIRECTORY + "reads_data_source_test3.bam");
    private static final File FIRST_TEST_SAM = new File(READS_DATA_SOURCE_TEST_DIRECTORY + "invalid_coord_sort_order.sam");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFile() {
        File nullFile = null;
        ReadsDataSource readsSource = new ReadsDataSource(nullFile);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFileList() {
        List<File> nullList = null;
        ReadsDataSource readsSource = new ReadsDataSource(nullList);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleEmptyFileList() {
        ReadsDataSource readsSource = new ReadsDataSource(new ArrayList<>());
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testHandleNonExistentFile() {
        ReadsDataSource readsSource = new ReadsDataSource(new File("/foo/bar/nonexistent.bam"));
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleUnindexedFileWithIntervals() {
        // Cannot initialize a reads source with intervals unless all files are indexed
        ReadsDataSource readsSource = new ReadsDataSource(new File(READS_DATA_SOURCE_TEST_DIRECTORY + "unindexed.bam"));
        readsSource.setIntervalsForTraversal(Arrays.asList(new SimpleInterval("1", 1, 5)));
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleUnindexedFileQuery() {
        // Construction should succeed, since we don't pass in any intervals, but the query should throw.
        ReadsDataSource readsSource = new ReadsDataSource(new File(READS_DATA_SOURCE_TEST_DIRECTORY + "unindexed.bam"));
        readsSource.query(new SimpleInterval("1", 1, 5));
    }

    @Test
    public void testDefaultSamReaderValidationStringency() {
        // Default validation stringency = SILENT results in no validation errors on invalid coordinate sort
        final ReadsDataSource readsSource = new ReadsDataSource(FIRST_TEST_SAM);
        //noinspection StatementWithEmptyBody
        for ( @SuppressWarnings("unused") final SAMRecord read : readsSource ) {
        }
    }

    @Test(expectedExceptions = SAMFormatException.class)
    public void testCustomSamReaderFactory() {
        // Custom SamReaderFactory with validation stringency = STRICT fails on invalid coordinate sort
        final ReadsDataSource readsSource = new ReadsDataSource(
                FIRST_TEST_SAM,
                SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT));
        //noinspection StatementWithEmptyBody
        for ( @SuppressWarnings("unused") final SAMRecord read : readsSource ) {
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
    public void testSingleFileCompleteTraversal( final File samFile, final List<String> expectedReadNames ) {
        ReadsDataSource readsSource = new ReadsDataSource(samFile);

        List<SAMRecord> reads = new ArrayList<>();
        for ( SAMRecord read : readsSource ) {
            reads.add(read);
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in complete traversal of " + samFile.getAbsolutePath());

        // Make sure we got the reads we expected in the right order
        for ( int readIndex = 0; readIndex < reads.size(); ++readIndex ) {
            Assert.assertEquals(reads.get(readIndex).getReadName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in complete traversal of " + samFile.getAbsolutePath() + " not equal to expected read");
        }

        readsSource.close();
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
    public void testSingleFileTraversalWithIntervals( final File samFile, final List<SimpleInterval> intervals, final List<String> expectedReadNames ) {
        ReadsDataSource readsSource = new ReadsDataSource(samFile);
        readsSource.setIntervalsForTraversal(intervals);

        List<SAMRecord> reads = new ArrayList<>();
        for ( SAMRecord read : readsSource ) {
            reads.add(read);
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFile.getAbsolutePath());

        // Make sure we got the reads we expected in the right order
        for ( int readIndex = 0; readIndex < reads.size(); ++readIndex ) {
            Assert.assertEquals(reads.get(readIndex).getReadName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFile.getAbsolutePath() + " not equal to expected read");
        }

        readsSource.close();
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
    public void testSingleFileQueryByInterval( final File samFile, final SimpleInterval interval, final List<String> expectedReadNames ) {
        ReadsDataSource readsSource = new ReadsDataSource(samFile);

        List<SAMRecord> reads = new ArrayList<>();
        Iterator<SAMRecord> queryIterator = readsSource.query(interval);
        while ( queryIterator.hasNext() ) {
            reads.add(queryIterator.next());
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in query by interval of " + samFile.getAbsolutePath());

        // Make sure we got the reads we expected in the right order
        for ( int readIndex = 0; readIndex < reads.size(); ++readIndex ) {
            Assert.assertEquals(reads.get(readIndex).getReadName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in query by interval of " + samFile.getAbsolutePath() + " not equal to expected read");
        }

        readsSource.close();
    }

    @DataProvider(name = "MultipleFilesCompleteTraversalData")
    public Object[][] getMultipleFilesCompleteTraversalData() {
        // Files, with expected read names in the expected order
        return new Object[][] {
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM), Arrays.<String>asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "t", "i", "j", "u", "v", "k") },
                { Arrays.<File>asList(SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.<String>asList("l", "m", "n", "o", "p", "q", "r", "s", "w", "t", "x", "u", "v", "y", "z") },
                { Arrays.<File>asList(FIRST_TEST_BAM, THIRD_TEST_BAM), Arrays.<String>asList("a", "b", "c", "d", "e", "f", "g", "h", "w", "x", "i", "j", "y", "k", "z") },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM), Arrays.<String>asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z") }
        };
    }

    @Test(dataProvider = "MultipleFilesCompleteTraversalData")
    public void testMultipleFilesCompleteTraversal(final List<File> samFiles, final List<String> expectedReadNames) {
        ReadsDataSource readsSource = new ReadsDataSource(samFiles);
        List<SAMRecord> reads = new ArrayList<>();

        for ( SAMRecord read : readsSource ) {
            reads.add(read);
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in complete traversal of " + samFiles);

        // Make sure we got the reads we expected in the right order
        for ( int readIndex = 0; readIndex < reads.size(); ++readIndex ) {
            Assert.assertEquals(reads.get(readIndex).getReadName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in complete traversal of " + samFiles + " not equal to expected read");
        }

        readsSource.close();
    }

    @DataProvider(name = "MultipleFilesTraversalWithIntervalsData")
    public Object[][] getMultipleFilesTraversalWithIntervalsData() {
        // Files, with intervals, and expected read names in the expected order
        return new Object[][] {
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 205, 207), new SimpleInterval("1", 400, 1000), new SimpleInterval("4", 500, 704)),
                  Arrays.<String>asList("a", "b", "l", "n", "d", "y", "k")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("4", 500, 704), new SimpleInterval("1", 400, 1000), new SimpleInterval("1", 205, 207)),
                  Arrays.<String>asList("a", "b", "l", "n", "d", "y", "k")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("2", 500, 600), new SimpleInterval("2", 2099, 2200), new SimpleInterval("3", 50, 100), new SimpleInterval("3", 300, 500)),
                  Arrays.<String>asList("f", "p", "g", "s", "w", "i", "j", "u")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 1, 300), new SimpleInterval("1", 100, 500), new SimpleInterval("1", 200, 600)),
                  Arrays.<String>asList("a", "b", "l", "c", "m", "n")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 11000, 12000), new SimpleInterval("3", 1000, 2000)),
                  Arrays.<String>asList()
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  Arrays.<SimpleInterval>asList(new SimpleInterval("1", 1, 16000), new SimpleInterval("2", 1, 16000), new SimpleInterval("3", 1, 16000), new SimpleInterval("4", 1, 16000)),
                  Arrays.<String>asList("a", "b", "l", "c", "m", "n", "d", "e", "o", "f", "p", "g", "h", "q", "r", "s", "w", "t", "x", "i", "j", "u", "v", "y", "k", "z")
                }
        };
    }

    @Test(dataProvider = "MultipleFilesTraversalWithIntervalsData")
    public void testMultipleFilesTraversalWithIntervals( final List<File> samFiles, final List<SimpleInterval> intervals, final List<String> expectedReadNames ) {
        ReadsDataSource readsSource = new ReadsDataSource(samFiles);
        readsSource.setIntervalsForTraversal(intervals);

        List<SAMRecord> reads = new ArrayList<>();
        for ( SAMRecord read : readsSource ) {
            reads.add(read);
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in traversal by intervals of " + samFiles);

        // Make sure we got the reads we expected in the right order
        for ( int readIndex = 0; readIndex < reads.size(); ++readIndex ) {
            Assert.assertEquals(reads.get(readIndex).getReadName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in traversal by intervals of " + samFiles + " not equal to expected read");
        }

        readsSource.close();
    }

    @DataProvider(name = "MultipleFilesQueryByIntervalData")
    public Object[][] getMultipleFilesQueryByIntervalData() {
        // Files, with a single query interval, and expected read names in the expected order
        return new Object[][] {
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("1", 285, 1000),
                  Arrays.<String>asList("c", "m", "n", "d")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("3", 200, 300),
                  Arrays.<String>asList("t", "x", "i")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("1", 9000, 11000),
                  Arrays.<String>asList("o")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("3", 1, 16000),
                  Arrays.<String>asList("w", "t", "x", "i", "j", "u", "v")
                },
                { Arrays.<File>asList(FIRST_TEST_BAM, SECOND_TEST_BAM, THIRD_TEST_BAM),
                  new SimpleInterval("2", 10000, 12000),
                  Arrays.<String>asList()
                }
        };
    }

    @Test(dataProvider = "MultipleFilesQueryByIntervalData")
    public void testMultipleFilesQueryByInterval( final List<File> samFiles, final SimpleInterval interval, final List<String> expectedReadNames ) {
        ReadsDataSource readsSource = new ReadsDataSource(samFiles);

        List<SAMRecord> reads = new ArrayList<>();
        Iterator<SAMRecord> queryIterator = readsSource.query(interval);
        while ( queryIterator.hasNext() ) {
            reads.add(queryIterator.next());
        }

        // Make sure we got the right number of reads
        Assert.assertEquals(reads.size(), expectedReadNames.size(), "Wrong number of reads returned in query by interval of " + samFiles);

        // Make sure we got the reads we expected in the right order
        for ( int readIndex = 0; readIndex < reads.size(); ++readIndex ) {
            Assert.assertEquals(reads.get(readIndex).getReadName(), expectedReadNames.get(readIndex), "Read #" + (readIndex + 1) + " in query by interval of " + samFiles + " not equal to expected read");
        }

        readsSource.close();
    }

}
