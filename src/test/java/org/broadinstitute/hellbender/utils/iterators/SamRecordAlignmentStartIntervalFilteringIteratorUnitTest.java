package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadQueryIterator;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit test class for {@link SamRecordAlignmentStartIntervalFilteringIterator}
 * Created by jonn on 7/23/19.
 */
public class SamRecordAlignmentStartIntervalFilteringIteratorUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    /**
     * Helper class to deal with converting artificial reads to {@link SAMRecord}s.
     */
    private class ArtificialReadQueryIteratorSamIteratorConverter implements CloseableIterator<SAMRecord> {

        private ArtificialReadQueryIterator iterator;
        private SAMFileHeader header;

        public ArtificialReadQueryIteratorSamIteratorConverter(final SAMFileHeader header,
                                                               final ArtificialReadQueryIterator it) {
            this.header = header;
            iterator = it;
        }

        @Override
        public void close() {}

        @Override
        public SAMRecord next() {
            return iterator.next().convertToSAMRecord(header);
        }

        @Override
        public boolean hasNext() {
            return iterator.hasNext();
        }
    }

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestIteratorEmptyIntervals() {
        return new Object[][] {
                {
                        Collections.emptyList().iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, 5),
                        Collections.emptyList(),
                }
        };
    }

    @DataProvider
    private Object[][] provideForTestIterator() {

        final int numReads = 5;
        final SAMFileHeader header = ArtificialReadUtils.mappedReadIterator(1, 5, numReads).getHeader();

        return new Object[][] {
                // Checks for increasing number of reads on one contig with other contigs containing reads:
                {
                        Collections.singletonList(
                            new SimpleInterval("1", 1, 1)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Collections.singletonList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 1, 2)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 1, 3)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 1, 4)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 1, 5)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                // Check on user interval that extends past end of contig:
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 1, 1000)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                // Checks on interval with at least 1 read starting before it:
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 2, 1000)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 3, 1000)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 4, 1000)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 5, 1000)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Collections.singletonList(
                                new SimpleInterval("1", 6, 1000)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Collections.emptyList(),
                },
                // Multi-interval tests:
                {
                        Arrays.asList(
                                new SimpleInterval("1", 1, 1),
                                new SimpleInterval("1", 2, 2),
                                new SimpleInterval("1", 3, 3),
                                new SimpleInterval("1", 4, 4),
                                new SimpleInterval("1", 5, 5)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Arrays.asList(
                                new SimpleInterval("1", 1, 1),
                                new SimpleInterval("1", 3, 3),
                                new SimpleInterval("1", 4, 4),
                                new SimpleInterval("1", 5, 5)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Arrays.asList(
                                new SimpleInterval("1", 1, 1),
                                new SimpleInterval("1", 2, 2),
                                new SimpleInterval("1", 3, 3),
                                new SimpleInterval("1", 5, 5)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(5), 0, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Arrays.asList(
                                new SimpleInterval("1", 2, 2),
                                new SimpleInterval("1", 4, 4)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                // Multi-Interval / Multi-Contig tests:
                {
                        Arrays.asList(
                                new SimpleInterval("1", 1, 1),
                                new SimpleInterval("2", 2, 2),
                                new SimpleInterval("3", 3, 3),
                                new SimpleInterval("4", 4, 4),
                                new SimpleInterval("5", 5, 5)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads + 2), 1, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 2 + 3), 2, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 3 + 4), 3, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 4 + 5), 4, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Arrays.asList(
                                new SimpleInterval("1", 1, 1),
                                new SimpleInterval("2", 2, 3),
                                new SimpleInterval("3", 3, 3),
                                new SimpleInterval("4", 4, 5),
                                new SimpleInterval("5", 5, 5)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads + 2), 1, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads + 3), 1, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 2 + 3), 2, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 3 + 4), 3, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 3 + 5), 3, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 4 + 5), 4, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
                {
                        Arrays.asList(
                                new SimpleInterval("1", 1, 2),
                                new SimpleInterval("1", 3, 4),
                                new SimpleInterval("3", 3, 3),
                                new SimpleInterval("4", 1, 2),
                                new SimpleInterval("5", 5, 5)
                        ).iterator(),
                        ArtificialReadUtils.mappedReadIterator(1, 5, numReads),
                        Arrays.asList(
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(1), 0, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(2), 0, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(3), 0, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(4), 0, 4, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 2 + 3), 2, 3, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 3 + 1), 3, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 3 + 2), 3, 2, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header),
                                ArtificialReadUtils.createArtificialRead(header, String.valueOf(numReads * 4 + 5), 4, 5, ArtificialReadUtils.DEFAULT_READ_LENGTH).convertToSAMRecord(header)
                        ),
                },
        };
    }

    //==================================================================================================================
    // Tests:

    @Test(dataProvider = "provideForTestIterator")
    public void testIterator( final Iterator<SimpleInterval>    intervalIterator,
                              final ArtificialReadQueryIterator artificialReadQueryIterator,
                              final List<SAMRecord>             expected ) {

        final ArtificialReadQueryIteratorSamIteratorConverter samIterator =
                new ArtificialReadQueryIteratorSamIteratorConverter(
                        artificialReadQueryIterator.getHeader(),
                        artificialReadQueryIterator);

        final SamRecordAlignmentStartIntervalFilteringIterator it =
                new SamRecordAlignmentStartIntervalFilteringIterator(
                        artificialReadQueryIterator.getHeader().getSequenceDictionary(),
                        intervalIterator,
                        samIterator );

        final List<SAMRecord> records = new ArrayList<>();
        while(it.hasNext()) {
            records.add(it.next());
        }

        Assert.assertEquals( records, expected );
    }

    @Test(dataProvider = "provideForTestIteratorEmptyIntervals",
            expectedExceptions = IllegalStateException.class)
    public void testIteratorEmptyIntervals( final Iterator<SimpleInterval>   intervalIterator,
                                            final ArtificialReadQueryIterator artificialReadQueryIterator,
                                            final List<SAMRecord>            expected ) {

        final ArtificialReadQueryIteratorSamIteratorConverter samIterator =
                new ArtificialReadQueryIteratorSamIteratorConverter(
                        artificialReadQueryIterator.getHeader(),
                        artificialReadQueryIterator);

        final SamRecordAlignmentStartIntervalFilteringIterator it =
                new SamRecordAlignmentStartIntervalFilteringIterator(
                        artificialReadQueryIterator.getHeader().getSequenceDictionary(),
                        intervalIterator,
                        samIterator );
    }
}


