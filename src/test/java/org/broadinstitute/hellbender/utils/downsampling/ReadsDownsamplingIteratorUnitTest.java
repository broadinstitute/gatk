package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ReadsDownsamplingIteratorUnitTest extends GATKBaseTest {

    // Toy downsampler that keeps only reads with names "B" or "C"
    private static class KeepReadsBAndCOnlyDownsampler extends ReadsDownsampler {
        private List<GATKRead> finalizedReads = new ArrayList<>();

        @Override
        public boolean requiresCoordinateSortOrder() {
            return false;
        }

        @Override
        public void signalNoMoreReadsBefore( GATKRead read ) {
            // no-op
        }

        @Override
        public void submit( GATKRead item ) {
            if ( item.getName() != null && (item.getName().equals("B") || item.getName().equals("C")) ) {
                finalizedReads.add(item);
            }
            else {
                incrementNumberOfDiscardedItems(1);
            }
        }

        @Override
        public boolean hasFinalizedItems() {
            return ! finalizedReads.isEmpty();
        }

        @Override
        public List<GATKRead> consumeFinalizedItems() {
            final List<GATKRead> toReturn = finalizedReads;
            finalizedReads = new ArrayList<>();
            return toReturn;
        }

        @Override
        public boolean hasPendingItems() {
            return false;
        }

        @Override
        public GATKRead peekFinalized() {
            return hasFinalizedItems() ? finalizedReads.get(0) : null;
        }

        @Override
        public GATKRead peekPending() {
            return null;
        }

        @Override
        public int size() {
            return finalizedReads.size();
        }

        @Override
        public void signalEndOfInput() {
            // no-op
        }

        @Override
        public void clearItems() {
            finalizedReads.clear();
        }
    };

    @DataProvider(name = "ReadsDownsamplingIteratorTestData")
    public Object[][] readsDownsamplingIteratorTestData() {
        return new Object[][] {
                // input reads, and expected read names after downsampling
                { Arrays.asList(readWithName("A")), Collections.emptyList() },
                { Arrays.asList(readWithName("A"), readWithName("B")), Arrays.asList("B") },
                { Arrays.asList(readWithName("A"), readWithName("A"), readWithName("A"), readWithName("B"), readWithName("A"), readWithName("A")), Arrays.asList("B") },
                { Arrays.asList(readWithName("A"), readWithName("A"), readWithName("A"), readWithName("A"), readWithName("A")), Collections.emptyList() },
                { Arrays.asList(readWithName("A"), readWithName("B"), readWithName("C")), Arrays.asList("B", "C") },
                { Arrays.asList(readWithName("B"), readWithName("B"), readWithName("A"), readWithName("B"), readWithName("C")), Arrays.asList("B", "B", "B", "C") },
                { Arrays.asList(readWithName("B"), readWithName("B"), readWithName("C")), Arrays.asList("B", "B", "C") },
                { Collections.emptyList(), Collections.emptyList() }
        };
    }

    @Test(dataProvider = "ReadsDownsamplingIteratorTestData")
    public void testReadsDownsamplingIterator( final List<GATKRead> reads, final List<String> expectedReadNames ) {
        final ReadsDownsamplingIterator downsamplingIter = new ReadsDownsamplingIterator(reads.iterator(), new KeepReadsBAndCOnlyDownsampler());
        final List<GATKRead> readsAfterDownsampling = new ArrayList<>();
        for ( final GATKRead read : downsamplingIter ) {
            readsAfterDownsampling.add(read);
        }

        Assert.assertEquals(readsAfterDownsampling.size(), expectedReadNames.size());
        for ( int i = 0; i < readsAfterDownsampling.size(); ++i ) {
            Assert.assertEquals(readsAfterDownsampling.get(i).getName(), expectedReadNames.get(i));
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullIter() {
        final ReadsDownsamplingIterator iter = new ReadsDownsamplingIterator(null, new KeepReadsBAndCOnlyDownsampler());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullDownsampler() {
        final ReadsDownsamplingIterator iter = new ReadsDownsamplingIterator(Collections.<GATKRead>emptyIterator(), null);
    }

    @Test(expectedExceptions = NoSuchElementException.class)
    public void testNextCalledWithNoItems() {
        final ReadsDownsamplingIterator iter = new ReadsDownsamplingIterator(Collections.<GATKRead>emptyIterator(), new KeepReadsBAndCOnlyDownsampler());
        Assert.assertFalse(iter.hasNext());
        // Should throw a NoSuchElementException
        iter.next();
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void testRemoveThrows() {
        final ReadsDownsamplingIterator iter = new ReadsDownsamplingIterator(Collections.<GATKRead>emptyIterator(), new KeepReadsBAndCOnlyDownsampler());
        iter.remove();
    }

    private GATKRead readWithName( final String name ) {
        return ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("10M"), name);
    }
}
