package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.mockito.Mockito.*;

public class AllLocusIteratorUnitTest extends GATKBaseTest {

    @DataProvider
    public Object[][] testIteratorReturnsPileupsForAllLociData() {
        // All provided pileups within the interval
        final SimpleInterval firstInterval = new SimpleInterval("1", 1, 25);
        final List<AlignmentContext> firstPileups = Arrays.asList(
                new AlignmentContext(new SimpleInterval("1", 5, 5),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()),
                new AlignmentContext(new SimpleInterval("1", 15, 15),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()),
                new AlignmentContext(new SimpleInterval("1", 19, 19),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()),
                new AlignmentContext(new SimpleInterval("1", 22, 22),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()));
        final List<Integer> firstLoci = Arrays.asList(5, 15, 19, 22);

        // Some provided pileups are before/after the interval
        final SimpleInterval secondInterval = new SimpleInterval("1", 10, 20);
        final List<AlignmentContext> secondPileups = Arrays.asList(
                new AlignmentContext(new SimpleInterval("1", 5, 5),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()),
                new AlignmentContext(new SimpleInterval("1", 15, 15),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()),
                new AlignmentContext(new SimpleInterval("1", 19, 19),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()),
                new AlignmentContext(new SimpleInterval("1", 22, 22),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()));
        final List<Integer> secondLoci = Arrays.asList(5, 15, 19, 22);

        // No provided pileups
        final SimpleInterval thirdInterval = new SimpleInterval("1", 10, 25);
        final List<AlignmentContext> thirdPileups = Collections.emptyList();
        final List<Integer> thirdLoci = Collections.emptyList();

        return new Object[][] {
                { firstInterval, firstPileups, firstLoci },
                { secondInterval, secondPileups, secondLoci },
                { thirdInterval, thirdPileups, thirdLoci }
        };
    }

    @Test(dataProvider = "testIteratorReturnsPileupsForAllLociData")
    public void testIteratorReturnsPileupsForAllLoci( final SimpleInterval providedInterval, final List<AlignmentContext> providedPileups, final List<Integer> providedLoci ) {
        final Iterator<AlignmentContext> locusIterator = new AllLocusIterator(providedInterval, providedPileups.iterator());

        final List<AlignmentContext> returnedPileups = new ArrayList<>();
        locusIterator.forEachRemaining(pileup -> returnedPileups.add(pileup));

        Assert.assertEquals(returnedPileups.size(), providedInterval.size(), "Wrong number of pileup objects returned");

        AlignmentContext previousPileup = null;
        int position = providedInterval.getStart();

        for ( final AlignmentContext pileup : returnedPileups ) {
            Assert.assertEquals(pileup.getContig(), providedInterval.getContig(), "Wrong contig in returned pileup");
            Assert.assertEquals(pileup.getStart(), position, "Wrong position in returned pileup");

            if ( providedLoci.contains(pileup.getStart()) ) {
                Assert.assertEquals(pileup.getBasePileup().size(), 10, "Didn't get back a pileup we provided explicitly");
            } else {
                Assert.assertEquals(pileup.getBasePileup().size(), 0, "Didn't get back an empty pileup where expected");
            }
            
            Assert.assertEquals(pileup.getEnd() - pileup.getStart() + 1, 1, "Returned pileup has wrong size");

            if ( previousPileup != null ) {
                Assert.assertEquals(previousPileup.getContig(), pileup.getContig(), "Returned pileup has wrong contig");
                Assert.assertEquals(previousPileup.getStart(), pileup.getStart() - 1, "Returned pileup has wrong location");
            }

            previousPileup = pileup;
            position++;
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testWrongContig() {
        final SimpleInterval interval = new SimpleInterval("1", 1, 25);
        final List<AlignmentContext> badPileups = Arrays.asList(
                new AlignmentContext(new SimpleInterval("2", 5, 5),
                        when(mock(ReadPileup.class).size()).thenReturn(10).getMock()));

        final Iterator<AlignmentContext> locusIterator = new AllLocusIterator(interval, badPileups.iterator());
    }
}
