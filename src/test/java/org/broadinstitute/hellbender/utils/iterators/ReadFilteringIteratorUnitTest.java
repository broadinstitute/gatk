package org.broadinstitute.hellbender.utils.iterators;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Iterator;

public class ReadFilteringIteratorUnitTest extends GATKBaseTest {

    private Iterator<GATKRead> makeReadsIterator() {
        return Arrays.asList(
                ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("101M")),
                ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("50M")),
                ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("101M")),
                ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("30M")),
                ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode("101M"))
        ).iterator();
    }

    @DataProvider(name = "FilteringIteratorTestData")
    public Object[][] filteringIteratorTestData() {
        final ReadFilter allowNoReadsFilter = new ReadFilter() {
            private static final long serialVersionUID = 1L;
            @Override public boolean test(final GATKRead read){return false;} };
        final ReadFilter allowLongReadsFilter = new ReadFilter() {
            private static final long serialVersionUID = 1L;
            @Override public boolean test(final GATKRead read){return read.getLength() > 100;} };

        return new Object[][] {
                { makeReadsIterator(), ReadFilterLibrary.ALLOW_ALL_READS, 5 },
                { makeReadsIterator(), allowNoReadsFilter, 0 },
                { makeReadsIterator(), allowLongReadsFilter, 3 }
        };
    }

    @Test(dataProvider = "FilteringIteratorTestData")
    public void testFilteringIterator( final Iterator<GATKRead> readsIter, final ReadFilter filter, final int expectedReadCount ) {
        final ReadFilteringIterator filteringIterator = new ReadFilteringIterator(readsIter, filter);
        int count = 0;
        for ( final GATKRead read : filteringIterator ) {
            ++count;
        }
        Assert.assertEquals(count, expectedReadCount, "Wrong number of reads from the ReadFilteringIterator");
    }
}
