package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.List;

import static org.mockito.Mockito.*;

public class ReadCachingIteratorUnitTest extends GATKBaseTest {

    @Test
    public void testReadCachingIterator() {
        @SuppressWarnings("unchecked")
        Iterator<GATKRead> mockIterator = (Iterator<GATKRead>)mock(Iterator.class);
        when(mockIterator.hasNext()).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(false);
        when(mockIterator.next()).thenReturn(ArtificialReadUtils.createArtificialRead("101M"));

        final ReadCachingIterator cachingIterator = new ReadCachingIterator(mockIterator);
        for ( int i = 1; i <= 3; ++i ) {
            cachingIterator.next();
        }
        List<GATKRead> cachedReads = cachingIterator.consumeCachedReads();
        Assert.assertEquals(cachedReads.size(), 3, "wrong number of reads cached");

        cachedReads = cachingIterator.consumeCachedReads();
        Assert.assertEquals(cachedReads.size(), 0, "wrong number of reads cached");

        for ( int i = 1; i <= 5; ++i ) {
            cachingIterator.next();
        }
        cachedReads = cachingIterator.consumeCachedReads();
        Assert.assertEquals(cachedReads.size(), 5, "wrong number of reads cached");

        cachedReads = cachingIterator.consumeCachedReads();
        Assert.assertEquals(cachedReads.size(), 0, "wrong number of reads cached");

        for ( int i = 1; i <= 2; ++i ) {
            cachingIterator.next();
        }
        cachedReads = cachingIterator.consumeCachedReads();
        Assert.assertEquals(cachedReads.size(), 2, "wrong number of reads cached");

        cachedReads = cachingIterator.consumeCachedReads();
        Assert.assertEquals(cachedReads.size(), 0, "wrong number of reads cached");

        Assert.assertFalse(cachingIterator.hasNext(), "iterator should be exhausted but isn't");
    }
}
