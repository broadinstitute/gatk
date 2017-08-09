package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

public class ReadCachingIterator implements Iterator<GATKRead> {
    
    private final Iterator<GATKRead> wrappedIter;
    private List<GATKRead> cache;
    private static final int INITIAL_CACHE_CAPACITY = 10000;

    public ReadCachingIterator(final Iterator<GATKRead> wrappedIter) {
        this.wrappedIter = wrappedIter;
        this.cache = new ArrayList<>(INITIAL_CACHE_CAPACITY);
    }

    @Override
    public boolean hasNext() {
        return wrappedIter.hasNext();
    }

    @Override
    public GATKRead next() {
        if ( ! hasNext() ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        final GATKRead nextRead = wrappedIter.next();
        cache.add(nextRead);
        return nextRead;
    }

    public List<GATKRead> consumeCachedReads() {
        final List<GATKRead> oldCache = cache;
        cache = new ArrayList<>(INITIAL_CACHE_CAPACITY);
        return oldCache;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported");
    }
}
