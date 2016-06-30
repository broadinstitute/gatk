package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Iterator wrapper around our generic {@link ReadsDownsampler)} interface. Wraps an iterator of reads,
 * and downsamples the reads from that iterator using the provided downsampler.
 *
 * Converts the push-style {@link ReadsDownsampler)} interface to a pull model.
 */
public final class ReadsDownsamplingIterator implements Iterator<GATKRead>, Iterable<GATKRead> {

    private final Iterator<GATKRead> nestedReadIterator;
    private final ReadsDownsampler downsampler;
    private Iterator<GATKRead> cachedDownsampledReads = null;
    private GATKRead nextRead = null;

    /**
     * @param iter wrapped iterator from which this iterator will pull reads to be downsampled
     * @param downsampler downsampler through which the reads from the wrapped iterator will be fed
     */
    public ReadsDownsamplingIterator( Iterator<GATKRead> iter, ReadsDownsampler downsampler ) {
        Utils.nonNull(iter, "iterator must not be null");
        Utils.nonNull(downsampler, "downsampler must not be null");

        this.nestedReadIterator = iter;
        this.downsampler = downsampler;

        advanceToNextRead();
    }

    @Override
    public boolean hasNext() {
        return nextRead != null;
    }

    @Override
    public GATKRead next() {
        if ( nextRead == null ) {
            throw new NoSuchElementException("next() called when there are no more items");
        }

        final GATKRead toReturn = nextRead;
        advanceToNextRead();

        return toReturn;
    }

    private void advanceToNextRead() {
        if ( readyToReleaseReads() || fillDownsampledReadsCache() ) {
            nextRead = cachedDownsampledReads.next();
        }
        else {
            nextRead = null;
        }
    }

    private boolean readyToReleaseReads() {
        return cachedDownsampledReads != null && cachedDownsampledReads.hasNext();
    }

    private boolean fillDownsampledReadsCache() {
        while ( nestedReadIterator.hasNext() && ! downsampler.hasFinalizedItems() ) {
            downsampler.submit(nestedReadIterator.next());
        }

        if ( ! nestedReadIterator.hasNext() ) {
            downsampler.signalEndOfInput();
        }

        final Collection<GATKRead> downsampledReads = downsampler.consumeFinalizedItems();
        cachedDownsampledReads = downsampledReads.iterator();

        return cachedDownsampledReads.hasNext();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove records via a ReadsDownsamplingIterator");
    }

    @Override
    public Iterator<GATKRead> iterator() {
        return this;
    }
}