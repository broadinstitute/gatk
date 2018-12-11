package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.iterators.PushToPullIterator;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;

/**
 * Iterator wrapper around our generic {@link ReadsDownsampler)} interface. Wraps an iterator of reads,
 * and downsamples the reads from that iterator using the provided downsampler.
 *
 * Converts the push-style {@link ReadsDownsampler)} interface to a pull model.
 */
public final class ReadsDownsamplingIterator extends PushToPullIterator<GATKRead> {

    /**
     * @param iter        wrapped iterator from which this iterator will pull reads to be downsampled
     * @param downsampler downsampler through which the reads from the wrapped iterator will be fed
     */
    public ReadsDownsamplingIterator(Iterator<GATKRead> iter, ReadsDownsampler downsampler) {
        super(iter, downsampler);
    }
}