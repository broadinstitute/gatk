package org.broadinstitute.hellbender.utils.downsampling;

import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * An extension of the basic downsampler API with reads-specific operations
 */
public abstract class ReadsDownsampler extends Downsampler<GATKRead> {

    /**
     * Does this downsampler require that reads be fed to it in coordinate order?
     *
     * @return true if reads must be submitted to this downsampler in coordinate order, otherwise false
     */
    public abstract boolean requiresCoordinateSortOrder();

    /**
     * Tell this downsampler that no more reads located before the provided read (according to
     * the sort order of the read stream) will be fed to it.
     *
     * Allows position-aware downsamplers to finalize pending reads earlier than they would
     * otherwise be able to, particularly when doing per-sample downsampling and reads for
     * certain samples are sparser than average.
     *
     * @param read the downsampler will assume that no reads located before this read will ever
     *             be submitted to it in the future
     */
    public abstract void signalNoMoreReadsBefore( final GATKRead read );
}
