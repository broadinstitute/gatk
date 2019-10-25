package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A MultiplePassReadWalker traverses input reads multiple times. To use this class, implement the
 * method {@link #traverseReads()}, calling {@link #forEachRead(GATKReadConsumer)} with a
 * {@link GATKReadConsumer} implementation once for each desired traversal of the input reads. The
 * {@link GATKReadConsumer} provided to {@link #forEachRead(GATKReadConsumer) will be called presented
 * with each read.
 */
public abstract class MultiplePassReadWalker extends ReadWalker {
    private CountingReadFilter countedFilter;
    private int passCount = 1;

    /**
     * Implemented by MultiplePassReadWalker-derived tools. The implementation should call {@code forEachRead}
     * once for each desired traversal of the input reads, passing a {@code GATKReadConsumer} for each call.
     * Input reads will be processed and the {@code GATKReadConsumer} will be presented with each {@code GATKRead}
     * processed during the traversal.
     */
    public abstract void traverseReads();

    /**
     * Implemented by MultiplePassReadWalker-derived tools. Tool implementers should provide an implementation
     * of this {@link FunctionalInterface} for each call to {@link #forEachRead(GATKReadConsumer)}. The provided
     * {@link FunctionalInterface} will be called once with read processed during the iteration. For more
     * information about the arguments supplied to this method,
     * @see ReadWalker#apply(GATKRead, ReferenceContext, FeatureContext)
     */
    @FunctionalInterface
    public interface GATKReadConsumer {
        void consume( GATKRead read, ReferenceContext reference, FeatureContext features );
    }

    /**
     * Provides a single traversal over all input reads. MultiplePassReadWalker-derived tools should call this
     * method from the {@link #traverseReads} implementation once for each requested traversal of the input reads,
     * providing the {@code readHandler}, which will be presented with each {@code GATKRead} processed during
     * the traversal.
     *
     * NOTE: Tool implementers should consider whether any state that might otherwise be retained across
     * traversals, such as that maintained by a random number generator or downsampler, should be reset between
     * traversals in order to ensure that each traversal is executed under the same conditions.
     *
     * @param readHandler handler for reads for the current reads iteration
     */
    public void forEachRead( final GATKReadConsumer readHandler) {
        if (passCount > 1) {
            countedFilter = makeReadFilter();
            resetReadsDataSource();
            logger.info(String.format("Starting traversal pass %d", passCount));
        }

        getTransformedReadStream(countedFilter).forEach( read -> {
            final SimpleInterval readInterval = getReadInterval(read);
            readHandler.consume(
                    read,
                    new ReferenceContext(reference, readInterval), // will be empty if reference or readInterval is null
                    new FeatureContext(features, readInterval));   // will be empty if features or readInterval is null
            progressMeter.update(readInterval);
        });

        logger.info(countedFilter.getSummaryLine());
        passCount++;
    }

    @Override
    public final void traverse() {
        countedFilter = makeReadFilter();
        traverseReads();
    }

    @Override
    public final void apply( final GATKRead read,
                             final ReferenceContext referenceContext,
                             final FeatureContext featureContext ) {
        throw new GATKException("apply can't be called in MultiplePassReadWalker");
    }
}
