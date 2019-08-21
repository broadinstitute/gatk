package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A MultiplePassReadWalker traverses input reads multiple times. To use this class, implement the
 * method @{code traverseReads()}, calling {@link #forEachRead(GATKReadConsumer)} as many times as required.
 * The functional object provided to forEachRead will be called with each read.
 */
public abstract class MultiplePassReadWalker extends ReadWalker {
    private CountingReadFilter countedFilter;
    private boolean firstPass = true;

    @FunctionalInterface
    public interface GATKReadConsumer {
        void consume( GATKRead read, ReferenceContext reference, FeatureContext features );
    }

    public abstract void traverseReads();

    /**
     * MultiplePassReadWalker-derived tools should call {@code forEachRead} once for each requested traversal
     * of the input reads. Reads will be processed and the {@code readHandler} will be presented with each
     * {@code GATKRead} processed during the traversal.
     * @param readHandler handler for reads for the current reads iteration
     */
    public void forEachRead( final GATKReadConsumer readHandler) {
        if (!firstPass) {
            resetReadsDataSource();
        }

        getTransformedReadStream(countedFilter).forEach( read -> {
            final SimpleInterval readInterval = getReadInterval(read);
            readHandler.consume(
                    read,
                    new ReferenceContext(reference, readInterval), // will be empty if reference or readInterval is null
                    new FeatureContext(features, readInterval));   // will be empty if features or readInterval is null
            progressMeter.update(readInterval);
        });

        firstPass = false;
    }

    @Override
    public final void traverse() {
        countedFilter = makeReadFilter();
        traverseReads();
        logger.info(countedFilter.getSummaryLine());
    }

    @Override
    public final void apply( final GATKRead read,
                             final ReferenceContext referenceContext,
                             final FeatureContext featureContext ) {
        throw new GATKException("apply can't be called in MultiplePassReadWalker");
    }
}
