package org.broadinstitute.hellbender.engine;


import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.stream.StreamSupport;

public abstract class TwoPassReadWalker extends ReadWalker {

    @Override
    public void traverse() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        final CountingReadFilter countedFilter = disable_all_read_filters ?
                new CountingReadFilter("Allow all", ReadFilterLibrary.ALLOW_ALL_READS ) :
                makeReadFilter();

        traverseReads(countedFilter, this::firstPassApply);
        logger.info("Finished First Pass");
        afterFirstPass();
        logger.info("Starting SecondPass");
        traverseReads(countedFilter, this::secondPassApply);

        logger.info(countedFilter.getSummaryLine());
    }

    private void traverseReads(CountingReadFilter countedFilter, GATKApply f) {
        StreamSupport.stream(reads.spliterator(), false)
                .filter(countedFilter)
                .forEach(read -> {
                    final SimpleInterval readInterval = getReadInterval(read);
                    f.consume(read,
                            new ReferenceContext(reference, readInterval), // Will create an empty ReferenceContext if reference or readInterval == null
                            new FeatureContext(features, readInterval));   // Will create an empty FeatureContext if features or readInterval == null

                    progressMeter.update(readInterval);
                });
    }

    @FunctionalInterface
    private interface GATKApply{
        void consume(GATKRead read, ReferenceContext reference, FeatureContext features);
    }

    abstract protected void firstPassApply(GATKRead read, ReferenceContext bytes, FeatureContext featureContext);

    abstract protected void secondPassApply(GATKRead read, ReferenceContext bytes, FeatureContext featureContext);

    protected void afterFirstPass() {
        initializeReads();

        if ( hasIntervals() ) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()));
        }
    }


    @Override
    final public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {}
}
