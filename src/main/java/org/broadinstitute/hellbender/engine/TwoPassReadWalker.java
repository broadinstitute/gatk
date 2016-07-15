package org.broadinstitute.hellbender.engine;


import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.walkers.rnaseq.SplitNCigarReads;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.stream.StreamSupport;


/**
 * A ReadWalker is a tool that processes a single read at a time from one or multiple sources of reads, with
 * optional contextual information from a reference and/or sets of variants/Features. The TwoPassReadsWalker iterates
 * through the reads twice, allowing for different processing to be performed on the reads each time.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * TwoPassReadWalker authors must implement the firstPassApply() and secondPassApply() methods to process each read,
 * and may optionally implement onTraversalStart() and/or onTraversalSuccess() and afterFirstPass() to perform
 * operations between passes. See the {@link SplitNCigarReads} walker for an example.
 */
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

    /**
     * Process an individual read (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible. Indicates what actions are to be taken on the first iteration through the reads.
     *
     * @param read current read
     * @param referenceContext Reference bases spanning the current read. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current read's interval
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current read. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    abstract protected void firstPassApply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext);

    /**
     * Similar to firstPassApply(), will perform processing on the reads that is intended for the second iteration
     * through the reads.
     *
     * @param read
     * @param referenceContext
     * @param featureContext
     */
    abstract protected void secondPassApply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext);

    /**
     * A method to be overridden. This gets called between iteration through the reads. Allows for tool authors to change
     * the state of the tool between the first and second pass.
     */
    protected void afterFirstPass() {
    }


    @Override
    final public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {}
}
