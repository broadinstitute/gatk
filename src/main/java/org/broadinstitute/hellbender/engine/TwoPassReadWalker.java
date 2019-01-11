package org.broadinstitute.hellbender.engine;


import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.tools.walkers.rnaseq.SplitNCigarReads;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;
import java.util.stream.StreamSupport;


/**
 * A ReadWalker is a tool that processes a single read at a time from one or multiple sources of reads, with
 * optional contextual information from a reference and/or sets of variants/Features. The TwoPassReadsWalker iterates
 * through the reads twice, allowing for different processing to be performed on the reads each time.
 *
 * WARNING: This traversal should only be used in the rare case that complex state must be maintained between two passes
 * over the reads. In most cases it's preferable to either split the two passes into separate walkers or to restructure
 * the algorithm so that it's not necessary to make multiple passes.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * TwoPassReadWalker authors must implement the {@link #firstPassApply} and {@link #secondPassApply} methods to process
 * each read.  These are analogous to and replace {@link ReadWalker#apply}.  Authors may optionally implement
 * {@link #onTraversalStart} and/or {@link #onTraversalSuccess} and {@link #afterFirstPass} to perform
 * operations between passes. See the {@link SplitNCigarReads} walker for an example.
 */
public abstract class TwoPassReadWalker extends ReadWalker {

    @Override
    public void traverse() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        final CountingReadFilter countedFilter = makeReadFilter();

        traverseReads(countedFilter, this::firstPassApply);
        logger.info("Finished first pass through the reads");
        afterFirstPass();
        // Need to reinitialize the reads and intervals so they are guaranteed to pass over a file
        initializeReads();
        setReadTraversalBounds();
        logger.info("Starting second pass through the reads");
        traverseReads(countedFilter, this::secondPassApply);
        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Process using the given filter and function.
     * @param countedFilter a filter to apply to all reads.
     * @param f function applied to each read, should produce some useful side effect
     */
    private void traverseReads(final CountingReadFilter countedFilter, final GATKApply f) {
        getTransformedReadStream(countedFilter)
                .forEach(read -> {
                    final SimpleInterval readInterval = getReadInterval(read);
                    f.consume(read,
                            new ReferenceContext(reference, readInterval), // Will create an empty ReferenceContext if reference or readInterval == null
                            new FeatureContext(features, readInterval));   // Will create an empty FeatureContext if features or readInterval == null

                    progressMeter.update(readInterval);
                });
    }

    /**
     * a common abstraction for first and second pass apply functions
     */
    @FunctionalInterface
    private interface GATKApply{
        void consume(GATKRead read, ReferenceContext reference, FeatureContext features);
    }

    /**
     * Process an individual read (with optional contextual information) on the first pass through the reads. Must be
     * implemented by tool authors.
     *
     * Since the {@link TwoPassReadWalker} is inherently stateful, any necessary state should be accumulated by this
     * method during the first pass.
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
     * Process an individual read (with optional contextual information) on the second pass through the reads. Must be
     * implemented by tool authors.
     *
     * The same reads and context will be presented in the same order as were seen during the first pass.
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
    abstract protected void secondPassApply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext);

    /**
     *  Called after all reads in the first pass have been handled by {@link #firstPassApply} and before any reads
     *  are processed by {@link #secondPassApply}.
     *
     *  Tool authors may override in order to update the state of the tool between passes.
     *  The default implementation does nothing.
     */
    protected void afterFirstPass() {}

    /**
     * Not called by {@link TwoPassReadWalker}.  Does nothing.
     *
     * See {@link #firstPassApply} and {@link #secondPassApply} instead.
     */
    @Override
    final public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {}
}
