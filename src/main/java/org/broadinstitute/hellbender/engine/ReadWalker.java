package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;

import java.util.stream.StreamSupport;

/**
 * A ReadWalker is a tool that processes a single read at a time from one or multiple sources of reads, with
 * optional contextual information from a reference and/or sets of variants/Features.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * ReadWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalDone(). See the PrintReadsWithReference walker for an example.
 */
public abstract class ReadWalker extends GATKTool {

    @Argument(fullName = "disable_all_read_filters", shortName = "f", doc = "Disable all read filters", common = false, optional = true)
    public boolean disable_all_read_filters = false;

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( hasIntervals() ) {
            reads.setIntervalsForTraversal(intervalsForTraversal);
        }
    }

    /**
     * Implementation of read-based traversal.
     * Subclasses can override to provide own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation creates filters using {@link #makeReadFilter}
     * and then iterates over all reads, applies the filter and hands the resulting reads to the {@link #apply}
     * function of the walker (along with additional contextual information, if present, such as reference bases).
     */
    @Override
    public void traverse() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        ReadFilter filter = disable_all_read_filters ? ReadFilterLibrary.ALLOW_ALL_READS : makeReadFilter();

        StreamSupport.stream(reads.spliterator(), false)
                .filter(filter)
                .forEach(read -> {
                    final SimpleInterval readInterval = read.getReadUnmappedFlag() ? null :
                                                                                     new SimpleInterval(read);
                    apply(read,
                          new ReferenceContext(reference, readInterval), // Will create an empty ReferenceContext if reference or readInterval == null
                          new FeatureContext(features, readInterval));   // Will create an empty FeatureContext if features or readInterval == null
                });
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads before calling {@link #apply}.
     * The default implementation uses the {@link org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary.WELLFORMED} filter with all default options.
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can extend to provide own filters (ie override and call super).
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter} composition methods.
     */
    public ReadFilter makeReadFilter(){
          return ReadFilterLibrary.WELLFORMED;
    }

    /**
     * Process an individual read (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * TODO: Determine whether and to what degree the GATK engine should provide a reduce operation
     * TODO: to complement this operation. At a minimum, we should make apply() return a value to
     * TODO: discourage statefulness in walkers, but how this value should be handled is TBD.
     * @param read current read
     * @param referenceContext Reference bases spanning the current read. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current read's interval
     *                         by invoking {@link org.broadinstitute.hellbender.engine.ReferenceContext#setWindow}
     *                         on this object before calling {@link org.broadinstitute.hellbender.engine.ReferenceContext#getBases}
     * @param featureContext Features spanning the current read. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply( SAMRecord read, ReferenceContext referenceContext, FeatureContext featureContext );

    /**
     * Shutdown data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
