package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadFilterArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A ReadWalker is a tool that processes a single read at a time from one or multiple sources of reads, with
 * optional contextual information from a reference and/or sets of variants/Features.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * ReadWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalSuccess(). See the PrintReadsWithReference walker for an example.
 */
public abstract class ReadWalker extends GATKTool {

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * This number controls the size of the cache for our FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 1_000;

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( hasIntervals() ) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()));
        }
    }

    @Override
    void initializeFeatures() {
        //We override this method to change lookahead of the cache
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
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
        final CountingReadFilter countedFilter = makeReadFilter();

        StreamSupport.stream(reads.spliterator(), false)
                .filter(countedFilter)
                .forEach(read -> {
                    final SimpleInterval readInterval = getReadInterval(read);
                    apply(read,
                          new ReferenceContext(reference, readInterval), // Will create an empty ReferenceContext if reference or readInterval == null
                          new FeatureContext(features, readInterval));   // Will create an empty FeatureContext if features or readInterval == null

                    progressMeter.update(readInterval);
                });

        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Returns an interval for the read.
     * Note: some walkers must be able to work on any read, including those whose coordinates do not form a valid SimpleInterval.
     * So here we check this condition and create null intervals for such reads.
     */
    private SimpleInterval getReadInterval(final GATKRead read) {
        return !read.isUnmapped() && SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd()) ? new SimpleInterval(read) : null;
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads before calling {@link #apply}.
     * The default implementation combines the default read filters for this tool (returned by
     * {@link org.broadinstitute.hellbender.engine.ReadWalker#getDefaultReadFilter} with any read filter command
     * line arguments specified by the user; wraps each filter in the resulting list with a CountingReadFilter;
     * and returns a single composite filter resulting from the list by and'ing them together.
     *
     * Default tool implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter
     * objects is strongly discouraged.
     *
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter}
     * composition methods.
     */
    public CountingReadFilter makeReadFilter(){
        // Merge the default filters with the users's command line read filter requests, then instantiate
        // and wrap the resulting individual filters with CountingReadFilter.
        final ReadFilterArgumentCollection rfc = readArguments.getReadFilterArgumentCollection();
        return rfc.getMergedCommandLineFilters(getDefaultReadFilters())
                .stream()
                .map(f -> new CountingReadFilter(
                        f.name(),
                        f.getFilterInstance(getHeaderForReads(), rfc)))
                .reduce(new CountingReadFilter(
                        "ALLOW_ALL",
                        ReadFilterArgumentCollection.CommandLineReadFilter.ALLOW_ALL.getFilterInstance(getHeaderForReads(), rfc)),
                        (f1, f2) -> f1.and(f2));
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} filter with all default options. Subclasses
     * can override to provide alternative filters.
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilterArgumentCollection.CommandLineReadFilter> getDefaultReadFilters() {
        return Arrays.asList(ReadFilterArgumentCollection.CommandLineReadFilter.WELLFORMED);
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
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current read. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext );

    /**
     * Shutdown data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalSuccess() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
