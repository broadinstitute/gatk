package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;

/**
 * A ReadWalker is a tool that processes a single read at a time from one or multiple sources of reads, with
 * optional contextual information from a reference and/or sets of variants/Features.
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * Reads will be:
 * - Transformed with {@link #makePreReadFilterTransformer()} before filtering.
 * - Filtered with {@link #makeReadFilter()} before post-transformers.
 * - Transformed with {@link #makePostReadFilterTransformer()} before processing.
 * - Passed to {@link #apply(GATKRead, ReferenceContext, FeatureContext)} for processing.
 *
 * ReadWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalSuccess(). See the PrintReadsWithReference walker for an example.
 */
public abstract class ReadWalker extends WalkerBase {

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public String getProgressMeterRecordLabel() { return "reads"; }

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

        setReadTraversalBounds();
    }

    /**
     * Initialize traversal bounds if intervals are specified
     */
    void setReadTraversalBounds() {
        if ( hasUserSuppliedIntervals() ) {
            final SAMSequenceDictionary dict = getHeaderForReads().getSequenceDictionary();
            final boolean traverseUnmapped =
                    intervalArgumentCollection.getTraversalParameters(dict).traverseUnmappedReads();
            reads.setTraversalBounds(new TraversalParameters(userIntervals, traverseUnmapped));
        }
    }

    @Override
    void initializeFeatures() {
        //We override this method to change lookahead of the cache
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                      getGenomicsDBOptions());
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of read-based traversal.
     *
     * The default implementation creates filters using {@link #makeReadFilter} and transformers using
     * {@link #makePreReadFilterTransformer()} {@link #makePostReadFilterTransformer()} and then iterates over all reads, applies
     * the pre-filter transformer, the filter, then the post-filter transformer and hands the resulting reads to the {@link #apply}
     * function of the walker (along with additional contextual information, if present, such as reference bases).
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the engine
     * package.
     */
    @Override
    public void traverse() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        final CountingReadFilter countedFilter = makeReadFilter();
        getTransformedReadStream(countedFilter)
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
    SimpleInterval getReadInterval(final GATKRead read) {
        return !read.isUnmapped() && SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd()) ? new SimpleInterval(read) : null;
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} filter with all default options. Subclasses
     * can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {link #getHeaderForReads}.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(new WellformedReadFilter());
    }

    /**
     * Reset the reads data source so the caller can iterate through the reads again.
     */
    public void resetReadsDataSource() {
        if (!reads.supportsSerialIteration()) {
            // Reinitialize the reads data source to prepare for another iteration.
            reads.close();
            initializeReads();
        }
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
