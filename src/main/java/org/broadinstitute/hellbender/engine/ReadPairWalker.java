package org.broadinstitute.hellbender.engine;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * A ReadPairWalker is a tool that processes a two reads  (with the same read name) at a time from one or multiple
 * sources of reads (queryname sorted).
 *
 * If multiple sources of reads are specified, they are merged together into a single sorted stream of reads.
 *
 * Reads will be:
 * - Transformed with {@link #makePreReadFilterTransformer()} before filtering.
 * - Filtered with {@link #makeReadFilter()} before post-transformers.
 * - Transformed with {@link #makePostReadFilterTransformer()} before processing.
 * - Passed to {@link #apply(Set)} for processing.
 *
 * ReadWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalSuccess(). See the PrintReadsWithReference walker for an example.
 */
public abstract class ReadPairWalker extends GATKTool {

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
        if (getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.queryname) {
            throw new GATKException("This walker requires input in queryname Order. Found:" + getHeaderForReads().getSortOrder());
        }

        setReadTraversalBounds();
    }

    /**
     * Initialize traversal bounds if intervals are specified
     */
    private void setReadTraversalBounds() {
        if ( hasUserSuppliedIntervals() ) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()));
        }
    }

    @Override
    void initializeFeatures() {
        //We override this method to change lookahead of the cache
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                      referenceArguments.getReferencePath());
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Implementation of read-based traversal.
     * Subclasses can override to provide own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation creates filters using {@link #makeReadFilter} and transformers using
     * {@link #makePreReadFilterTransformer()} {@link #makePostReadFilterTransformer()} and then iterates over all reads, applies
     * the pre-filter transformer, the filter, then the post-filter transformer and hands the resulting reads to the {@link #apply}
     * function of the walker (along with additional contextual information, if present, such as reference bases).
     */
    @Override
    public void traverse() {
        // Process each read in the input stream.
        // Supply reference bases spanning each read, if a reference is available.
        final CountingReadFilter countedFilter = makeReadFilter();
        final Set<GATKRead> readSet = new HashSet<>();
        String[] currentReadname = {null};
        getTransformedReadStream(countedFilter)
                .forEachOrdered(read -> {
                    if (!read.getName().equals(currentReadname[0])){
                        apply(readSet);
                        readSet.removeIf((A)->true);
                        currentReadname[0] = read.getName();
                    }
                    readSet.add(read);
                    progressMeter.update(getReadInterval(read));
                });
        if(!readSet.isEmpty()) {
            apply(readSet);
        }
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
        return Arrays.asList(new WellformedReadFilter(),
                new ReadFilterLibrary.NotSecondaryAlignmentReadFilter(),
                new ReadFilterLibrary.NotSupplementaryAlignmentReadFilter());
    }

    /**
     * Process an individual read (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * TODO: Determine whether and to what degree the GATK engine should provide a reduce operation
     * TODO: to complement this operation. At a minimum, we should make apply() return a value to
     * TODO: discourage statefulness in walkers, but how this value should be handled is TBD.
     * @param reads current read set
     */
    public abstract void apply(Set<GATKRead> reads);

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
