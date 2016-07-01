package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReadFilterArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.iterators.IntervalOverlappingIterator;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A LocusWalker is a tool that processes reads that overlap a single position in a reference at a time from
 * one or multiple sources of reads, with optional contextual information from a reference and/or sets of
 * variants/Features.
 *
 * LocusWalker authors must implement the apply() method to process each position, and may optionally implement
 * onTraversalStart(), onTraversalSuccess() and/or closeTool().
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class LocusWalker extends GATKTool {

    /**
     * Should the LIBS keep unique reads? Tools that do should override to return {@code true}.
     */
    protected boolean keepUniqueReadListInLibs() {
        return false;
    }

    /**
     * LocusWalkers requires read sources
     */
    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * Does this tool require deletions in the AlignmentContext? Tools that don't should override to return {@code false}.
     *
     * @return {@code true} if this tool requires deletions, {@code false} otherwise
     */
    public boolean includeDeletions() {
        return true;
    }

    /**
     * Does this tool require Ns in the AlignmentContext? Tools that do should override to return {@code true}.
     *
     * @return {@code true} if this tool requires Ns, {@code false} otherwise
     */
    public boolean includeNs() {
        return false;
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
    public CountingReadFilter makeReadFilter() {
        final ReadFilterArgumentCollection rfc = readArguments.getReadFilterArgumentCollection();
        // Merge the default filters with the users's command line read filter requests, then wrap the
        // resulting individual filters with CountingReadFilter.
        return rfc.getMergedCommandLineFilters(getDefaultReadFilters()).stream()
                .map(f -> new CountingReadFilter(
                        f.name(),
                        f.getFilterInstance(getHeaderForReads(), readArguments.getReadFilterArgumentCollection())))
                .reduce(new CountingReadFilter(
                                "ALLOW_ALL",
                                ReadFilterArgumentCollection.CommandLineReadFilter.ALLOW_ALL.getFilterInstance(getHeaderForReads(), null)),
                        (f1, f2) -> f1.and(f2));
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} and {@link ReadFilterLibrary.MAPPED}filter
     * with all default options. Subclasses can override to provide alternative filters.
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilterArgumentCollection.CommandLineReadFilter> getDefaultReadFilters() {
        final List<ReadFilterArgumentCollection.CommandLineReadFilter> defaultFilters = new ArrayList<>();
        defaultFilters.add(ReadFilterArgumentCollection.CommandLineReadFilter.WELLFORMED);
        return defaultFilters;
    }

    /**
     * Get the information about how to downsample the reads. By default {@link DownsamplingMethod#NONE} is returned.
     * Subclasses should override it to provide their own downsampling methods.
     *
     * @return the downsampling method for the reads
     */
    public DownsamplingMethod getDownsamplingMethod() {
        // TODO: change when downsampling command line options are implemented
        return DownsamplingMethod.NONE;
    }

    /**
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        if ( hasIntervals() ) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()));
        }
    }

    /**
     * Implementation of locus-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation iterates over all positions in the reference covered by reads for all samples in the read groups, using
     * the downsampling method provided by {@link #getDownsamplingMethod()}
     * and including deletions only if {@link #includeDeletions()} returns {@code true}.
     */
    @Override
    public void traverse() {
        final SAMFileHeader header = getHeaderForReads();
        // get the samples from the read groups
        final Set<String> samples = header.getReadGroups().stream()
                                          .map(SAMReadGroupRecord::getSample)
                                          .collect(Collectors.toSet());
        CountingReadFilter countedFilter = makeReadFilter();
        // get the LIBS
        LocusIteratorByState libs = new LocusIteratorByState(new ReadFilteringIterator(reads.iterator(), countedFilter), getDownsamplingMethod(), includeDeletions(), includeNs(), keepUniqueReadListInLibs(), samples, header);
        // prepare the iterator
        Spliterator<AlignmentContext> iterator = (hasIntervals()) ? new IntervalOverlappingIterator<>(libs, intervalsForTraversal, header.getSequenceDictionary()).spliterator() : libs.spliterator();
        // iterate over each alignment, and apply the function
        StreamSupport.stream(iterator, false)
            .forEach(alignmentContext -> {
                        final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext);
                        apply(alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(features, alignmentInterval));
                        progressMeter.update(alignmentInterval);
                }
            );
        logger.info(countedFilter.getSummaryLine());
    }

    /**
     * Process an individual AlignmentContext (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param alignmentContext current alignment context
     * @param referenceContext Reference bases spanning the current locus. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current locus
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current locus. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext);

    /**
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalSuccess() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
