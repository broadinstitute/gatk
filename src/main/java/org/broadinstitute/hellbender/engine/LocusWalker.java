package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentContextIteratorBuilder;
import org.broadinstitute.hellbender.utils.locusiterator.LIBSDownsamplingInfo;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * A LocusWalker is a tool that processes reads that overlap a single position in a reference at a time from
 * one or multiple sources of reads, with optional contextual information from a reference and/or sets of
 * variants/Features.
 *
 * Reads will be:
 * - Transformed with {@link #makePreReadFilterTransformer()} before filtering.
 * - Filtered with {@link #makeReadFilter()} before post-transformers.
 * - Transformed with {@link #makePostReadFilterTransformer()} before processing.
 * - Include them into the {@link AlignmentContext} before passing to {@link #apply(AlignmentContext, ReferenceContext, FeatureContext)}.
 *
 * LocusWalker authors must implement the apply() method to process each position, and may optionally implement
 * onTraversalStart(), onTraversalSuccess() and/or closeTool().
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public abstract class LocusWalker extends GATKTool {

    @Argument(fullName = "maxDepthPerSample", shortName = "maxDepthPerSample", doc = "Maximum number of reads to retain per sample per locus. Reads above this threshold will be downsampled. Set to 0 to disable.", optional = true)
    protected int maxDepthPerSample = defaultMaxDepthPerSample();

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

    @Override
    public String getProgressMeterRecordLabel() { return "loci"; }

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
     * Does this tool emit information for uncovered loci? Tools that do should override to return {@code true}.
     *
     * NOTE:  Typically, this should only be used when intervals are specified.
     * NOTE:  If MappedReadFilter is removed, then emitting empty loci will fail.
     * NOTE:  If there is no available sequence dictionary and this is set to true, there should be a failure.  Please
     *  consider requiring reads and/or references for all tools that wish to set this to {@code true}.
     *
     * @return {@code true} if this tool requires uncovered loci information to be emitted, {@code false} otherwise
     */
    public boolean emitEmptyLoci() {
        return false;
    }


    /**
     * Returns default value for the {@link #maxDepthPerSample} parameter, if none is provided on the command line.
     * Default implementation returns 0 (no downsampling by default).
     */
    protected int defaultMaxDepthPerSample() {
        return 0;
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} and {@link ReadFilterLibrary.MappedReadFilter} filter
     * with all default options. Subclasses can override to provide alternative filters.
     *
     * Note: this method is called before command line parsing begins, and thus before a SAMFileHeader is
     * available through {link #getHeaderForReads}.
     *
     * @return List of individual filters to be applied for this tool.
     */
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }

    /** Returns the downsampling info using {@link #maxDepthPerSample} as target coverage. */
    protected final LIBSDownsamplingInfo getDownsamplingInfo() {
        if (maxDepthPerSample < 0) {
            throw new CommandLineException.BadArgumentValue("maxDepthPerSample",
                    String.valueOf(maxDepthPerSample),
                    "should be a positive number");
        }
        return (maxDepthPerSample == 0) ? LocusIteratorByState.NO_DOWNSAMPLING : new LIBSDownsamplingInfo(true, maxDepthPerSample);
    }

    /**
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        if ( hasUserSuppliedIntervals() ) {
            reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()));
        }
    }

    /**
     * Implementation of locus-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation iterates over all positions in the reference covered by reads (filtered and transformed)
     * for all samples in the read groups, using the downsampling method provided by {@link #getDownsamplingInfo()}
     * and including deletions only if {@link #includeDeletions()} returns {@code true}.
     */
    @Override
    public void traverse() {
        final SAMFileHeader header = getHeaderForReads();
        // get the samples from the read groups
        final Set<String> samples = header.getReadGroups().stream()
                                          .map(SAMReadGroupRecord::getSample)
                                          .collect(Collectors.toSet());
        final CountingReadFilter countedFilter = makeReadFilter();
        // get the filter and transformed iterator
        final Iterator<GATKRead> readIterator = getTransformedReadStream(countedFilter).iterator();

        final AlignmentContextIteratorBuilder alignmentContextIteratorBuilder = new AlignmentContextIteratorBuilder();
        alignmentContextIteratorBuilder.setDownsamplingInfo(getDownsamplingInfo());
        alignmentContextIteratorBuilder.setEmitEmptyLoci(emitEmptyLoci());
        alignmentContextIteratorBuilder.setIncludeDeletions(includeDeletions());
        alignmentContextIteratorBuilder.setKeepUniqueReadListInLibs(keepUniqueReadListInLibs());
        alignmentContextIteratorBuilder.setIncludeNs(includeNs());

        final Iterator<AlignmentContext> iterator = alignmentContextIteratorBuilder.build(
                readIterator, header, userIntervals, getBestAvailableSequenceDictionary(),
                hasReference());

        // iterate over each alignment, and apply the function
        iterator.forEachRemaining(alignmentContext -> {
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

    /**
     *  The emit empty loci parameter comes with several pitfalls when used incorrectly.  Here we check and either give
     *   warnings or errors.
     */
    protected void validateEmitEmptyLociParameters() {
        if (emitEmptyLoci()) {
            if (getBestAvailableSequenceDictionary() == null) {
                throw new UserException.MissingReference("Could not create a sequence dictionary nor find a reference.  Therefore, emitting empty loci is impossible and this tool cannot be run.  The easiest fix here is to specify a reference dictionary.");
            }
            if (!hasReference() && !hasUserSuppliedIntervals()) {
                logger.warn("****************************************");
                logger.warn("* Running this tool without a reference nor intervals can yield unexpected results, since it will emit results for loci with no reads.  A sequence dictionary has been found.  The easiest way avoid this message is to specify a reference.");
                logger.warn("****************************************");
            }
        }
    }
}
