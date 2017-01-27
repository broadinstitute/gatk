package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.IntervalOverlappingIterator;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.locusiterator.LIBSDownsamplingInfo;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A MultiPassLocusWalker is a tool that processes reads that overlap a single position in a reference at a time from
 * one or multiple sources of reads, with optional contextual information from a reference and/or sets of
 * variants/Features. Each source of reads is processed one at a time, one after another in serial fashion, unlike
 * LocusWalker.
 *
 * MultiPassLocusWalker authors must implement the apply() method to process each position, and may optionally implement
 * onTraversalStart(), onTraversalSuccess() and/or closeTool().
 */
public abstract class MultiPassLocusWalker extends GATKTool {

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
     * Returns default value for the {@link #maxDepthPerSample} parameter, if none is provided on the command line.
     * Default implementation returns 0 (no downsampling by default).
     */
    protected int defaultMaxDepthPerSample() {
        return 0;
    }

    /**
     * Returns the default list of CommandLineReadFilters that are used for this tool. The filters returned
     * by this method are subject to selective enabling/disabling by the user via the command line. The
     * default implementation uses the {@link WellformedReadFilter} and {@link ReadFilterLibrary.MappedReadFilter}filter
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

    @Override
    void initializeReads() {
        reads = readArguments.getReadPaths().isEmpty() ? null : initializeReadSource(0);
    }

    private ReadsDataSource initializeReadSource( final int readInputIndex ) {
        final List<Path> readPaths = readArguments.getReadPaths();
        Utils.validIndex(readInputIndex, readPaths.size());

        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(readArguments.getReadValidationStringency());
        if (hasReference()) { // pass in reference if available, because CRAM files need it
            factory = factory.referenceSequence(referenceArguments.getReferenceFile());
        }
        else if (hasCramInput()) {
            throw new UserException.MissingReference("A reference file is required when using CRAM files.");
        }

        final List<Path> readIndexPath = readArguments.getReadIndexPaths() == null ? null : Arrays.asList(readArguments.getReadIndexPaths().get(readInputIndex));
        return new ReadsDataSource(Arrays.asList(readPaths.get(readInputIndex)), readIndexPath,
                factory, cloudPrefetchBuffer, (cloudIndexPrefetchBuffer < 0 ? cloudPrefetchBuffer : cloudIndexPrefetchBuffer));
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
     * the downsampling method provided by {@link #getDownsamplingInfo()}
     * and including deletions only if {@link #includeDeletions()} returns {@code true}.
     */
    @Override
    public void traverse() {
        final List<Path> readPaths = readArguments.getReadPaths();

        for ( int readInputIndex = 0; readInputIndex < readPaths.size(); ++readInputIndex ) {
            // Use the Path as the input label for now.
            // TODO: use the logical label assigned to each input on the command line, such as "tumor" or "normal",
            // TODO: once that feature exists
            final String readInputLabel = readPaths.get(readInputIndex).toString();

            final CountingReadFilter countedFilter = makeReadFilter();
            final SAMFileHeader header = getHeaderForReads();
            // get the samples from the read groups
            final Set<String> samples = header.getReadGroups().stream()
                    .map(SAMReadGroupRecord::getSample)
                    .collect(Collectors.toSet());

            // get the LIBS
            final LocusIteratorByState libs = new LocusIteratorByState(new ReadFilteringIterator(reads.iterator(), countedFilter), getDownsamplingInfo(), keepUniqueReadListInLibs(), samples, header, includeDeletions(), includeNs());
            // prepare the iterator
            Spliterator<AlignmentContext> iterator = (hasIntervals()) ? new IntervalOverlappingIterator<>(libs, intervalsForTraversal, header.getSequenceDictionary()).spliterator() : libs.spliterator();
            // iterate over each alignment, and apply the function
            StreamSupport.stream(iterator, false)
                    .forEach(alignmentContext -> {
                                final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext);
                                apply(readInputLabel, alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(features, alignmentInterval));
                                progressMeter.update(alignmentInterval);
                            }
                    );

            logger.info(countedFilter.getSummaryLine());
            advanceToNextReadSource(readInputIndex);
        }
    }

    private void advanceToNextReadSource( final int currentReadInputIndex ) {
        reads.close();
        if ( currentReadInputIndex + 1 < readArguments.getReadPaths().size() ) {
            reads = initializeReadSource(currentReadInputIndex + 1);
            if ( hasIntervals() ) {
                reads.setTraversalBounds(intervalArgumentCollection.getTraversalParameters(getHeaderForReads().getSequenceDictionary()));
            }
        }
    }

    /**
     * Process an individual AlignmentContext (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param inputLabel Label for the current read input
     * @param alignmentContext current alignment context
     * @param referenceContext Reference bases spanning the current locus. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current locus
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current locus. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply(final String inputLabel, AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext);

    /**
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalSuccess() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
