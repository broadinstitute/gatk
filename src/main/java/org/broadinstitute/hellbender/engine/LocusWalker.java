package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
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
 * onTraversalStart() and/or onTraversalDone().
 *
 * @author Daniel Gómez-Sánchez (magicDGS)
 */
public abstract class LocusWalker extends GATKTool {

    @Argument(fullName = "disable_all_read_filters", shortName = "f", doc = "Disable all read filters", common = false, optional = true)
    public boolean disableAllReadFilters = false;

    /**
     * Should the LIBS keep unique reads?
     */
    protected boolean KEEP_UNIQUE_READ_LIST_IN_LIBS = false;

    /**
     * LocusWalkers requires read sources
     */
    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * Does this tool require deletions in the AlignmentContext? Tools that do should override to return true.
     *
     * @return {@code true} if this tool requires deletions, {@code false} otherwise
     */
    public boolean includeDeletions() {
        return true;
    }

    /**
     * Returns the read filter (simple or composite) that will be applied to the reads in each window.
     *
     * The default implementation uses the {@link org.broadinstitute.hellbender.engine.filters.WellformedReadFilter} filter with all default options,
     * as well as the {@link ReadFilterLibrary#MAPPED} filter.
     *
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can override to provide their own filters
     * Multiple filters can be composed by using {@link org.broadinstitute.hellbender.engine.filters.ReadFilter} composition methods.
     */
    public CountingReadFilter makeReadFilter(){
        return new CountingReadFilter("Wellformed", new WellformedReadFilter(getHeaderForReads()))
                .and(new CountingReadFilter("Mapped", ReadFilterLibrary.MAPPED));
    }

    /**
     * Get the information about how to downsample the reads. By default {@link org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod#NONE} is returned.
     * Subclasses should override it to provide own downsampling methods.
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
            reads.setIntervalsForTraversal(intervalsForTraversal);
        }
    }

    /**
     * Implementation of locus-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     *
     * The default implementation iterates over all positions in the reference for all samples in the read groups, using
     * the downsampling method provided by {@link #getDownsamplingMethod()}
     * and including deletions if {@link #includeDeletions()} returns {@code true}.
     *
     * {@link #intervalsForTraversal} are not used in this implementation
     */
    @Override
    public void traverse() {
        final SAMFileHeader header = getHeaderForReads();
        // get the samples from the read groups
        final Set<String> samples = header.getReadGroups().stream()
                                          .map(SAMReadGroupRecord::getSample)
                                          .collect(Collectors.toSet());
        CountingReadFilter countedFilter = disableAllReadFilters ?
                new CountingReadFilter("Allow all", ReadFilterLibrary.ALLOW_ALL_READS ) :
                makeReadFilter();
        // TODO: wrap reads.iterator() with ReadFilteringIterator and countedFilter
        // get the LIBS
        LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(), getDownsamplingMethod(), includeDeletions(), KEEP_UNIQUE_READ_LIST_IN_LIBS, samples, header);
		// prepare the iterator
		Spliterator<AlignmentContext> iterator = (hasIntervals()) ? new IntervalOverlappingIterator(libs, intervalsForTraversal).spliterator() : libs.spliterator();
        // iterate over each alignment, and apply the function
        StreamSupport.stream(iterator, false)
            .forEach(alignmentContext -> {
                        final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext.getLocation());
                        apply(alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(features, alignmentInterval));
                        progressMeter.update(alignmentInterval);
                }
            );
    }

	/**
	 * Wrapper for {@link org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState} and intervals to iterate
	 */
	protected class IntervalOverlappingIterator implements Iterable<AlignmentContext>, Iterator<AlignmentContext> {

		private final Iterator<AlignmentContext> iterator;

		private final Iterator<SimpleInterval> intervals;

		private SimpleInterval currentInterval;

		private AlignmentContext next;

		/**
		 * Constructor
		 *
		 * @param iterator		underlying iterator
		 * @param intervals sorted list of intervals to traverse
		 */
		public IntervalOverlappingIterator(Iterator<AlignmentContext> iterator, List<SimpleInterval> intervals) {
			this.iterator = iterator;
			this.intervals = intervals.iterator();
			currentInterval = this.intervals.next();
			advance();
		}

		@Override
		public Iterator<AlignmentContext> iterator() {
			return this;
		}

		@Override
		public boolean hasNext() {
			return next != null;
		}

		@Override
		public AlignmentContext next() {
			if(!hasNext()) {
				throw new NoSuchElementException();
			}
			AlignmentContext toReturn = next;
			advance();
			return toReturn;
		}

		/**
		 * Advance to the next AlignmentContext
		 */
		private void advance() {
			// all sources are finished
			if(!iterator.hasNext() || currentInterval == null) {
				next = null;
			} else {
				next = iterator.next();
				// if the next AlignmentContext is not in the current interval
				if(!currentInterval.overlaps(next.getLocation())) {
					// advance the interval and try with the next one
					currentInterval = (intervals.hasNext()) ? intervals.next() : null;
					advance();
				}
			}
		}

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
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
