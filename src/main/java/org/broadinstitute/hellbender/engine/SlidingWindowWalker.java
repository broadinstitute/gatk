package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import java.util.ArrayList;

/**
 * An SlidingWindowWalker is a tool that processes a single window over the genome (constructed with a window-size and window-step)
 * at a time, with the ability to query optional overlapping sources of reads, reference data, and/or variants/features.
 *
 * SlidingWindow authors must implement the apply() method to process each position, and may optionally implement
 * onTraversalStart() and/or onTraversalDone().
 *
 * @author Daniel Gómez-Sánchez (magicDGS)
 */
public abstract class SlidingWindowWalker extends GATKTool {

	@Argument(fullName = "windowSize", shortName = "windowSize", doc = "Window size for iterate over the variants", common = false, optional = false)
	public int windowSize;

	@Argument(fullName = "windowStep", shortName = "windowStep", doc = "Window step for iterate over the variants", common = false, optional = false)
	public int windowStep;

	/**
	 * This is the array of intervals that configure the windows included in the tool
	 */
	protected final ArrayList<SimpleInterval> windows = new ArrayList<>();

	/**
	 * Initialize data sources for traversal.
	 *
	 * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
	 */
	@Override
	protected final void onStartup() {
		super.onStartup();
		// the tool needs a sequence dictionary for get the windows
		if(getBestAvailableSequenceDictionary() == null) {
			throw new UserException("Tool " + getClass().getSimpleName() + " requires some source for sequence dictionary, but none were provided");
		}
		// construct the windows for the available sequence dictionary
		constructWindows();
	}

	/**
	 * Construct the windows to iterate from.
	 */
	private void constructWindows() {
		logger.info("Computing windows for the sequence");
		for(final SAMSequenceRecord seq: getBestAvailableSequenceDictionary().getSequences()) {
			for(int i = 1; i <= seq.getSequenceLength(); i+=windowStep) {
				final SimpleInterval currentWindow = new SimpleInterval(seq.getSequenceName(), i, i+windowSize);
				if(hasIntervals()) {
					// only adding windows that overlaps with the intervals for traverse
					if(intervalsForTraversal.stream().filter(interval -> interval.overlaps(currentWindow)).findAny().isPresent()) {
						windows.add(currentWindow);
					}
				} else {
					windows.add(currentWindow);
				}
			}
		}
		logger.info("In total ", windows.size(), " windows will be analyzed regarding the intervals provided");
	}

	/**
	 * Implementation of window-based traversing.
	 * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
	 *
	 * The default implementation iterates over every interval stored in {@link #windows} and pass all the available information
	 * to {@link #apply}
	 */
	@Override
	public void traverse() {
		// iterate over the windows
		for ( final SimpleInterval w : windows ) {
			apply(w,
				new ReadsContext(reads, w),
				new ReferenceContext(reference, w),
				new FeatureContext(features, w));
			progressMeter.update(w);
		}
	}

	/**
	 * Process an individual window. Must be implemented by tool authors.
	 * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
	 * as possible.
	 *
	 * @param window current window being processed
	 * @param readsContext Reads overlapping the current window. Will be an empty, but non-null, context object
	 *                     if there is no backing source of reads data (in which case all queries on it will return
	 *                     an empty array/iterator)
	 * @param referenceContext Reference bases spanning the current window. Will be an empty, but non-null, context object
	 *                         if there is no backing source of reference data (in which case all queries on it will return
	 *                         an empty array/iterator). Can request extra bases of context around the current read's interval
	 *                         by invoking {@link ReferenceContext#setWindow}
	 *                         on this object before calling {@link ReferenceContext#getBases}
	 * @param featureContext Features spanning the current window. Will be an empty, but non-null, context object
	 *                       if there is no backing source of Feature data (in which case all queries on it will return an
	 *                       empty List).
	 */
	public abstract void apply(SimpleInterval window, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext);


	/**
	 * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
	 */
	@Override
	protected final void onShutdown() {
		// Overridden only to make final so that concrete tool implementations don't override
		super.onShutdown();
	}
}
