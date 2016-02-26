package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;

import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * A LocusWalker is a tool that processes reads that overlaps a single position in a reference at a time from
 * one or multiple sources of reads, with optional contextual information from a reference and/or sets of
 * variants/Features.
 *
 * LocusWalker authors must implement the apply() method to process each position, and may optionally implement
 * onTraversalStart() and/or onTraversalDone().
 *
 * @author Daniel Gómez-Sánchez (magicDGS)
 */
public abstract class LocusWalker extends GATKTool {

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
		// Overridden only to make final so that concrete tool implementations don't override
		super.onStartup();
	}

	/**
	 * Implementation of AlignmentContext traversal.
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
		final SAMFileHeader header = reads.getHeader();
		// get the samples from the read groups
		final Set<String> samples = header.getReadGroups().stream()
										  .map(SAMReadGroupRecord::getSample)
										  .collect(Collectors.toSet());
		// get the LIBS
		LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(), getDownsamplingMethod(), includeDeletions(), false, samples, header);
		// iterate over each alignment, and apply the function
		StreamSupport.stream(libs.spliterator(), false)
			.forEach(alignmentContext -> {
				final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext.getLocation());
				apply(alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(features, alignmentInterval));
				progressMeter.update(alignmentInterval);
				}
			);
	}

	/**
	 * Process an individual AlignmentContext (with optional contextual information). Must be implemented by tool authors.
	 * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
	 * as possible.
	 *
	 * @param alignmentContext current alignment context
	 * @param referenceContext Reference bases spanning the current read. Will be an empty, but non-null, context object
	 *                         if there is no backing source of reference data (in which case all queries on it will return
	 *                         an empty array/iterator). Can request extra bases of context around the current read's interval
	 *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
	 * @param featureContext Features spanning the current read. Will be an empty, but non-null, context object
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
