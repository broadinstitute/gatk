package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.collections4.SetUtils;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * An implementation of {@link LocusWalker} that supports arbitrary interval side inputs.
 *
 * The class works as follows: Before traversal {@link #getIntervalObjectsToQueryOver()} will be called to generate a global
 * list of Locatable that will be stored in memory for the duration of the traversal. For each site {@link #apply(AlignmentContext, ReferenceContext, FeatureContext, Set)}} will be called
 * with a set consisting of all the provided overlapping intervals. The first time any Locus overlaps with the traversal
 * {@link #onIntervalStart(Locatable)} will be called once. Once a the traversal locus no longer overlaps a Locatable, {@link #onIntervalEnd(Locatable)}
 * will be called once to perform any necessary post-processing. Otherwise See {@link LocusWalker} for general implementation details.
 *
 * NOTE: If there are Locatables provided by {@link #getIntervalObjectsToQueryOver()} that are never covered by the traversal of
 * the tool, {@link #onIntervalStart(Locatable)} and {@link #onIntervalEnd(Locatable)} will not be called on those intervals.
 */
public abstract class LocusWalkerByInterval extends LocusWalker {

    private OverlapDetector<Locatable> intervalsToTrack = null;
    private Set<Locatable> previousIntervals = new LinkedHashSet<>();

    /**
     * Implementation of locus-based traversal.
     *
     * This implementation behaves similarly to {@link LocusWalker#traverse()} in that it iterates over all positions in the reference
     * covered by filtered and transformed reads including deletions only if {@link #includeDeletions()} returns {@code true}.
     *
     * This method also keeps track of interval objects provided by {@link #getIntervalObjectsToQueryOver()} and constructs a
     * global overlaps detector for all of the intervals which is used by the {@link #apply(AlignmentContext, ReferenceContext, FeatureContext)}
     * method to track which locatable still have active hooks. This method also makes sure to close out the list of previous intervals
     * when traversal has completed so that writers can be populated.
     */
    @Override
    public void traverse() {
        final CountingReadFilter countedFilter = makeReadFilter();
        final Iterator<AlignmentContext> iterator = getAlignmentContextIterator(countedFilter);

        intervalsToTrack = OverlapDetector.create(getIntervalObjectsToQueryOver());

        // iterate over each alignment, and apply the function
        iterator.forEachRemaining(alignmentContext -> {
                    final SimpleInterval alignmentInterval = new SimpleInterval(alignmentContext);
                    apply(alignmentContext, new ReferenceContext(reference, alignmentInterval), new FeatureContext(features, alignmentInterval));
                    progressMeter.update(alignmentInterval);
                }
        );
        for (Locatable l : previousIntervals) {
            onIntervalEnd(l);
        }
        logger.info(countedFilter.getSummaryLine());
    }

    @Override
    // A locusWalkerByInterval requires intervals be specified
    public final boolean requiresIntervals() {
        return true;
    }

    /**
     * Tool-specified list of Locatable objects (which have been read into memory) that will have overlaps queried at each locus
     *
     * @return A list of Locatable Objects that will be checked for overlap with each genomic location we traverse
     */
    public abstract List<Locatable> getIntervalObjectsToQueryOver();

    @Override
    public final void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        Set<Locatable> currentIntervals = intervalsToTrack.getOverlaps(alignmentContext);
        Set<Locatable> passedIntervals = SetUtils.difference(previousIntervals, currentIntervals);
        Set<Locatable> newIntervals = SetUtils.difference(currentIntervals, previousIntervals);
        previousIntervals = currentIntervals;

        // Close out existing intervals that have been passed
        for(Locatable l : passedIntervals) {
            onIntervalEnd(l);
        }
        // Open new intervals that have been reached
        for(Locatable l : newIntervals) {
            onIntervalStart(l);
        }

        apply(alignmentContext, referenceContext, featureContext, currentIntervals);
    }

    /**
     * Process an individual AlignmentContext (with optional contextual information). Must be implemented by tool authors.
     * Will provide a set of all intervals overlapping the given locus.
     *
     * @param alignmentContext current alignment context
     * @param referenceContext Reference bases spanning the current locus. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current locus
     *                         by invoking {@link ReferenceContext#setWindow} on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current locus. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     * @param activeIntervals Locatables from the set provided by getIntervalObjectsToQueryOver() spanning the current locus.
     */
    public abstract void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext, Set<Locatable> activeIntervals);

    /**
     * Perform any initialization needed the first time a provided interval is seen.
     *
     * @param activeInterval Locatable an interval from the set provided by getIntervalObjectsToQueryOver() to the walker to be initialized
     */
    public abstract void onIntervalStart(Locatable activeInterval);

    /**
     * Perform any closing operations needed once the provided getIntervalObjectsToQueryOver() interval has been passed by the tool
     *
     * NOTE: This will only ever be called on an interval that has been covered by {@link #onIntervalStart(Locatable)}
     *
     * @param activeInterval Locatable provided to the walker to be closed
     */
    public abstract void onIntervalEnd(Locatable activeInterval);
}
