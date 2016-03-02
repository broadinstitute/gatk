package org.broadinstitute.hellbender.utils.locusiterator;

import htsjdk.samtools.CigarOperator;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.Downsampler;
import org.broadinstitute.hellbender.utils.downsampling.LevelingDownsampler;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * ReadStateManager for a single sample
 */
final class PerSampleReadStateManager implements Iterable<AlignmentStateMachine> {
    private static final Logger logger = LogManager.getLogger(PerSampleReadStateManager.class);
    private static final boolean CAPTURE_DOWNSAMPLING_STATS = false;

    /**
     * A list (potentially empty) of alignment state machines.
     *
     * The state machines must be ordered by the alignment start of their underlying reads, with the
     * lowest alignment starts on the left, and the largest on the right
     */
    private List<AlignmentStateMachine> readStatesByAlignmentStart = new LinkedList<>();

    private final Downsampler<LinkedList<AlignmentStateMachine>> levelingDownsampler;
    private final int downsamplingTarget;

    /**
     * The number of sites where downsampling has been invoked
     */
    private int nSitesNeedingDownsampling = 0;

    /**
     * The number of sites we've visited
     */
    private int nSites = 0;

    /**
     * Create a new PerSampleReadStateManager with downsampling parameters as requested by LIBSDownsamplingInfo
     * @param info the downsampling params we want to use
     */
    public PerSampleReadStateManager(final LIBSDownsamplingInfo info) {
        Utils.nonNull(info);
        this.downsamplingTarget = info.isPerformDownsampling() ? info.getToCoverage() : -1;
        this.levelingDownsampler = info.isPerformDownsampling()
                ? new LevelingDownsampler<>(info.getToCoverage())
                : null;
    }

    /**
     * Group the underlying readStatesByAlignmentStart into a list of list of alignment state machines,
     * where each list contains machines with a unique genome site.  The outer list is ordered
     * by alignment start.
     *
     * For example, if the flat list has alignment starts [10, 10, 11, 12, 12, 13] then
     * the resulting grouping will be [[10, 10], [11], [12, 12], [13]].
     *
     * @return a non-null list of lists
     */
    private List<LinkedList<AlignmentStateMachine>> groupByAlignmentStart() {
        final LinkedList<LinkedList<AlignmentStateMachine>> grouped = new LinkedList<>();

        AlignmentStateMachine last = null;
        for ( final AlignmentStateMachine stateMachine : readStatesByAlignmentStart ) {
            if ( last == null || stateMachine.getGenomeOffset() != last.getGenomeOffset() ) {
                // we've advanced to a place where the state machine has a different state,
                // so start a new list
                grouped.add(new LinkedList<>());
                last = stateMachine;
            }
            grouped.getLast().add(stateMachine);
        }

        return grouped;
    }

    /**
     * Flattens the grouped list of list of alignment state machines into a single list in order
     * @return a non-null list contains the state machines
     */
    private LinkedList<AlignmentStateMachine> flattenByAlignmentStart(final List<LinkedList<AlignmentStateMachine>> grouped) {
        final LinkedList<AlignmentStateMachine> flat = new LinkedList<>();
        for ( final List<AlignmentStateMachine> l : grouped ) {
            flat.addAll(l);
        }
        return flat;
    }

    /**
     * Assumes it can just keep the states linked lists without making a copy
     * @param states the new states to add to this manager
     * @return The change in the number of states, after including states and potentially downsampling.  Note
     * that this return result might be negative, if downsampling is enabled, as we might drop
     * more sites than have been added by the downsampler
     */
    public int addStatesAtNextAlignmentStart(final List<AlignmentStateMachine> states) {
        Utils.nonNull(states);
        if ( states.isEmpty() ) {
            return 0;
        }

        readStatesByAlignmentStart.addAll(states);
        int nStatesAdded = states.size();

        if ( isDownsampling() && readStatesByAlignmentStart.size() > downsamplingTarget ) {
            // only go into the downsampling branch if we are downsampling and the coverage > the target
            captureDownsamplingStats();
            levelingDownsampler.submit(groupByAlignmentStart());
            levelingDownsampler.signalEndOfInput();

            nStatesAdded -= levelingDownsampler.getNumberOfDiscardedItems();

            // use returned List directly rather than make a copy, for efficiency's sake
            readStatesByAlignmentStart = flattenByAlignmentStart(levelingDownsampler.consumeFinalizedItems());
            levelingDownsampler.resetStats();
        }

        return nStatesAdded;
    }

    /**
     * Is downsampling enabled for this manager?
     * @return true if we are downsampling, false otherwise
     */
    private boolean isDownsampling() {
        return levelingDownsampler != null;
    }

    /**
     * Get the leftmost alignment state machine, or null if the read states is empty
     * @return a potentially null AlignmentStateMachine
     */
    public AlignmentStateMachine getFirst() {
        return isEmpty() ? null : readStatesByAlignmentStart.get(0);
    }

    /**
     * Capture some statistics about the behavior of the downsampling, but only if CAPTURE_DOWNSAMPLING_STATS is true
     */
    private void captureDownsamplingStats() {
        if ( CAPTURE_DOWNSAMPLING_STATS ) {
            nSites++;
            final int loc = getFirst().getGenomePosition();
            String message = "Pass through";
            final boolean downsampling = size() > downsamplingTarget;
            if ( downsampling ) {
                nSitesNeedingDownsampling++;
                message = "Downsampling";
            }

            if ( downsampling || nSites % 10000 == 0 ) {
                logger.info(String.format("%20s at %s: coverage=%d, max=%d, fraction of downsampled sites=%.2e",
                        message, loc, size(), downsamplingTarget, (1.0 * nSitesNeedingDownsampling / nSites)));
            }
        }
    }

    /**
     * Is there at least one alignment for this sample in this manager?
     * @return true if there's at least one alignment, false otherwise
     */
    public boolean isEmpty() {
        return readStatesByAlignmentStart.isEmpty();
    }

    /**
     * Get the number of read states currently in this manager
     * @return the number of read states
     */
    public int size() {
        return readStatesByAlignmentStart.size();
    }

    /**
     * Advances all read states forward by one element, removing states that are
     * no long aligned to the current position.
     * @return the number of states we're removed after advancing
     */
    public int updateReadStates() {
        int nRemoved = 0;
        final Iterator<AlignmentStateMachine> it = iterator();
        while (it.hasNext()) {
            final AlignmentStateMachine state = it.next();
            final CigarOperator op = state.stepForwardOnGenome();
            if (op == null) {
                // we discard the read only when we are past its end AND indel at the end of the read (if any) was
                // already processed. Keeping the read state that returned null upon stepForwardOnGenome() is safe
                // as the next call to stepForwardOnGenome() will return null again AND will clear hadIndel() flag.
                it.remove();                                                // we've stepped off the end of the object
                nRemoved++;
            }
        }

        return nRemoved;
    }

    /**
     * Iterate over the AlignmentStateMachine in this manager in alignment start order.
     * @return a valid iterator
     */
    @Override
    public Iterator<AlignmentStateMachine> iterator() {
        return readStatesByAlignmentStart.iterator();
    }
}
