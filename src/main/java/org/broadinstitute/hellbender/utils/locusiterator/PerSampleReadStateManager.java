/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.locusiterator;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import htsjdk.samtools.CigarOperator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.downsampling.Downsampler;
import org.broadinstitute.gatk.utils.downsampling.LevelingDownsampler;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * ReadStateManager for a single sample
 *
 * User: depristo
 * Date: 1/13/13
 * Time: 12:28 PM
 */
@Invariant({
        "readStartsAreWellOrdered()",
        "! isDownsampling() || downsamplingTarget > 0",
        "nSites >= 0",
        "nSitesNeedingDownsampling >= 0",
        "nSitesNeedingDownsampling <= nSites"
})
final class PerSampleReadStateManager implements Iterable<AlignmentStateMachine> {
    private final static Logger logger = Logger.getLogger(ReadStateManager.class);
    private final static boolean CAPTURE_DOWNSAMPLING_STATS = false;

    /**
     * A list (potentially empty) of alignment state machines.
     *
     * The state machines must be ordered by the alignment start of their underlying reads, with the
     * lowest alignment starts on the left, and the largest on the right
     */
    private LinkedList<AlignmentStateMachine> readStatesByAlignmentStart = new LinkedList<AlignmentStateMachine>();

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
     * @param LIBSDownsamplingInfo the downsampling params we want to use
     */
    public PerSampleReadStateManager(final LIBSDownsamplingInfo LIBSDownsamplingInfo) {
        this.downsamplingTarget = LIBSDownsamplingInfo.isPerformDownsampling() ? LIBSDownsamplingInfo.getToCoverage() : -1;
        this.levelingDownsampler = LIBSDownsamplingInfo.isPerformDownsampling()
                ? new LevelingDownsampler<LinkedList<AlignmentStateMachine>, AlignmentStateMachine>(LIBSDownsamplingInfo.getToCoverage())
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
    @Ensures("result != null")
    private List<LinkedList<AlignmentStateMachine>> groupByAlignmentStart() {
        final LinkedList<LinkedList<AlignmentStateMachine>> grouped = new LinkedList<LinkedList<AlignmentStateMachine>>();

        AlignmentStateMachine last = null;
        for ( final AlignmentStateMachine stateMachine : readStatesByAlignmentStart ) {
            if ( last == null || stateMachine.getGenomeOffset() != last.getGenomeOffset() ) {
                // we've advanced to a place where the state machine has a different state,
                // so start a new list
                grouped.add(new LinkedList<AlignmentStateMachine>());
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
    @Ensures("result != null")
    private LinkedList<AlignmentStateMachine> flattenByAlignmentStart(final List<LinkedList<AlignmentStateMachine>> grouped) {
        final LinkedList<AlignmentStateMachine> flat = new LinkedList<AlignmentStateMachine>();
        for ( final List<AlignmentStateMachine> l : grouped )
            flat.addAll(l);
        return flat;
    }

    /**
     * Test that the reads are ordered by their alignment starts
     * @return true if well ordered, false otherwise
     */
    private boolean readStartsAreWellOrdered() {
        int lastStart = -1;
        for ( final AlignmentStateMachine machine : readStatesByAlignmentStart ) {
            if ( lastStart > machine.getRead().getAlignmentStart() )
                return false;
            lastStart = machine.getRead().getAlignmentStart();
        }
        return true;
    }

    /**
     * Assumes it can just keep the states linked lists without making a copy
     * @param states the new states to add to this manager
     * @return The change in the number of states, after including states and potentially downsampling.  Note
     * that this return result might be negative, if downsampling is enabled, as we might drop
     * more sites than have been added by the downsampler
     */
    @Requires("states != null")
    public int addStatesAtNextAlignmentStart(final LinkedList<AlignmentStateMachine> states) {
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
        return isEmpty() ? null : readStatesByAlignmentStart.getFirst();
    }

    /**
     * Capture some statistics about the behavior of the downsampling, but only if CAPTURE_DOWNSAMPLING_STATS is true
     */
    @Requires("isDownsampling()")
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

            if ( downsampling || nSites % 10000 == 0 )
                logger.info(String.format("%20s at %s: coverage=%d, max=%d, fraction of downsampled sites=%.2e",
                        message, loc, size(), downsamplingTarget, (1.0 * nSitesNeedingDownsampling / nSites)));
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
    @Ensures("result >= 0")
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
    @Ensures("result != null")
    public Iterator<AlignmentStateMachine> iterator() {
        return readStatesByAlignmentStart.iterator();
    }
}
