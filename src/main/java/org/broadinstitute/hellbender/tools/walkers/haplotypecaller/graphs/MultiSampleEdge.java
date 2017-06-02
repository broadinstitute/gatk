package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.PriorityQueue;

/**
 * Edge class for connecting nodes in the graph that tracks some per-sample information.
 * <p>
 * This class extends BaseEdge with the additional functionality of tracking the maximum
 * multiplicity seen within any single sample.  The workflow for using this class is:
 * </p>
 * <pre>
 * {@code
 *      MultiSampleEdge e = new MultiSampleEdge(ref, 1)
 *      e.incMultiplicity(1)              // total is 2, per sample is 2, max per sample is 1
 *      e.getPruningMultiplicity()        // = 1
 *      e.flushSingleSampleMultiplicity() // total is 2, per sample is 0, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.incMultiplicity(3)              // total is 5, per sample is 3, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.flushSingleSampleMultiplicity() // total is 5, per sample is 0, max per sample is 3
 *      e.getPruningMultiplicity()        // = 3
 * }
 * </pre>
 */
public final class MultiSampleEdge extends BaseEdge {
    private int currentSingleSampleMultiplicity;
    private final int singleSampleCapacity;
    private final PriorityQueue<Integer> singleSampleMultiplicities;

    /**
     * Create a new MultiSampleEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     * @param multiplicity the number of observations of this edge in this sample
     * @param singleSampleCapacity the max number of samples to track edge multiplicities
     */
    public MultiSampleEdge(final boolean isRef, final int multiplicity, final int singleSampleCapacity) {
        super(isRef, multiplicity);

        Utils.validateArg( singleSampleCapacity > 0, () -> "singleSampleCapacity must be > 0 but found: " + singleSampleCapacity);
        singleSampleMultiplicities = new PriorityQueue<>(singleSampleCapacity);
        singleSampleMultiplicities.add(multiplicity);
        currentSingleSampleMultiplicity = multiplicity;
        this.singleSampleCapacity = singleSampleCapacity;
    }

    @Override
    public MultiSampleEdge copy() {
        return new MultiSampleEdge(isRef(), getMultiplicity(), singleSampleCapacity); // TODO -- should I copy values for other features?
    }

    /**
     * update the single sample multiplicities by adding the current single sample multiplicity to the priority queue, and
     * reset the current single sample multiplicity to 0.
     */
    public void flushSingleSampleMultiplicity() {
        singleSampleMultiplicities.add(currentSingleSampleMultiplicity);
        if( singleSampleMultiplicities.size() == singleSampleCapacity + 1 ) {
            singleSampleMultiplicities.poll(); // remove the lowest multiplicity from the list
        } else if( singleSampleMultiplicities.size() > singleSampleCapacity + 1 ) {
            throw new IllegalStateException("Somehow the per sample multiplicity list has grown too big: " + singleSampleMultiplicities);
        }
        currentSingleSampleMultiplicity = 0;
    }

    @Override
    public void incMultiplicity(final int incr) {
        super.incMultiplicity(incr);
        currentSingleSampleMultiplicity += incr;
    }

    @Override
    public int getPruningMultiplicity() {
        return singleSampleMultiplicities.peek();
    }

    @Override
    public String getDotLabel() {
        return super.getDotLabel() + '/' + getPruningMultiplicity();
    }

    @VisibleForTesting
    int getCurrentSingleSampleMultiplicity() {
        return currentSingleSampleMultiplicity;
    }
}
