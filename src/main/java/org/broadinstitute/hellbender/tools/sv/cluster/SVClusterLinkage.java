package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public abstract class SVClusterLinkage<T extends SVLocatable> {

    /**
     * Returns whether two given items cluster.
     * @param a first item
     * @param b second item
     */
    abstract boolean areClusterable(final T a, final T b);

    /**
     * Returns the maximum feasible starting position of any other item with the given item. That is, given item A and
     * X = getMaxClusterableStartingPosition(A), then for any item B on the current contig,
     * Y = start(B) > X => clusterTogether(A, B) == false. Note that this is an upper-bound, but tighter estimates
     * can greatly improve performance.
     * @param item item in question
     * @return max feasible clusterable start coordinate on the current contig
     */
    abstract int getMaxClusterableStartingPosition(final T item);

    /**
     * Compute max feasible starting position of any other item for all items in the given collection. Note the items
     * must all have the same starting contig.
     */
    public int getMaxClusterableStartingPosition(final Collection<T> items) {
        final List<String> contigA = items.stream().map(T::getContigA).distinct().collect(Collectors.toList());
        if (contigA.size() > 1) {
            throw new IllegalArgumentException("Items start on multiple contigs");
        }
        return items.stream().mapToInt(item -> getMaxClusterableStartingPosition(item)).max().getAsInt();
    }

    /**
     * Checks for minimum fractional sample overlap of the two sets. Defaults to true if both sets are empty.
     */
    protected static boolean hasSampleSetOverlap(final Set<String> samplesA, final Set<String> samplesB, final double minSampleOverlap) {
        final int denom = Math.max(samplesA.size(), samplesB.size());
        if (denom == 0) {
            return true;
        }
        final double sampleOverlap = getSampleSetOverlap(samplesA, samplesB) / (double) denom;
        return sampleOverlap >= minSampleOverlap;
    }

    /**
     * Returns number of overlapping items
     */
    protected static int getSampleSetOverlap(final Collection<String> a, final Set<String> b) {
        return (int) a.stream().filter(b::contains).count();
    }

    /**
     * Returns true if there is sufficient fractional carrier sample overlap in the two records.
     */
    protected static boolean hasSampleOverlap(final SVCallRecord a, final SVCallRecord b, final double minSampleOverlap) {
        if (minSampleOverlap > 0) {
            final Set<String> samplesA = a.getCarrierSamples();
            final Set<String> samplesB = b.getCarrierSamples();
            return hasSampleSetOverlap(samplesA, samplesB, minSampleOverlap);
        } else {
            return true;
        }
    }

}
