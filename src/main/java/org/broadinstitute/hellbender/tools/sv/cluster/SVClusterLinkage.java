package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

public interface SVClusterLinkage<T extends SVLocatable> {

    /**
     * Returns whether two given items cluster.
     * @param a first item
     * @param b second item
     */
    boolean areClusterable(final T a, final T b);

    /**
     * Returns the maximum feasible starting position of any other item with the given item. That is, given item A and
     * X = getMaxClusterableStartingPosition(A), then for any item B on the current contig,
     * Y = start(B) > X => clusterTogether(A, B) == false. Note that this is an upper-bound, but tighter estimates
     * can greatly improve performance.
     * @param item item in question
     * @return max feasible clusterable start coordinate on the current contig
     */
    int getMaxClusterableStartingPosition(final T item);

    /**
     * Compute max feasible starting position of any other item for all items in the given collection. Note the items
     * must all have the same starting contig.
     */
    default int getMaxClusterableStartingPosition(final Collection<T> items) {
        final List<String> contigA = items.stream().map(T::getContigA).distinct().collect(Collectors.toList());
        if (contigA.size() > 1) {
            throw new IllegalArgumentException("Items start on multiple contigs");
        }
        return items.stream().mapToInt(item -> getMaxClusterableStartingPosition(item)).max().getAsInt();
    }

}
