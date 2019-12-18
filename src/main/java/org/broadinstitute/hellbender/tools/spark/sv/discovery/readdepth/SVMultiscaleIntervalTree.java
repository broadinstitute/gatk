package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class SVMultiscaleIntervalTree<V> {
    private final List<SVIntervalTree<V>> trees;

    public SVMultiscaleIntervalTree(final List<SVIntervalTree<V>> trees) {
        Utils.nonNull(trees, "Tree collection cannot be null");
        this.trees = trees;
    }

    public Collection<SVIntervalTree.Entry<V>> overlappers(final SVInterval interval ) {
        int maxOverlap = 0;
        Collection<SVIntervalTree.Entry<V>> maxOverlappers = Collections.emptyList();
        for (final SVIntervalTree<V> tree : trees) {
            final Collection<SVIntervalTree.Entry<V>> overlappers = Lists.newArrayList(tree.overlappers(interval));
            int sum = 0;
            for (final SVIntervalTree.Entry<V> entry : overlappers) {
                sum += entry.getInterval().overlapLen(interval);
            }
            if (sum == interval.getLength()) {
                return overlappers;
            }
            if (sum > maxOverlap) {
                maxOverlap = sum;
                maxOverlappers = overlappers;
            }
        }
        return maxOverlappers;
    }

    public boolean hasOverlapper(final SVInterval interval) {
        for (final SVIntervalTree<V> tree : trees) {
            if (tree.hasOverlapper(interval)) return true;
        }
        return false;
    }
}
