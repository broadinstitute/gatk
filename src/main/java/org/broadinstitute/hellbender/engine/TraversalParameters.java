package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;

/**
 * A simple container class for parameters controlling which records get returned during traversals.
 * Holds a List of intervals (may be empty) and a flag controlling whether unmapped records should be
 * returned.
 */
public class TraversalParameters implements Serializable{

    private static final long serialVersionUID = 1l;

    private final List<SimpleInterval> intervalsForTraversal;
    private final boolean traverseUnmapped;

    /**
     * @param intervalsForTraversal List of intervals for traversal. Only records overlapping these intervals will be returned.
     * @param traverseUnmapped True if unmapped records should be traversed
     */
    public TraversalParameters( final List<SimpleInterval> intervalsForTraversal, final boolean traverseUnmapped ) {
        this.intervalsForTraversal = intervalsForTraversal != null ? intervalsForTraversal : Collections.emptyList();
        this.traverseUnmapped = traverseUnmapped;
    }

    /**
     * @return List of intervals for traversal. Only records overlapping these intervals will be returned.
     */
    public List<SimpleInterval> getIntervalsForTraversal() {
        return Collections.unmodifiableList(intervalsForTraversal);
    }

    /**
     * @return True if unmapped records should be traversed
     */
    public boolean traverseUnmappedReads() {
        return traverseUnmapped;
    }
}
