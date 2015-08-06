package org.broadinstitute.hellbender.utils.fragments;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Useful helper class to represent the results of the reads -> fragment calculation.
 *
 * Contains singleton -- objects whose underlying reads do not overlap their mate pair
 * Contains overlappingPairs -- objects whose underlying reads do overlap their mate pair
 *
 * User: ebanks, depristo
 * Date: Jan 10, 2011
 */
public class FragmentCollection<T> {
    Collection<T> singletons;
    Collection<List<T>> overlappingPairs;

    public FragmentCollection(final Collection<T> singletons, final Collection<List<T>> overlappingPairs) {
        this.singletons = singletons == null ? Collections.<T>emptyList() : singletons;
        this.overlappingPairs = overlappingPairs == null ? Collections.<List<T>>emptyList() : overlappingPairs;
    }

    /**
     * Gets the T elements not containing overlapping elements, in no particular order
     *
     * @return
     */
    public Collection<T> getSingletonReads() {
        return singletons;
    }

    /**
     * Gets the T elements containing overlapping elements, in no particular order
     *
     * @return
     */
    public Collection<List<T>> getOverlappingPairs() {
        return overlappingPairs;
    }
}
