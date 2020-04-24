package org.broadinstitute.hellbender.engine;

@FunctionalInterface
interface IntervalComparator2 {
    int compare(final int start1, final int end1, final int start2, final int end2);
}
