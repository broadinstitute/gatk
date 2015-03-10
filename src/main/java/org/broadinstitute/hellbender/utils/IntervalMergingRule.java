package org.broadinstitute.hellbender.utils;

/**
 * a class we use to determine the merging rules for intervals passed to the GATK
 */
public enum IntervalMergingRule {
    ALL, // we merge both overlapping intervals and abutting intervals
    OVERLAPPING_ONLY // We merge intervals that are overlapping, but NOT ones that only abut each other
}
