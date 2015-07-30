package org.broadinstitute.hellbender.utils.samples;

/**
 * Categorical sample trait for association and analysis
 *
 * Samples can have unknown status, be affected or unaffected by the
 * categorical trait, or they can be marked as actually having an
 * other trait value (stored in an associated value in the Sample class)
 *
 * @author Mark DePristo
 * @since Sept. 2011
 */
public enum Affection {
    /** Status is unknown */
    UNKNOWN,
    /** Suffers from the disease */
    AFFECTED,
    /** Unaffected by the disease */
    UNAFFECTED,
    /** An "other" trait: value of the trait is stored elsewhere and is an arbitrary string */
    OTHER
}
