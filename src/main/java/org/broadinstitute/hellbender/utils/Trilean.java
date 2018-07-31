package org.broadinstitute.hellbender.utils;

/**
 * An enumeration to represent true, false, or unknown.
 */
public enum Trilean {
    TRUE, FALSE, UNKNOWN;

    public static Trilean of (final boolean booleanValue) {
        return booleanValue ? Trilean.TRUE : Trilean.FALSE;
    }
}