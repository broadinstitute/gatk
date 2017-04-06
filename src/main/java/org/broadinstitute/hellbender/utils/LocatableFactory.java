package org.broadinstitute.hellbender.utils;

/**
 * Common interface for Locatable factories.
 * <p>Creates locatables from its contig, start and stop locations</p>
 * @param <L> the locatable type.
 */
@FunctionalInterface
public interface LocatableFactory<L> {
    L create(final String contig, final int start, final int stop);
}
