package org.broadinstitute.hellbender.utils.functional;

/**
 * Created by davidben on 8/19/16.
 */
@FunctionalInterface
public interface IntBiConsumer {
    void accept(final int alleleIndex, final int alleleCount);
}
