package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;

import java.util.Collection;

public interface SVCollapser<T extends SVLocatable> {

    /**
     * Constructs a single record representative of the given variants.
     * @param items events to collapse
     * @return a call approximating the average event of the input
     */
    T collapse(final Collection<T> items);
}
