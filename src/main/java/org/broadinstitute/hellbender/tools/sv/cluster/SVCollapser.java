package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVLocatable;

public interface SVCollapser<T extends SVLocatable, C extends BasicOutputCluster<T>> {

    /**
     * Constructs a single record representative of the given variants.
     * @param items cluster to collapse
     * @return a call approximating the average event of the input
     */
    T collapse(C items);

}
