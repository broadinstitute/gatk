package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import org.broadinstitute.hellbender.engine.filters.VariantFilter;

public interface StructuralVariantFilter extends VariantFilter {

    /**
     * @return name of filter for use in filtered record
     */
    String getName();
}
