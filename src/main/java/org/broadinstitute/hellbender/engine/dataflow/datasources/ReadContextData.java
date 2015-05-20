package org.broadinstitute.hellbender.engine.dataflow.datasources;


import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.Serializable;

public class ReadContextData implements Serializable {
    private final ReferenceBases referenceBases;
    private final Iterable<Variant> variants;

    public ReadContextData( final ReferenceBases referenceBases, final Iterable<Variant> variants ) {
        this.referenceBases = referenceBases;
        this.variants = variants;
    }

    public ReferenceBases getOverlappingReferenceBases() {
        return referenceBases;
    }

    public Iterable<Variant> getOverlappingVariants() {
        return variants;
    }
}
