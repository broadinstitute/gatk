package org.broadinstitute.hellbender.engine.dataflow;


import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

public class ReadContextData {
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
