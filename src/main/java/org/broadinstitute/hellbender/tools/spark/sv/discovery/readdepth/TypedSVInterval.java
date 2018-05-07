package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

public final class TypedSVInterval {
    private final SVInterval interval;
    private final StructuralVariantType type;

    public TypedSVInterval(final SVInterval interval, final StructuralVariantType type) {
        this.interval = interval;
        this.type = type;
    }

    public int getContig() {
        return interval.getContig();
    }

    public int getLength() {
        return interval.getLength();
    }

    public int getStart() {
        return interval.getStart();
    }

    public int getEnd() {
        return interval.getEnd();
    }

    public StructuralVariantType getType() {
        return type;
    }
}
