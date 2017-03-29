package org.broadinstitute.hellbender.utils;

public final class GenotypeCounts {

    private final double ref;
    private final double het;
    private final double hom;

    public GenotypeCounts(final double ref, final double het, final double hom ){
        this.ref = ref;
        this.het = het;
        this.hom = hom;
    }

    public double getRefs() {
        return ref;
    }

    public double getHets() {
        return het;
    }

    public double getHoms() {
        return hom;
    }

}
