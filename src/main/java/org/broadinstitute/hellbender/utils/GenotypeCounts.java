package org.broadinstitute.hellbender.utils;

public final class GenotypeCounts {

    private final int ref;
    private final int het;
    private final int hom;

    public GenotypeCounts(final int ref, final int het, final int hom ){
        this.ref = ref;
        this.het = het;
        this.hom = hom;
    }

    public int getRefs() {
        return ref;
    }

    public int getHets() {
        return het;
    }

    public int getHoms() {
        return hom;
    }

}
