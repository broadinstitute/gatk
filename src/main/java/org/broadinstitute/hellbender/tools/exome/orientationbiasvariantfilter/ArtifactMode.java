package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import java.io.Serializable;

public class ArtifactMode implements Serializable {

    static final long serialVersionUID = 3373373374112L;

    private byte ref;
    private byte alt;

    public ArtifactMode(final byte ref, final byte alt){
        this.ref = ref;
        this.alt = alt;

    }

    public byte getRef() {
        return ref;
    }

    public void setRef(byte ref) {
        this.ref = ref;
    }

    public byte getAlt() {
        return alt;
    }

    public void setAlt(byte alt) {
        this.alt = alt;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        final ArtifactMode that = (ArtifactMode) obj;
        return (this.getRef() == that.getRef()) && (this.getAlt() == (that.getAlt()));
    }
}
