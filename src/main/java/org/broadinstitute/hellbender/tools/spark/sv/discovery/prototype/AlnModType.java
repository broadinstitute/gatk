package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

public enum AlnModType {
    NONE, UNDERGONE_OVERLAP_REMOVAL, EXTRACTED_FROM_LARGER_ALIGNMENT, FROM_SPLIT_GAPPED_ALIGNMENT;

    public enum ModTypeString {
        O, H, E, S;
    }
    @Override
    public String toString() {
        switch (this) {
            case NONE: return ModTypeString.O.name();
            case UNDERGONE_OVERLAP_REMOVAL: return ModTypeString.H.name();
            case EXTRACTED_FROM_LARGER_ALIGNMENT: return ModTypeString.E.name();
            case FROM_SPLIT_GAPPED_ALIGNMENT: return ModTypeString.S.name();
            default: throw new IllegalArgumentException();
        }
    }
}
