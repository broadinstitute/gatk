package org.broadinstitute.hellbender.tools.gvs.common;

public enum GQStateEnum {
    VARIANT("v", null, 7),
    STAR("*", null, 8),
    ZERO("0", 0, 0),
    TEN("1", 10, 1),
    TWENTY("2", 20, 2),
    THIRTY("3", 30, 3),
    FORTY("4", 40, 4),
    FIFTY("5", 50, 5),
    SIXTY("6", 60, 6),
    MISSING("m", null, 9),
    UNKNOWN("u", null, 10),
    NONE("");

    String value;
    Integer referenceGQ;
    Integer compressedValue;

    GQStateEnum(String value) {
        this(value, null);
    }

    GQStateEnum(String value, Integer referenceGQ) {
        this(value, referenceGQ, null);
    }

    GQStateEnum(String value, Integer referenceGQ, Integer compressedValue) {
        this.value = value;
        this.referenceGQ = referenceGQ;
        this.compressedValue = compressedValue;
    }

    public String getValue() {
        return value;
    }

    public Integer getReferenceGQ() {
        return referenceGQ;
    }

    public Integer getCompressedValue() { return compressedValue; }

}
