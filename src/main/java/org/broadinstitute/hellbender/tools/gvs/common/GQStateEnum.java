package org.broadinstitute.hellbender.tools.gvs.common;

public enum GQStateEnum {
    VARIANT("v"),
    STAR("*"),
    ZERO("0", 0),
    TEN("1", 10),
    TWENTY("2", 20),
    THIRTY("3", 30),
    FORTY("4", 40),
    FIFTY("5", 50),
    SIXTY("6", 60),
    MISSING("m"),
    UNKNOWN("u"),
    NONE("");

    String value;
    Integer referenceGQ;

    GQStateEnum(String value) {
        this(value, null);
    }

    GQStateEnum(String value, Integer referenceGQ) {
        this.value = value;
        this.referenceGQ = referenceGQ;
    }

    public String getValue() {
        return value;
    }

    public Integer getReferenceGQ() {
        return referenceGQ;
    }
}
