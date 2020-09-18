package org.broadinstitute.hellbender.tools.variantdb.arrays;

public enum GT_encoding {
    HOM_REF("R"),
    HET0_1("X"),
    HOM_VAR("A"),
    HET1_2("Y"),
    HOM_ALT2("B"),
    MISSING(".");

    String value;
    GT_encoding(String v) {
        value = v;
    }
    String getValue() {
        return value;
    }

    public static GT_encoding getGTEncodingFromValue(String value) {
        GT_encoding response = MISSING;
        if (value != null) {
            if (value.equalsIgnoreCase(HOM_REF.value)) {
                response = HOM_REF;
            } else if (value.equalsIgnoreCase(HET0_1.value)) {
                response = HET0_1;
            } else if (value.equalsIgnoreCase(HOM_VAR.value)) {
                response = HOM_VAR;
            } else if (value.equalsIgnoreCase(HET1_2.value)) {
                response = HET1_2;
            } else if (value.equalsIgnoreCase(HOM_ALT2.value)) {
                response = HOM_ALT2;
            } else if (value.equalsIgnoreCase(MISSING.value)) {
                response = MISSING;
            }
        }
        return response;
    }
}
