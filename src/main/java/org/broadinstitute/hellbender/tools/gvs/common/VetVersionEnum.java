package org.broadinstitute.hellbender.tools.gvs.common;

public enum VetVersionEnum {
    V1(1),
    V2(2);

    int index;

    VetVersionEnum(int index) {
        this.index = index;
    }
}
