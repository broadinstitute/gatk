package org.broadinstitute.hellbender.tools.variantdb;

public class ProbeInfo {
    long probeId;

    String contig;
    long position;
    String ref;
    String alleleA;
    String alleleB;
    String name;

    public ProbeInfo(long probeId, String name, String contig, long position, String ref, String alleleA, String alleleB) {
        this.probeId = probeId;
        this.name = name;
        this.contig = contig;
        this.position = position;
        this.ref = ref;
        this.alleleA = alleleA;
        this.alleleB = alleleB;     
    }

}
