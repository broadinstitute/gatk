package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.filters.UMIReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class UMI {
    // The umis are stored here in the order they are stored in the read.
    // e.g. if the read UMI is AAT-CGT, then umi1 = AAT, umi2 = CGT
    public final String umi1;
    public final String umi2;

    // Lexicographical ordering of the two (duplex) UMIs;
    public String umiSmall;
    public String umiLarge;

    public UMI(final GATKRead read){
        final String readUMI = read.getAttributeAsString(UMIReadFilter.UMI_TAG);
        this.umi1 = readUMI.split("-", 2)[0];
        this.umi2 = readUMI.split("-", 2)[1];
        if (this.umi1.compareTo(this.umi2) > 0) {
            umiLarge = this.umi1;
            umiSmall = this.umi2;
        } else {
            umiLarge = this.umi2;
            umiSmall = this.umi1;
        }

    }

    public static Pair<String, String> getUMI(final GATKRead read){
        final String readUMI = read.getAttributeAsString(UMIReadFilter.UMI_TAG);
        final String umi1 = readUMI.split("-", 2)[0];
        final String umi2 = readUMI.split("-", 2)[1];

        return new ImmutablePair<>(umi1, umi2);
    }

    public boolean equalsReadUMI(final GATKRead read){
        return false;
    }

    // Strict match e.g. equalsExactly("AGT-GCT", "GCT-AGT") returns false
    public boolean equalsExactly(final UMI that) {
        return that.umi1.equals(this.umi1) && that.umi2.equals(this.umi2);
    }

    // Check for whether the two umis came from the same molecule
    // e.g. equalsModuloOrder("AGT-GCT", "GCT-AGT") returns true
    //      equalsModuloOrder("AGT-GCT", "AGT-GCT") also returns true
    public boolean equalsModuloOrder(final UMI that) {
        // Sort the respective UMIs in lexicographical order and compare
        return that.umiSmall.equals(this.umiSmall) && that.umiLarge.equals(this.umiLarge);
    }

    public String getStandardizedUMI(){
        if (umi1.compareTo(umi2) > 0){
            // umi1 is lexicographically greater e.g. umi1 = TAT, umi2 = AAC
            return umi2 + "-" + umi1;
        } else {
            return umi1 + "-" + umi2;
        }
    }

}