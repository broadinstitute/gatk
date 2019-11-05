package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.filters.UMIReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class UMI {
    final String umi1;
    final String umi2;

    public UMI(final GATKRead read){
        final String readUMI = read.getAttributeAsString(UMIReadFilter.UMI_TAG);
        this.umi1 = readUMI.split("-", 2)[0];
        this.umi2 = readUMI.split("-", 2)[1];
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

    public boolean equals(final UMI that) {
        return true;
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
