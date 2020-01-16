package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.stream.Collectors;

public class ReadAndHaplotype {
    final GATKRead read;
    final Haplotype haplotype;

    public ReadAndHaplotype(final GATKRead read, final Haplotype haplotype){
        this.read = read;
        this.haplotype = haplotype;
    }

    public GATKRead getRead() {
        return read;
    }

    public Haplotype getHaplotype() {
        return haplotype;
    }

    public static List<Haplotype> getHaplotypeList(final List<ReadAndHaplotype> list){
        return list.stream().map(rah -> rah.getHaplotype()).collect(Collectors.toList());
    }
}
