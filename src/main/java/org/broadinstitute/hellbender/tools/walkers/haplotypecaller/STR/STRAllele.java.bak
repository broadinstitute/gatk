package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.Allele;

/**
 * Represents an str-allele.
 */
public final class STRAllele extends Allele {

    final STRAlleleSet set;
    final int repeatCount;
    final int index;

    STRAllele(final STRAlleleSet set, final int index, final int repeatCount) {
        super(set.composeRepeatBases(repeatCount) , set.referenceRepeatCount == repeatCount);
        this.set = set;
        this.index = index;
        this.repeatCount = repeatCount;
    }

    public final STRAlleleSet getAlleleSet() {
        return set;
    }

    public final int getRepeatCount() {
        return repeatCount;
    }
}
