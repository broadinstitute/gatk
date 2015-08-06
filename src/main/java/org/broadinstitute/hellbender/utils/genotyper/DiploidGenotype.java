package org.broadinstitute.hellbender.utils.genotyper;

import org.broadinstitute.hellbender.utils.BaseUtils;

public enum DiploidGenotype {
    AA ('A', 'A'),
    AC ('A', 'C'),
    CC ('C', 'C'),
    AG ('A', 'G'),
    CG ('C', 'G'),
    GG ('G', 'G'),
    AT ('A', 'T'),
    CT ('C', 'T'),
    GT ('G', 'T'),
    TT ('T', 'T');

    public byte base1, base2;

    @Deprecated
    private DiploidGenotype(final char base1, final char base2) {
        this((byte)base1, (byte)base2);
    }

    private DiploidGenotype(final byte base1, final byte base2) {
        this.base1 = base1;
        this.base2 = base2;
    }

    public boolean isHomRef(final byte r) {
        return isHom() && r == base1;
    }

    public boolean isHomVar(final byte r) {
        return isHom() && r != base1;
    }

    public boolean isHetRef(final byte r) {
        if ( base1 == r ) {
            return r != base2;
        } else {
            return base2 == r;
        }
    }

    public boolean isHom() {
        return ! isHet();
    }

    public boolean isHet() {
        return base1 != base2;
    }

    /**
     * create a diploid genotype, given a character to make into a hom genotype
     * @param hom the character to turn into a hom genotype, i.e. if it is A, then returned will be AA
     * @return the diploid genotype
     */
    public static DiploidGenotype createHomGenotype(final byte hom) {
        final int index = BaseUtils.simpleBaseToBaseIndex(hom);
        if ( index == -1 ) {
            throw new IllegalArgumentException(hom + " is not a valid base character");
        }
        return conversionMatrix[index][index];
    }

    /**
     * create a diploid genotype, given 2 chars which may not necessarily be ordered correctly
     * @param base1 base1
     * @param base2 base2
     * @return the diploid genotype
     */
    public static DiploidGenotype createDiploidGenotype(final byte base1, final byte base2) {
        final int index1 = BaseUtils.simpleBaseToBaseIndex(base1);
        if ( index1 == -1 ) {
            throw new IllegalArgumentException(base1 + " is not a valid base character");
        }
        final int index2 = BaseUtils.simpleBaseToBaseIndex(base2);
        if ( index2 == -1 ) {
            throw new IllegalArgumentException(base2 + " is not a valid base character");
        }
        return conversionMatrix[index1][index2];
    }

    /**
     * create a diploid genotype, given 2 base indexes which may not necessarily be ordered correctly
     * @param baseIndex1 base1
     * @param baseIndex2 base2
     * @return the diploid genotype
     */
    public static DiploidGenotype createDiploidGenotype(final int baseIndex1, final int baseIndex2) {
        if ( baseIndex1 == -1 ) {
            throw new IllegalArgumentException(baseIndex1 + " does not represent a valid base character");
        }
        if ( baseIndex2 == -1 ) {
            throw new IllegalArgumentException(baseIndex2 + " does not represent a valid base character");
        }
        return conversionMatrix[baseIndex1][baseIndex2];
    }

    private static final DiploidGenotype[][] conversionMatrix = {
            { DiploidGenotype.AA, DiploidGenotype.AC, DiploidGenotype.AG, DiploidGenotype.AT },
            { DiploidGenotype.AC, DiploidGenotype.CC, DiploidGenotype.CG, DiploidGenotype.CT },
            { DiploidGenotype.AG, DiploidGenotype.CG, DiploidGenotype.GG, DiploidGenotype.GT },
            { DiploidGenotype.AT, DiploidGenotype.CT, DiploidGenotype.GT, DiploidGenotype.TT }
    };
}