package org.broadinstitute.hellbender.utils.genotyper;

import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents the valid genotypes of a diploid sample.
 */
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

    public final byte base1;
    public final byte base2;

    private static final DiploidGenotype[][] conversionMatrix = {
            { AA, AC, AG, AT },
            { AC, CC, CG, CT },
            { AG, CG, GG, GT },
            { AT, CT, GT, TT }
    };

    private DiploidGenotype(final char base1, final char base2) {
        this((byte)base1, (byte)base2);
    }

    private DiploidGenotype(final byte base1, final byte base2) {
        this.base1 = base1;
        this.base2 = base2;
    }

    public boolean isHomRef(final byte ref) {
        return isHom() && ref == base1;
    }

    public boolean isHomVar(final byte ref) {
        return isHom() && ref != base1;
    }

    public boolean isHetRef(final byte ref) {
        if ( base1 == ref ) {
            return ref != base2;
        } else {
            return base2 == ref;
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
        Utils.validIndex(baseIndex1, values().length);
        Utils.validIndex(baseIndex2, values().length);
        return conversionMatrix[baseIndex1][baseIndex2];
    }
}