package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import java.io.Serializable;

/**
 * This class provides a key made of the String pair (sex genotype, contig).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SexGenotypeContigPairKey implements Serializable {

    private static final long serialVersionUID = -1710038102655907133L;

    private final String sexGenotype;
    private final String contig;

    public SexGenotypeContigPairKey(@Nonnull final String sexGenotype, @Nonnull final String contig) {
        this.sexGenotype = Utils.nonNull(sexGenotype, "The sex genotype string identifier must be non-null");
        this.contig = Utils.nonNull(contig, "The contig string identifier must be non-null");
    }

    public static SexGenotypeContigPairKey of(@Nonnull final String sexGenotype, @Nonnull final String contig) {
        return new SexGenotypeContigPairKey(sexGenotype, contig);
    }

    public String getSexGenotype() {
        return sexGenotype;
    }

    public String getContig() {
        return contig;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SexGenotypeContigPairKey that = (SexGenotypeContigPairKey) o;

        if (sexGenotype.equals(that.sexGenotype) && contig.equals(that.contig)) {
            return true;
        } else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        return 31 * sexGenotype.hashCode() + contig.hashCode();
    }
}
