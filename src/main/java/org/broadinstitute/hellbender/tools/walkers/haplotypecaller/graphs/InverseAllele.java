package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import htsjdk.variant.variantcontext.Allele;

/**
 * Utility class for defining a "not" allele concept that is used to score haplotypes that are not supporting the allele.
 * In the context of AlleleFiltering we basically score an allele (set of haplotypes that support the allele) versus
 * the not allele: set of haplotypes that do not support this allele
 *
 * @author Ilya Soifer &lt;ilya.soifer@ultimagen.com
 */

public class InverseAllele extends Allele {
    final static public long serialVersionUID = 1L;

    private final Allele internalAllele;
    private final boolean referenceStatus;
    private InverseAllele(final Allele allele, boolean isReference) {
        super(allele, false);
        this.internalAllele = allele;
        referenceStatus = isReference;
    }

    // InverseAllele of inverseAllele. By definition it is the allele. In Allele filtering code we normally genotype
    // three genotypes: Hom. Allele, Het Allele/InverseAllele, Hom InverseAllele
    // Due to the way genotyping functions work one of the alleles has to be considered reference
    public static Allele of(final Allele allele, boolean refFlag){
        if (allele instanceof InverseAllele) {
            return ((InverseAllele)allele).internalAllele;
        } else {
            return new InverseAllele(allele, refFlag);
        }
    }

    public byte [] getBases(){
        return getDisplayString().getBytes();
    }
    @Override
    public boolean isReference() {
        return referenceStatus;
    }

    @Override
    public boolean isSymbolic() {
        return true;
    }

    @Override
    public String getDisplayString() {
        return "~" + super.getDisplayString();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof InverseAllele) )
            return false;

        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }

        final InverseAllele that = (InverseAllele) o;

        return internalAllele.equals(that.internalAllele);
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + internalAllele.hashCode();
        return result;
    }
}
