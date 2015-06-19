package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;

/**
 * Minimal interface for random access to a collection of Alleles.
 */
public interface AlleleList<A extends Allele> {

    public int alleleCount();

    public int alleleIndex(final A allele);

    public A alleleAt(final int index);
}
