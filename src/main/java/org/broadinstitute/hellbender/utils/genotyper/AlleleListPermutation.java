package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.collections.Permutation;

/**
 * Marks allele list permutation implementation classes.
 */
public interface AlleleListPermutation<A extends Allele> extends Permutation<A>, AlleleList<A> {
}