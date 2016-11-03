package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.iterators.PushToPullIterator;

import java.util.*;

/**
 * Turns an iterator of {@link VariantContext} into one which merges GVCF blocks.
 */
public class GVCFBlockMergingIterator extends PushToPullIterator<VariantContext> {

    public GVCFBlockMergingIterator(Iterator<VariantContext> variants, final List<Integer> gqPartitions, final int defaultPloidy){
       super(variants, new GVCFBlockCombiner(gqPartitions, defaultPloidy));
    }

}
