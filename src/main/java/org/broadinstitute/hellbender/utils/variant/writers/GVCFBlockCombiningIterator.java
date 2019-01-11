package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.iterators.PushToPullIterator;

import java.util.*;

/**
 * Turns an iterator of {@link VariantContext} into one which combines GVCF blocks.
 */
public class GVCFBlockCombiningIterator extends PushToPullIterator<VariantContext> {

    public GVCFBlockCombiningIterator(Iterator<VariantContext> variants, final List<Number> gqPartitions, final int defaultPloidy){
       super(variants, new GVCFBlockCombiner(gqPartitions, defaultPloidy));
    }

}
