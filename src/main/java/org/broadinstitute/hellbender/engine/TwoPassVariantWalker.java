package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.exceptions.GATKException;

public abstract class TwoPassVariantWalker extends MultiplePassVariantWalker {
    @Override
    final protected int numberOfPasses() { return 2; }

    @Override
    final protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        if (n == 0) {
            firstPassApply(variant, readsContext, referenceContext, featureContext);
        } else if (n == 1) {
            secondPassApply(variant, readsContext, referenceContext, featureContext);
        } else {
            throw new GATKException.ShouldNeverReachHereException("This two-pass walker should never reach (zero-indexed) pass " + n);
        }
    }


    @Override
    final protected void afterNthPass(final int n) {
        if (n == 0) {
            afterFirstPass();
        } else if (n > 1) {
            throw new GATKException.ShouldNeverReachHereException("This two-pass walker should never reach (zero-indexed) pass " + n);
        }
    }

    protected abstract void firstPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext);

    /**
     * Process the data collected during the first pass. This method is called between the two traversals
     */
    protected abstract void afterFirstPass();

    protected abstract void secondPassApply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext);

}
