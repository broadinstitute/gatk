package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

@CommandLineProgramProperties(
        usage = "Walks over the input data set, calculating the number of variants seen.",
        usageShort = "Count variants",
        programGroup = VariantProgramGroup.class
)
public final class CountVariants extends VariantWalker{
    private long count;

    @Override
    public void onTraversalStart() {
        count = 0;
    }

    @Override
    public void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        count++;
    }

    @Override
    public Long onTraversalDone() {
        return count;
    }
}
