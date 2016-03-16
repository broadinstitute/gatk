package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

@CommandLineProgramProperties(
        summary = "Walks over the input data set, calculating the number of variants seen.",
        oneLineSummary = "Count variants in a VCF file",
        programGroup = VariantProgramGroup.class
)
public final class CountVariants extends VariantWalker{
    private long count = 0;

    @Override
    public void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext ) {
        count++;
    }

    @Override
    public Object onTraversalSuccess() {
        return count;
    }
}
