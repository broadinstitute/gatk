package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.VariantWalker;

@CommandLineProgramProperties(
        usage = "Walks over the input data set, calculating the number of variants seen.",
        usageShort = "Count variants",
        programGroup = VariantProgramGroup.class
)

public final class CountVariants extends VariantWalker{
    private long count = 0;

    @Override
    public void apply(VariantContext variant) {
        count++;
    }

    @Override
    public Object onTraversalDone() {
        return count;
    }
}
