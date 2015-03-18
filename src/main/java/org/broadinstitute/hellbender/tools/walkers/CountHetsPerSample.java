package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import java.util.HashMap;

@CommandLineProgramProperties(
        usage = "Count number of hets per sample.",
        usageShort = "Count hets per sample",
        programGroup = VariantProgramGroup.class
)

public final class CountHetsPerSample extends VariantWalker {
    private HashMap<String, Long> sampleCounts = new HashMap<>(); // FIXME when do member initializers run?

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        for (final Genotype g : variant.getGenotypes()) {
            if (g.isHet()) {
                sampleCounts.merge(g.getSampleName(), 1L, (v1, v2) -> v1 + v2);
            }
        }
    }

    @Override
    public Object onTraversalDone() { return sampleCounts; }
}
