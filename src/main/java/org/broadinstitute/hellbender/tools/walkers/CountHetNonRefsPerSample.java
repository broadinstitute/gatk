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
import java.util.Map;

@CommandLineProgramProperties(
        usage = "Count number of het non-refs per sample.",
        usageShort = "Count het non-refs per sample",
        programGroup = VariantProgramGroup.class
)
public final class CountHetNonRefsPerSample extends VariantWalker {
    private final Map<String, Long> sampleCounts = new HashMap<>();

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        variant.getGenotypes().stream()
                .filter(Genotype::isAvailable)
                .forEach(g -> sampleCounts.merge(g.getSampleName(), 1L, (v1, v2) -> v1 + v2) );
    }

    @Override
    public Object onTraversalDone() {
        return sampleCounts;
    }
}
