package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

@CommandLineProgramProperties(
        usage = "Calculate the proportion of samples that have missing data (and hence no genotype) for a given site.",
        usageShort = "Calculate proportion of samples missing genotype for a given site",
        programGroup = VariantProgramGroup.class
)
public final class PercentMissing extends VariantWalker {

    @Argument(fullName = "minimum GQ", shortName = "minGQ", doc="Variants with QG less than minGQ are considered missing")
    public int minGQ; //default should be 0, but would be nice to make this explicit

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        Long numMissing = variant.getGenotypes().stream()
                                 .filter(g -> isMissing(g))
                                 .count();
        double percentMissing = (double) numMissing / variant.getNSamples();
        System.out.println(percentMissing);
    }

    private boolean isMissing(Genotype g) {
        return g.isNoCall() || g.getGQ() < minGQ;
    }


    @Override
    public Object onTraversalDone() {
        return "WORKED!";
    }
}
