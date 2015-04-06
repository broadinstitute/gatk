package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Jon on 4/6/15.
 */

@CommandLineProgramProperties(
        usage = "Compute the MLE for het dispersion.",
        usageShort = "Compute MLE for het dispersion",
        programGroup = VariantProgramGroup.class
)
public final class HetDispersionPerSample extends VariantWalker {
    private final Map<String, Double> sampleMLEs = new HashMap<>();
    private final Map<String, ArrayList<Pair<Integer, Integer>>> sampleDepthAlts = new HashMap<>();

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (!variant.isBiallelic())
            return;
        variant.getGenotypes().stream()
                .filter(Genotype::isHet)
                .forEach(g -> {
                    ArrayList<Pair<Integer, Integer>> gDepthAlts = sampleDepthAlts.computeIfAbsent(g.getSampleName(),
                            s -> new ArrayList<>());
                    gDepthAlts.add(Pair.of(g.getDP(),
                            g.getAD()[1]));
                });
    }

    @Override
    public Object onTraversalDone() {
        return sampleDepthAlts;
    }
}
