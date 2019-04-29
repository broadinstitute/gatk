package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Collections;
import java.util.List;

/**
 * Stratifies the eval RODs by the allele frequency of the alternate allele
 *
 * Uses a constant 0.005 frequency grid, and projects the AF INFO field value.  Requires
 * that AF be present in every ROD, otherwise this stratification throws an exception
 */
@DocumentedFeature(groupName= "tools", summary = "Stratify by  eval RODs by the allele frequency of the alternate allele") // , extraDocs = {VariantEval.class})
public class AlleleFrequency extends VariantStratifier {

    public enum StratifyingScale {
        LINEAR,
        LOGARITHMIC
    }

    private StratifyingScale scale;
    private boolean useCompAFStratifier;

    private static int log_lim = 30; // go form -30 to 30

    @Override
    public void initialize() {
        scale = getVariantEvalWalker().getAFScale();
        useCompAFStratifier = getVariantEvalWalker().getCompAFStratifier();

        switch (scale) {
            case LINEAR:
                for (double a = 0.000; a <= 1.005; a += 0.005) {
                    states.add(String.format("%.3f", a));
                }
                break;
            case LOGARITHMIC:
                for (int a = -log_lim; a <= log_lim; a += 1) {
                    states.add(String.format("%d", a));
                }

        }
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval != null) {
            try {
                Double allele_fraction  = Collections.max(eval.getAttributeAsDoubleList("AF", 0.0));
                if (useCompAFStratifier) {
                    if (comp != null) {
                        allele_fraction = Collections.max(comp.getAttributeAsDoubleList("AF", 0.0));
                    } else {
                        allele_fraction = 0.0; // any site that isn't in the comp should be the expected lowest allele fraction
                    }
                }
                switch(scale) {
                    case LINEAR:
                        return Collections.singletonList((Object)String.format("%.3f", (5.0 * MathUtils.roundToNDecimalPlaces(allele_fraction / 5.0, 3))));
                    case LOGARITHMIC:
                        allele_fraction += Math.pow(10.0, -6); // never have AF of 0.0
                        // using logit function
                        Float score = (float)( -10 * Math.log10((allele_fraction/(1-allele_fraction))));
                        return Collections.singletonList(String.format("%d", Math.min(log_lim, Math.max(-log_lim, Math.round(score)))));
                }

            } catch (Exception e) {
                return Collections.emptyList();
            }
        }

        return Collections.emptyList();
    }
}
