package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Collections;
import java.util.List;

/**
 * Stratifies the eval RODs by the allele frequency of the alternate allele
 *
 * Either uses a constant 0.005 frequency grid, and projects the AF INFO field value or logit scale from -30 to 30.
 * Requires that AF be present in every ROD, otherwise this stratification throws an exception
 */
public class AlleleFrequency extends VariantStratifier {
    public enum StratifyingScale {
        LINEAR,
        LOGARITHMIC
    }

    private static final int LOGIT_BIN_WIDTH = 1;

    private static int logLimit = 30; // go from -30 to 30 using logit function

    public AlleleFrequency(VariantEvalEngine engine) {
        super(engine);

        switch (getEngine().getAFScale()) {
            case LINEAR:
                for (double a = 0.000; a <= 1.005; a += 0.005) {
                    states.add(String.format("%.3f", a));
                }
                break;
            case LOGARITHMIC:
                for (int a = -logLimit; a <= logLimit; a += LOGIT_BIN_WIDTH) {
                    states.add(String.format("%d", a));
                }

        }
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null) {
            try {
                double alleleFrequency = Collections.max(eval.getAttributeAsDoubleList("AF", 0.0));
                if (getEngine().getCompAFStratifier()) {
                    if (comp != null) {
                        alleleFrequency = Collections.max(comp.getAttributeAsDoubleList("AF", 0.0));
                    } else {
                        alleleFrequency = 0.0; // any site that isn't in the comp should be the expected lowest allele fraction
                    }
                }
                switch(getEngine().getAFScale()) {
                    case LINEAR:
                        return Collections.singletonList(String.format("%.3f", (5.0 * MathUtils.roundToNDecimalPlaces(alleleFrequency / 5.0, 3))));
                    case LOGARITHMIC:
                        return Collections.singletonList(String.format("%d", getLogitBucket(alleleFrequency + Math.pow(10.0, -6))));
                }

            } catch (Exception e) {
                return Collections.emptyList();
            }
        }

        return Collections.emptyList();
    }

    private int getLogitBucket(double alleleFrequency) {
        // calculate logit value
        Float score = (float)( -10 * Math.log10((alleleFrequency/(1-alleleFrequency))));

        // find correct bucket
        return Math.min(logLimit, Math.max(-logLimit, Math.round(score)));
    }
}

