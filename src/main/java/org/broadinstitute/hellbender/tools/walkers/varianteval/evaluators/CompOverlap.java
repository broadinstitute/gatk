package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;

@Analysis(description = "The overlap between eval and comp sites")
public class CompOverlap extends VariantEvaluator implements StandardEval {
    @DataPoint(description = "number of eval variant sites", format = "%d")
    public long nEvalVariants = 0;

    @DataPoint(description = "number of eval sites outside of comp sites", format = "%d")
    public long novelSites = 0;

    @DataPoint(description = "number of eval sites at comp sites", format = "%d")
    public long nVariantsAtComp = 0;

    @DataPoint(description = "percentage of eval sites at comp sites", format = "%.2f" )
    public double compRate = 0.0;

    @DataPoint(description = "number of concordant sites", format = "%d")
    public long nConcordant = 0;

    @DataPoint(description = "the concordance rate", format = "%.2f")
    public double concordantRate = 0.0;

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public long nNovelSites() { return nEvalVariants - nVariantsAtComp; }
    public double compRate() { return rate(nVariantsAtComp, nEvalVariants); }
    public double concordanceRate() { return rate(nConcordant, nVariantsAtComp); }

    public void finalizeEvaluation() {
        compRate = 100 * compRate();
        concordantRate = 100 * concordanceRate();
        novelSites = nNovelSites();
    }

    /**
     * Returns true if every allele in eval is also in comp
     *
     * @param eval  eval context
     * @param comp db context
     * @return true if eval and db are discordant
     */
    public boolean discordantP(VariantContext eval, VariantContext comp) {
        for (Allele a : eval.getAlleles()) {
            if (!comp.hasAllele(a, true))
                return true;
        }

        return false;
    }

    public void update2(VariantContext eval, VariantContext comp, final ReferenceContext referenceContext, final ReadsContext readsContext, final FeatureContext featureContext) {
        boolean evalIsGood = eval != null && eval.isPolymorphicInSamples();
        boolean compIsGood = comp != null && comp.isNotFiltered();

        if (evalIsGood) nEvalVariants++;           // count the number of eval events

        if (compIsGood && evalIsGood) {
            nVariantsAtComp++;

            if (!discordantP(eval, comp)) {    // count whether we're concordant or not with the comp value
                nConcordant++;
            }
        }
    }
}
