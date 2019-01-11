package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;

@Analysis(name = "PrintMissingComp", description = "count the number of comp SNP sites that are not in eval")
public class PrintMissingComp extends VariantEvaluator {
    @DataPoint(description = "number of comp SNP sites outside of eval sites", format = "%d")
    public long nMissing = 0;

    public String getName() {
        return "PrintMissingComp";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public void update2(VariantContext eval, VariantContext comp, final ReferenceContext referenceContext, final ReadsContext readsContext, final FeatureContext featureContext) {
        final boolean compIsGood = comp != null && comp.isNotFiltered() && comp.isSNP();
        final boolean evalIsGood = eval != null && eval.isSNP();

        if ( compIsGood && !evalIsGood ) {
            nMissing++;
        }
    }
}