package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

@Analysis(name = "PrintMissingComp", description = "count the number of comp SNP sites that are not in eval")
public class PrintMissingComp extends VariantEvaluator {
    public PrintMissingComp(VariantEvalEngine engine) {
        super(engine);
    }

    @DataPoint(description = "number of comp SNP sites outside of eval sites", format = "%d")
    public long nMissing = 0;

    public String getName() {
        return "PrintMissingComp";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    @Override
    public void update2(final VariantContext eval, final VariantContext comp, final VariantEvalContext context) {
        final boolean compIsGood = comp != null && comp.isNotFiltered() && comp.isSNP();
        final boolean evalIsGood = eval != null && eval.isSNP();

        if ( compIsGood && !evalIsGood ) {
            nMissing++;
        }
    }
}