package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.Iterator;

/**
 * Created by bimber on 5/17/2017.
 */
@Analysis(description = "Counts called and filtered genotypes across samples")
public class GenotypeFilterSummary extends VariantEvaluator {
    public GenotypeFilterSummary(VariantEvalEngine engine) {
        super(engine);
    }

    @DataPoint(description = "Number of Called Genotypes", format = "%d")
    public long nCalledNotFiltered = 0;

    @DataPoint(description = "Number of No-Call Genotypes", format = "%d")
    public long nNoCallOrFiltered = 0;

    @Override
    public int getComparisonOrder() {
        return 1;
    }

    @Override
    public void update1(final VariantContext eval, final VariantEvalContext context) {
        Iterator<Genotype> it = eval.getGenotypes().iterator();
        while (it.hasNext()){
            Genotype g = it.next();
            if (g.isCalled() && !g.isFiltered()){
                nCalledNotFiltered++;
            }
            else if (g.isNoCall() || g.isFiltered()){
                nNoCallOrFiltered++;
            }
        }
    }
}
