package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

@Analysis(description = "Ti/Tv Variant Evaluator")
public class TiTvVariantEvaluator extends VariantEvaluator implements StandardEval {
    @DataPoint(description = "number of transition loci", format = "%d")
    public long nTi = 0;
    @DataPoint(description = "number of transversion loci", format = "%d")
    public long nTv = 0;
    @DataPoint(description = "the transition to transversion ratio", format = "%.2f")
    public double tiTvRatio = 0.0;
    @DataPoint(description = "number of comp transition sites", format = "%d")
    public long nTiInComp = 0;
    @DataPoint(description = "number of comp transversion sites", format = "%d")
    public long nTvInComp = 0;
    @DataPoint(description = "the transition to transversion ratio for comp sites", format = "%.2f")
    public double TiTvRatioStandard = 0.0;
    @DataPoint(description = "number of derived transition loci", format = "%d")
    public long nTiDerived = 0;
    @DataPoint(description = "number of derived transversion loci", format = "%d")
    public long nTvDerived = 0;
    @DataPoint(description = "the derived transition to transversion ratio", format = "%.2f")
    public double tiTvDerivedRatio = 0.0;

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    public TiTvVariantEvaluator(VariantEvalEngine engine) {
        super(engine);
    }

    public void updateTiTv(VariantContext vc, boolean updateStandard) {
        if (vc != null && vc.isSNP() && vc.isBiallelic() && vc.isPolymorphicInSamples()) {
            if ( GATKVariantContextUtils.isTransition(vc)) {
                if (updateStandard) nTiInComp++;
                else nTi++;
            } else {                                
                if (updateStandard) nTvInComp++;
                else nTv++;
            }

            if (vc.hasAttribute("ANCESTRALALLELE")) {
                final String aaStr = vc.getAttributeAsString("ANCESTRALALLELE", "null").toUpperCase();
                if ( ! aaStr.equals(".") ) {
                    switch ( BaseUtils.SNPSubstitutionType(aaStr.getBytes()[0], vc.getAlternateAllele(0).getBases()[0] ) ) {
                        case TRANSITION: nTiDerived++; break;
                        case TRANSVERSION: nTvDerived++; break;
                        default: break;
                    }
                }
            }
        }
    }

    @Override
    public void update2(final VariantContext eval, final VariantContext comp, final VariantEvalContext context) {
        if (eval != null)
            updateTiTv(eval, false);
        if (comp != null)
            updateTiTv(comp, true);
    }

    @Override
    public void finalizeEvaluation() {
        // the ti/tv ratio needs to be set (it's not calculated per-variant).
        this.tiTvRatio = rate(nTi,nTv);
        this.tiTvDerivedRatio = rate(nTiDerived,nTvDerived);
        this.TiTvRatioStandard = rate(nTiInComp, nTvInComp);
    }
}
