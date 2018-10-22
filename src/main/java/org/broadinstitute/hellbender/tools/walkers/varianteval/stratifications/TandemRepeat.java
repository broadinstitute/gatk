package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.Arrays;
import java.util.List;

/**
 * Stratifies the eval RODs into sites that are tandem repeats
 */
public class TandemRepeat extends VariantStratifier {
    private final static List<Object> JUST_ALL = Arrays.asList((Object)"all");
    private final static List<Object> ALL = Arrays.asList((Object)"all", (Object)"is.repeat", (Object)"not.repeat");
    private final static List<Object> REPEAT = Arrays.asList((Object)"all", (Object)"is.repeat");
    private final static List<Object> NOT_REPEAT = Arrays.asList((Object)"all", (Object)"not.repeat");

    @Override
    public void initialize() {
        states.addAll(ALL);
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if ( eval == null || ! eval.isIndel() )
            return ALL;

        //NOTE: GATK3 ReferenceContext always used a window of 50 total bases, no matter the variant length
        byte[] bases = referenceContext.getBases(0, 50 - referenceContext.getInterval().size() + 1);
        if ( GATKVariantContextUtils.isTandemRepeat(eval, bases) ) {
            return REPEAT;
        } else {
            return NOT_REPEAT;
        }
    }
}
