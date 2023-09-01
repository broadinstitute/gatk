package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.Arrays;
import java.util.List;

/**
 * Stratifies the evals into sites that are tandem repeats
 */
public class TandemRepeat extends VariantStratifier {
    private final static List<Object> ALL = Arrays.asList((Object)"all", (Object)"is.repeat", (Object)"not.repeat");
    private final static List<Object> REPEAT = Arrays.asList((Object)"all", (Object)"is.repeat");
    private final static List<Object> NOT_REPEAT = Arrays.asList((Object)"all", (Object)"not.repeat");

    public TandemRepeat(VariantEvalEngine engine) {
        super(engine);

        states.addAll(ALL);
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if ( eval == null || ! eval.isIndel() )
            return ALL;

        //NOTE: GATK3 ReferenceContext always used a window of 50 total bases, no matter the variant length.
        // This seems a little arbitrary to only look at the start + 50 (as opposed to also backward), but behavior was established in GATK3
        final byte[] bases = context.getReferenceContext().getBases(new SimpleInterval(eval.getContig(), eval.getStart(), eval.getStart() + 50));
        if ( GATKVariantContextUtils.isTandemRepeat(eval, bases) ) {
            return REPEAT;
        } else {
            return NOT_REPEAT;
        }
    }
}
