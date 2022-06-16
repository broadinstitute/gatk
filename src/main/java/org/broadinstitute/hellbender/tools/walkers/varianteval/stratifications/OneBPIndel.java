package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.Arrays;
import java.util.List;

/**
 * Stratifies the eval RODs into sites where the indel is 1 bp in length and those where the event is 2+.
 * all non indel events go into all bins, so that SNP counts can be used as contrasts in eval modules.
 */
public class OneBPIndel extends VariantStratifier {
    private final static List<Object> ALL = Arrays.asList((Object)"all", (Object)"one.bp", (Object)"two.plus.bp");
    private final static List<Object> ONE_BP = Arrays.asList((Object)"all", (Object)"one.bp");
    private final static List<Object> TWO_PLUS_BP = Arrays.asList((Object)"all", (Object)"two.plus.bp");

    public OneBPIndel(VariantEvalEngine engine) {
        super(engine);

        states.addAll(ALL);
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null && eval.isIndel()) {
            for ( int l : eval.getIndelLengths() )
                if ( Math.abs(l) > 1 )
                    return TWO_PLUS_BP; // someone is too long
            return ONE_BP; // all lengths are one
        } else
            return ALL;
    }
}
