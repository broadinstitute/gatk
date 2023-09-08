package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.Collections;
import java.util.List;

/**
 * Stratifies the eval RODs by the indel size
 *
 * Indel sizes are stratified from sizes -100 to +100. Sizes greater than this are lumped in the +/- 100 bin
 * This stratification ignores multi-allelic indels (whose size is not defined uniquely)
 */
public class IndelSize extends VariantStratifier {
    static final int MAX_INDEL_SIZE = 100;

    public IndelSize(VariantEvalEngine engine) {
        super(engine);

        for ( int a=-MAX_INDEL_SIZE; a <=MAX_INDEL_SIZE; a++ ) {
            states.add(a);
        }
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval != null && eval.isIndel() && eval.isBiallelic()) {
            try {
                int eventLength = 0;
                if ( eval.isSimpleInsertion() ) {
                    eventLength = eval.getAlternateAllele(0).length();
                } else if ( eval.isSimpleDeletion() ) {
                    eventLength = -eval.getReference().length();
                }

                if (eventLength > MAX_INDEL_SIZE)
                    eventLength = MAX_INDEL_SIZE;
                else if (eventLength < -MAX_INDEL_SIZE)
                    eventLength = -MAX_INDEL_SIZE;

                return Collections.singletonList(eventLength);
            } catch (Exception e) {
                return Collections.emptyList();
            }
        }

        return Collections.emptyList();
    }

    @Override
    public String getFormat() {
        return "%d";
    }
}
