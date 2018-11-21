package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

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

    @Override
    public void initialize() {
        for( int a=-MAX_INDEL_SIZE; a <=MAX_INDEL_SIZE; a++ ) {
            states.add(a);
        }
    }

    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
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

                return Collections.singletonList((Object)eventLength);
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
