package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;

/**
 * Stratifies by the FILTER type(s) for each line, with PASS used for passing
 */
public class FilterType extends VariantStratifier {
    @Override
    public void initialize() {
        Set<String> filterNames = new HashSet<>();
        for (FeatureInput<VariantContext> eval : getVariantEvalWalker().getEvals()) {
            VCFHeader header = (VCFHeader)getVariantEvalWalker().getHeaderForFeatures(eval);
            for (VCFFilterHeaderLine line : header.getFilterLines()) {
                filterNames.add(line.getID());
            }
        }

        states.addAll(filterNames);
        states.add("PASS");
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String familyName) {
        if (eval == null){
            return Collections.emptyList();
        }

        ArrayList<Object> relevantStates = new ArrayList<Object>();
        if (eval.isFiltered()){
            relevantStates.addAll(eval.getFilters());
        }
        else {
            relevantStates.add("PASS");
        }

        return relevantStates;
    }
}
