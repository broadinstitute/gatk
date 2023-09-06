package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.*;

/**
 * Stratifies by the FILTER type(s) for each line, with PASS used for passing
 */
public class FilterType extends VariantStratifier {
    public FilterType(VariantEvalEngine engine) {
        super(engine);

        final Set<String> filterNames = new HashSet<>();
        for (FeatureInput<VariantContext> eval : getEngine().getVariantEvalArgs().getEvals()) {
            final VCFHeader header = (VCFHeader)getEngine().getHeaderForFeatures(eval);
            for (VCFFilterHeaderLine line : header.getFilterLines()) {
                filterNames.add(line.getID());
            }
        }

        states.addAll(filterNames);
        states.add("PASS");
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        if (eval == null){
            return Collections.emptyList();
        }

        final ArrayList<Object> relevantStates = new ArrayList<>();
        if (eval.isFiltered()){
            relevantStates.addAll(eval.getFilters());
        }
        else {
            relevantStates.add("PASS");
        }

        return relevantStates;
    }
}
