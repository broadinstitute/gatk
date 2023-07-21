package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.SnpEffUtil;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies by nonsense, missense, silent, and all annotations in the input ROD, from the INFO field annotation.
 */
public class FunctionalClass extends VariantStratifier {

    private static final String FUNCTIONAL_CLASS_KEY = "SNPEFF_FUNCTIONAL_CLASS";

    public enum FunctionalType {
        silent,
        missense,
        nonsense
    }

    public FunctionalClass(VariantEvalEngine engine) {
        super(engine);

        states.add("all");
        for ( FunctionalType type : FunctionalType.values() )
            states.add(type.name());
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        final ArrayList<Object> relevantStates = new ArrayList<>();

        relevantStates.add("all");

        if (eval != null && eval.isVariant()) {
            FunctionalType type = null;

            if (eval.hasAttribute("refseq.functionalClass")) {
                try {
                    type = FunctionalType.valueOf(eval.getAttributeAsString("refseq.functionalClass", null));
                } catch ( Exception e ) {} // don't error out if the type isn't supported
            } else if (eval.hasAttribute("refseq.functionalClass_1")) {
                int annotationId = 1;
                String key;

                do {
                    key = String.format("refseq.functionalClass_%d", annotationId);

                    String newtypeStr = eval.getAttributeAsString(key, null);
                    if ( newtypeStr != null && !newtypeStr.equalsIgnoreCase("null") ) {
                        try {
                            final FunctionalType newType = FunctionalType.valueOf(newtypeStr);
                            if ( type == null ||
                                    ( type == FunctionalType.silent && newType != FunctionalType.silent ) ||
                                    ( type == FunctionalType.missense && newType == FunctionalType.nonsense ) ) {
                                type = newType;
                            }
                        } catch ( Exception e ) {} // don't error out if the type isn't supported
                    }

                    annotationId++;
                } while (eval.hasAttribute(key));

            } else if ( eval.hasAttribute(FUNCTIONAL_CLASS_KEY) ) {
                try {
                    final SnpEffUtil.EffectFunctionalClass snpEffFunctionalClass = SnpEffUtil.EffectFunctionalClass.valueOf(eval.getAttribute(FUNCTIONAL_CLASS_KEY).toString());
                    if ( snpEffFunctionalClass == SnpEffUtil.EffectFunctionalClass.NONSENSE )
                        type = FunctionalType.nonsense;
                    else if ( snpEffFunctionalClass == SnpEffUtil.EffectFunctionalClass.MISSENSE )
                        type = FunctionalType.missense;
                    else if ( snpEffFunctionalClass == SnpEffUtil.EffectFunctionalClass.SILENT )
                        type = FunctionalType.silent;
                }
                catch ( Exception e ) {} // don't error out if the type isn't supported
            }

            if ( type != null ) {
                relevantStates.add(type.name());
            }
        }

        return relevantStates;
    }
}
