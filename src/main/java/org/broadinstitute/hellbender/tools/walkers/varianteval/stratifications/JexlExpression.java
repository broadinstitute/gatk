package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Stratifies the eval RODs by user-supplied JEXL expressions
 *
 * https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions for more details
 */
public class JexlExpression extends VariantStratifier implements StandardStratification {
    // needs to know the jexl expressions
    private Set<SortableJexlVCMatchExp> jexlExpressions;

    public JexlExpression(VariantEvalEngine engine) {
        super(engine);

        jexlExpressions = getEngine().getJexlExpressions();

        states.add("none");
        for ( SortableJexlVCMatchExp jexlExpression : jexlExpressions ) {
            states.add(jexlExpression.name);
        }
    }

    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        final ArrayList<Object> relevantStates = new ArrayList<>();
        relevantStates.add("none");

        for ( SortableJexlVCMatchExp jexlExpression : jexlExpressions ) {
            if (eval != null && VariantContextUtils.match(eval, jexlExpression)) {
                relevantStates.add(jexlExpression.name);
            }
        }

        return relevantStates;
    }
}
