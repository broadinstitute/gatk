package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEvalEngine;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.VariantEvalContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;

/**
 * CpG is a stratification module for VariantEval that divides the input data by within/not within a CpG site
 *
 * <p>
 * It is a three-state stratification:
 * <ul>
 *     <li>The locus is a CpG site ("CpG")
 *     <li>The locus is not a CpG site ("non_CpG")
 *     <li>The locus is either a CpG or not a CpG site ("all")
 * </ul>
 * A CpG site is defined as a site where the reference base at a locus is a C and the adjacent reference base in the 3' direction is a G.
 */
public class CpG extends VariantStratifier {
    public CpG(VariantEvalEngine engine) {
        super(engine);

        states.add("all");
        states.add("CpG");
        states.add("non_CpG");
    }

    @Override
    public List<Object> getRelevantStates(final VariantEvalContext context, final VariantContext comp, final String compName, final VariantContext eval, final String evalName, final String sampleName, final String familyName) {
        boolean isCpG = false;
        if (context.getReferenceContext() != null && context.getReferenceContext().getBases() != null) {
            String fwRefBases = new String(context.getReferenceContext().getBases(new SimpleInterval(context.getReferenceContext().getContig(), context.getReferenceContext().getStart(), context.getReferenceContext().getStart() + 1)));

            //NOTE: this matches the GATK3 behavior; however, it could be argued we should consider additional cases as valid CpGs?
            //if (leftFlank.equalsIgnoreCase("CG") || leftFlank.equalsIgnoreCase("GC") || rightFlank.equalsIgnoreCase("CG") || rightFlank.equalsIgnoreCase("GC")) {
            if (fwRefBases.startsWith("CG")) {
                isCpG = true;
            }
        }

        final ArrayList<Object> relevantStates = new ArrayList<>(2);
        relevantStates.add("all");
        relevantStates.add(isCpG ? "CpG" : "non_CpG");

        return relevantStates;
    }
}
