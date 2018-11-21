package org.broadinstitute.hellbender.tools.walkers.varianteval.stratifications;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

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
    @Override
    public void initialize() {
        states.add("all");
        states.add("CpG");
        states.add("non_CpG");
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        boolean isCpG = false;
        if (referenceContext != null && referenceContext.getBases() != null) {
            String fwRefBases = new String(referenceContext.getBases(0, 1));

            //NOTE: this matches the GATK3 behavior; however, it could be argued we should consider additional cases as valid CpGs?
            //if (leftFlank.equalsIgnoreCase("CG") || leftFlank.equalsIgnoreCase("GC") || rightFlank.equalsIgnoreCase("CG") || rightFlank.equalsIgnoreCase("GC")) {
            if (fwRefBases.startsWith("CG")) {
                isCpG = true;
            }
        }

        ArrayList<Object> relevantStates = new ArrayList<Object>(2);
        relevantStates.add("all");
        relevantStates.add(isCpG ? "CpG" : "non_CpG");

        return relevantStates;
    }
}
