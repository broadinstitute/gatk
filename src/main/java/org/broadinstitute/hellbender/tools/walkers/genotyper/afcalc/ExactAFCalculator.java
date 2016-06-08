package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.StreamSupport;

/**
 * Uses the Exact calculation of Heng Li
 */
abstract class ExactAFCalculator extends AFCalculator {

    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first

    /**
     * Unpack GenotypesContext into arraylist of double values
     * @param GLs            Input genotype context
     * @param includeDummy   //TODO: Does anyone have any clue what this is?????
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    public static List<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy) {
        final List<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);
        if (includeDummy) {
            genotypeLikelihoods.add((new double[]{0.0, 0.0, 0.0}));
        }

        StreamSupport.stream(GLs.iterateInSampleNameOrder().spliterator(), false)
                .filter(Genotype::hasLikelihoods)
                .map(gt -> gt.getLikelihoods().getAsVector())
                .filter(GATKVariantContextUtils::isInformative)
                .forEach(genotypeLikelihoods::add);

        return genotypeLikelihoods;
    }
}