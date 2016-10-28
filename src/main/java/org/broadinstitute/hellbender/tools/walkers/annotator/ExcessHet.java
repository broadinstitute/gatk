
package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.GenotypeCounts;
import org.broadinstitute.hellbender.utils.GenotypeUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Phred-scaled p-value for exact test of excess heterozygosity.
 * Using implementation from
 * Wigginton JE, Cutler DJ, Abecasis GR. A Note on Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 2005;76(5):887-893.
 */
public final class ExcessHet extends InfoFieldAnnotation implements StandardAnnotation {
    private static final double MIN_NEEDED_VALUE = 1.0E-16;
    public static final double PHRED_SCALED_MIN_P_VALUE = -10.0 * Math.log10(MIN_NEEDED_VALUE);

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || !vc.isVariant()) {
            return Collections.emptyMap();
        }
        final Pair<Integer, Double> sampleCountEH = calculateEH(vc, genotypes);
        final int sampleCount = sampleCountEH.getLeft();
        final double eh =  sampleCountEH.getRight();

        if (sampleCount < 1) {
            return Collections.emptyMap();
        }
        return Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", eh));
    }

    @VisibleForTesting
    Pair<Integer, Double> calculateEH(final VariantContext vc, final GenotypesContext genotypes) {
        final GenotypeCounts t = GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes);

        final int refCount = t.getRefs();
        final int hetCount = t.getHets();
        final int homCount = t.getHoms();
        // number of samples that have likelihoods
        final int sampleCount = (int) genotypes.stream().filter(g->GenotypeUtils.isDiploidWithLikelihoods(g)).count();

        final double pval = exactTest(hetCount, refCount, homCount);

        //If the actual phredPval would be infinity we will probably still filter out just a very large number
        //Since the method does not guarantee precision for any p-value smaller than 1e-16, we can return the phred scaled version
        if (pval < 10e-60) {
            return Pair.of(sampleCount, PHRED_SCALED_MIN_P_VALUE);
        }
        final double phredPval = -10.0 * Math.log10(pval);

        return Pair.of(sampleCount, phredPval);
    }

    /**
     * Note that this method is not accurate for very small p-values. Beyond 1.0E-16 there is no guarantee that the
     * p-value is accurate, just that it is in fact smaller than 1.0E-16 (and therefore we should filter it). It would
     * be more computationally expensive to calculate accuracy beyond a given threshold. Here we have enough accuracy
     * to filter anything below a p-value of 10E-6.
     *
     * @param hetCount Number of observed hets (n_ab)
     * @param refCount Number of observed homRefs (n_aa)
     * @param homCount Number of observed homVars (n_bb)
     * @return Right sided p-value or the probability of getting the observed or higher number of hets given the sample
     * size (N) and the observed number of allele a (rareCopies)
     */
    @VisibleForTesting
    static double exactTest(final int hetCount, final int refCount, final int homCount) {
        Utils.validateArg(hetCount >= 0, "Het count cannot be less than 0");
        Utils.validateArg(refCount >= 0, "Ref count cannot be less than 0");
        Utils.validateArg(homCount >= 0, "Hom count cannot be less than 0");

        //Split into observed common allele and rare allele
        final int obsHomR;
        final int obsHomC;
        if (refCount < homCount) {
            obsHomR = refCount;
            obsHomC = homCount;
        } else {
            obsHomR = homCount;
            obsHomC = refCount;
        }

        final int rareCopies = 2 * obsHomR + hetCount;
        final int N = hetCount + obsHomC + obsHomR;

        //If the probability distribution has only 1 point, then the mid p-value is .5
        if (rareCopies <= 1) {
            return .5;
        }

        final double[] probs = new double[rareCopies + 1];

        //Find (something close to the) mode for the midpoint
        int mid = (int) Math.floor(rareCopies * (2.0 * N - rareCopies) / (2.0 * N - 1.0));
        if ((mid % 2) != (rareCopies % 2)) {
            mid++;
        }

        probs[mid] = 1.0;
        double mysum = 1.0;

        //Calculate probabilities from midpoint down
        int currHets = mid;
        int currHomR = (rareCopies - mid) / 2;
        int currHomC = N - currHets - currHomR;

        while (currHets >= 2) {
            final double potentialProb = probs[currHets] * currHets * (currHets - 1.0) / (4.0 * (currHomR + 1.0) * (currHomC + 1.0));
            if (potentialProb < MIN_NEEDED_VALUE) {
                break;
            }

            probs[currHets - 2] = potentialProb;
            mysum = mysum + probs[currHets - 2];

            //2 fewer hets means one additional homR and homC each
            currHets = currHets - 2;
            currHomR = currHomR + 1;
            currHomC = currHomC + 1;
        }

        //Calculate probabilities from midpoint up
        currHets = mid;
        currHomR = (rareCopies - mid) / 2;
        currHomC = N - currHets - currHomR;

        while (currHets <= rareCopies - 2) {
            final double potentialProb = probs[currHets] * 4.0 * currHomR * currHomC / ((currHets + 2.0) * (currHets + 1.0));
            if (potentialProb < MIN_NEEDED_VALUE) {
                break;
            }

            probs[currHets + 2] = potentialProb;
            mysum = mysum + probs[currHets + 2];

            //2 more hets means 1 fewer homR and homC each
            currHets = currHets + 2;
            currHomR = currHomR - 1;
            currHomC = currHomC - 1;
        }

        final double rightPval = probs[hetCount] / (2.0 * mysum);
        //Check if we observed the highest possible number of hets
        if (hetCount == rareCopies) {
            return rightPval;
        }
        return rightPval + StatUtils.sum(Arrays.copyOfRange(probs, hetCount + 1, probs.length)) / mysum;
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.EXCESS_HET_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}
