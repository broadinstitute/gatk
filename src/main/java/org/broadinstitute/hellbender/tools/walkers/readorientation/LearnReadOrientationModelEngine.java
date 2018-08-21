package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Histogram;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.MathArrays;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;


public class LearnReadOrientationModelEngine {
    // We consider the likelihood converged to its maximum when the difference falls below this threshold
    private final double convergenceThreshold;

    // If the EM does not converge in a few steps we should suspect that something went wrong
    private final int maxEMIterations;

    private final String referenceContext;

    private final Nucleotide refAllele;

    private final Histogram<Integer> refHistogram;

    private final List<Histogram<Integer>> altDepthOneHistograms;

    private final List<AltSiteRecord> altDesignMatrix;

    /**
     * N by K matrix of posterior probabilities of latent variable z, where N is the number of alt sites,
     * evaluated at the current estimates of the mixture weights
     */
    private final RealMatrix altResponsibilities;

    private final Map<Triple<Integer, Nucleotide, ReadOrientation>, double[]> responsibilitiesOfAltDepth1Sites;

    /**
     * MAX_COVERAGE by K matrix of responsibilities of a ref site (i.e. ALT Depth = 0, ALT F1R2 = 0)
     * for ref sites with coverage 1, 2,..., maxDepth
     */
    private final RealMatrix refResponsibilities;

    private final int numAltExamples;

    private final int numRefExamples;

    private final int numExamples;

    // K-dimensional vector of effective sample counts for each class of z, weighted by the the altResponsibilities. For a fixed k,
    // we sum up the counts over all alleles. N_k in the docs.
    private RealVector effectiveCounts = new ArrayRealVector(F1R2FilterConstants.NUM_STATES);

    // K-dimensional vector of beta-binomial parameters alpha and beta for
    // the conditional distribution over the alt count m given artifact state z
    private final static Map<ArtifactState, BetaDistributionShape> alleleFractionPseudoCounts = getPseudoCountsForAlleleFraction();

    // K-dimensional vector of beta-binomial parameters alpha and beta for
    // the conditional distribution over the alt F1R2 count x given alt count m
    private final static Map<ArtifactState, BetaDistributionShape> altF1R2FractionPseudoCounts = getPseudoCountsForAltF1R2Fraction();

    private final MutableInt numIterations = new MutableInt();

    private final Logger logger;

    private int maxDepth;

    /** Fixed hyperparameters for betabinomials ***/
    //
    private final static double ALT_PSEUDOCOUNT = 1.0;

    private final static double REF_PSEUDOCOUNT = 9.0;

    // For home ref sites assume Q35
    private final static double PSEUDOCOUNT_OF_HOM_LIKELY = 10000.0;

    private final static double PSEUDOCOUNT_OF_HOM_UNLIKELY = 3.0;

    // To account for potential copy number variation we use a relatively broad allele fraction distribution.
    // In the future we will use the segmentation info as done in the germline filter and the contamination model
    private final static double BALANCED_HET_PSEUDOCOUNT = 5.0;

    // Give a sharper beta prior than the allele fraction
    private final static double BALANCED_F1R2_PRIOR = 10;

    // These variables define the distribution of allele fraction in a somatic variant. It should be learned from
    // the data in the future. In the meantime, one should tweak this by hand when e.g. applying the read orientation
    // filter on low allele fraction samples such as the blood biopsy
    private final static double PSEUDOCOUNT_OF_SOMATIC_ALT = 2.0;

    private final static double PSEUDOCOUNT_OF_SOMATIC_REF = 5.0;

    private final static double PSEUDOCOUNT_OF_LIKELY_OUTCOME = 100.0;

    private final static double PSEUDOCOUNT_OF_RARE_OUTCOME = 1.0;
    /*** END hyperparameters ***/

    /**
     * Contract: the reference contexts must be combined with its reverse complements prior to instantiating this class
     */
    public LearnReadOrientationModelEngine(final Histogram<Integer> refHistogram, final List<Histogram<Integer>> altDepthOneHistograms,
                                           final List<AltSiteRecord> altDesignMatrixForContext,
                                           final double convergenceThreshold, final int maxEMIterations,
                                           final int maxDepth, final Logger logger) {
        this.refHistogram = Utils.nonNull(refHistogram);
        this.altDepthOneHistograms = Utils.nonNull(altDepthOneHistograms);
        this.altDesignMatrix = Utils.nonNull(altDesignMatrixForContext);
        this.referenceContext = refHistogram.getValueLabel();
        Utils.validate(referenceContext.length() == F1R2FilterConstants.REFERENCE_CONTEXT_SIZE,
                String.format("reference context must have length %d but got %s", F1R2FilterConstants.REFERENCE_CONTEXT_SIZE, referenceContext));
        Utils.validate(F1R2FilterConstants.CANONICAL_KMERS.contains(referenceContext),
                referenceContext + " is not in the set of canonical kmers");
        this.numAltExamples = altDesignMatrix.size() + altDepthOneHistograms.stream().mapToInt(h -> (int) h.getSumOfValues()).sum();
        this.numRefExamples = (int) refHistogram.getSumOfValues();
        this.numExamples = numAltExamples + numRefExamples;

        // Responsibilities of ref sites with equal depth are the same so we can compute it for each depth and
        // multiply by the number of counts for that depth
        this.refResponsibilities = new Array2DRowRealMatrix(maxDepth, F1R2FilterConstants.NUM_STATES);

        this.altResponsibilities = new Array2DRowRealMatrix(altDesignMatrix.size(), F1R2FilterConstants.NUM_STATES);

        // Store responsibilities for each depth and the F1R2/F2R1 of the one alt read
        this.responsibilitiesOfAltDepth1Sites = new HashMap<>();
        this.refAllele = F1R2FilterUtils.getMiddleBase(referenceContext);
        this.convergenceThreshold = convergenceThreshold;
        this.maxEMIterations = maxEMIterations;
        this.maxDepth = maxDepth;
        this.logger = logger;
    }

    // Learn the prior probabilities for the artifact states by the EM algorithm
    public ArtifactPrior learnPriorForArtifactStates() {
        // Initialize the prior for artifact
        double[] statePrior = getFlatPrior(refAllele);
        double l2Distance;

        do {
            final double[] oldStatePrior = Arrays.copyOf(statePrior, F1R2FilterConstants.NUM_STATES);

            // Responsibilities are updated by side effect to save space
            takeEstep(statePrior);
            statePrior = takeMstep();

            // TODO: make sure EM increases the likelihood
            // newLikelihood >= oldLikelihood : "M step must increase the likelihood";
            l2Distance = MathArrays.distance(oldStatePrior, statePrior);

            numIterations.increment();
        } while (l2Distance > convergenceThreshold && numIterations.intValue() < maxEMIterations);

        if (numIterations.intValue() == maxEMIterations){
            logger.info(String.format("Context %s: with %s ref and %s alt examples, EM failed to converge within %d steps",
                    referenceContext, numRefExamples, numAltExamples, maxEMIterations));
        } else {
            logger.info(String.format("Context %s: with %s ref and %s alt examples, EM converged in %d steps",
                    referenceContext, numRefExamples, numAltExamples, numIterations.intValue()));
        }

        return new ArtifactPrior(referenceContext, statePrior, numExamples, numAltExamples);
    }

    /**
     * Given the current estimates of artifact prior probabilities, compute the responsibilities, which are
     * the posterior probabilities of artifact states, for each data point
     **/
    private void takeEstep(final double[] artifactPriors) {
        /**
         * Compute the responsibilities of ref examples.
         * Given that for moderate to high depths we will always have P(HOM REF) = 1, this is largely overkill
         *
         * Ref sites with the same depth have the same alt and alt F1R2 depths (i.e. zero) so avoid repeated computations
         */
        for (int i = 0; i < maxDepth; i++) {
            final int depth = i + 1;
            refResponsibilities.setRow(i, computeResponsibilities(refAllele, refAllele, 0, 0, depth, artifactPriors, false));
        }

        // Compute the responsibilities of alt sites
        for (int n = 0; n < altDesignMatrix.size(); n++) {
            final AltSiteRecord example = altDesignMatrix.get(n);
            final int depth = example.getDepth();
            final int altDepth = example.getAltCount();
            final int altF1R2 = example.getAltF1R2();

            altResponsibilities.setRow(n, computeResponsibilities(refAllele, example.getAltAllele(), altDepth, altF1R2, depth, artifactPriors, false));
        }

        // Compute the responsibilities of alt sites with depth=1
        for (int i = 0; i < maxDepth; i++){
            final int depth = i+1;
            for (Nucleotide altAllele : Nucleotide.STANDARD_BASES){
                for (ReadOrientation orientation : ReadOrientation.values()){
                    if (altAllele == refAllele){
                        continue;
                    }

                    final int f1r2Depth = orientation == ReadOrientation.F1R2 ? 1 : 0;

                    final Triple<Integer, Nucleotide, ReadOrientation> key = createKey(depth, altAllele, orientation);
                    responsibilitiesOfAltDepth1Sites.put(key, computeResponsibilities(refAllele, altAllele, 1, f1r2Depth, depth, artifactPriors, false));
                }
            }
        }
    }

    /**
     * Given the posterior distributions over the artifact states (ie responsibilities) under the current estimates for the prior,
     * update the prior weights such that they maximize the lower bound on the marginal likelihood P(Data).
     * We do so by maximizing the expectation of the complete data likelihood with respect to the posterior
     * for the artifact states from the E-step
     */
    private double[] takeMstep() {
        // First we compute the effective counts of each state, N_k in the docs. We do this separately over alt and ref sites
        final double[] effectiveAltCountsFromDesignMatrix = MathUtils.sumArrayFunction(0, altDesignMatrix.size(), n -> altResponsibilities.getRow(n));
        double[] effectiveAltCountsFromHistograms = new double[F1R2FilterConstants.NUM_STATES];

        for (Histogram<Integer> histogram : altDepthOneHistograms){
            final Triple<String, Nucleotide, ReadOrientation> triplet = F1R2FilterUtils.labelToTriplet(histogram.getValueLabel());
            final Nucleotide altAllele = triplet.getMiddle();
            final ReadOrientation orientation = triplet.getRight();


            final double[] effectiveAltCountsFromHistogram = MathUtils.sumArrayFunction(0, maxDepth,
                    i -> MathArrays.scale(histogram.get(i + 1).getValue(), responsibilitiesOfAltDepth1Sites.get(createKey(i+1, altAllele, orientation))));
            effectiveAltCountsFromHistograms = MathArrays.ebeAdd(effectiveAltCountsFromHistograms, effectiveAltCountsFromHistogram);
        }

        final double[] effectiveAltCounts = MathArrays.ebeAdd(effectiveAltCountsFromDesignMatrix, effectiveAltCountsFromHistograms);

        // TODO: at some depth, the responsibilities must be 1 for z = HOM_REF and 0 for everything else, we could probably save some time there
        // Over ref sites, we have a histogram of sites over different depths. At each depth we simply multiply the responsibilities by the number of sites,
        // and sum them over all of depths. Because we cut off the depth histogram at {@code MAX_COVERAGE}, we underestimate the ref effective counts by design
        final double[] effectiveRefCounts = MathUtils.sumArrayFunction(0, maxDepth,
                i -> MathArrays.scale(refHistogram.get(i + 1).getValue(), refResponsibilities.getRow(i)));

        effectiveCounts = new ArrayRealVector(MathArrays.ebeAdd(effectiveAltCounts, effectiveRefCounts));
        return effectiveCounts.mapMultiply(1.0/numExamples).toArray();
    }

    /**
     * Return normalized probabilities
     */
    public static double[] computeResponsibilities(final Nucleotide refAllele, final Nucleotide altAllele,
                                                   final int altDepth, final int f1r2AltCount, final int depth,
                                                   final double[] artifactPrior, final boolean givenNotHomRef) {
        final double[] log10UnnormalizedResponsibilities = new double[F1R2FilterConstants.NUM_STATES];
        final List<ArtifactState> refToRefArtifacts = ArtifactState.getRefToRefArtifacts(refAllele);

        for (ArtifactState state : ArtifactState.values()){
            final int stateIndex = state.ordinal();
            if (refToRefArtifacts.contains(state)) {
                // This state is really just hom ref so give it zero probability and skip
                log10UnnormalizedResponsibilities[stateIndex] = Double.NEGATIVE_INFINITY;
                continue;
            }

            if (ArtifactState.artifactStates.contains(state) && state.getAltAlleleOfArtifact() != altAllele) {
                // The indicator function is 0
                log10UnnormalizedResponsibilities[stateIndex] = Double.NEGATIVE_INFINITY;
                continue;
            }

            // If we get here, we have a non-artifact state i.e. { germline het, hom ref, hom var, somatic het }
            // or an artifact state whose transitions match the observed alt allele (e.g. alt allele = A, z = F1R2_A, F2R1_A)
            log10UnnormalizedResponsibilities[stateIndex] = computePosterior(altDepth, f1r2AltCount, depth, artifactPrior[stateIndex],
                    alleleFractionPseudoCounts.get(state), altF1R2FractionPseudoCounts.get(state));
        }

        if (givenNotHomRef){
            log10UnnormalizedResponsibilities[ArtifactState.HOM_REF.ordinal()] = Double.NEGATIVE_INFINITY;
        }

        return MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);
    }


    /**
     * Compute the posterior probability of the state z given data. The caller is responsible for not calling
     * this method on inconsistent states e.g. z = F1R2_C where the reference context is ACT
     */
    private static double computePosterior(final int altDepth, final int altF1R2Depth, final int depth,
                                           final double statePrior, final BetaDistributionShape afPseudoCounts,
                                           final BetaDistributionShape f1r2PseudoCounts){
        Utils.validateArg(MathUtils.isAProbability(statePrior), String.format("statePrior must be a probability but got %f", statePrior));

        return Math.log10(statePrior) +
                MathUtils.log10BetaBinomialProbability(altDepth, depth, afPseudoCounts.getAlpha(), afPseudoCounts.getBeta()) +
                MathUtils.log10BetaBinomialProbability(altF1R2Depth, altDepth, f1r2PseudoCounts.getAlpha(), f1r2PseudoCounts.getBeta());
    }

    /**
     * For each state z, define allele fraction distribution f_z and alt f1r2 fraction theta_z
     * They are both beta binomial distributions and are therefore parameterized by the pseudocounts alpha and beta
     */
    private static Map<ArtifactState, BetaDistributionShape> getPseudoCountsForAlleleFraction(){
        final Map<ArtifactState, BetaDistributionShape> alleleFractionPseudoCounts = new HashMap<>(ArtifactState.values().length);

        // The allele fraction distribution, which is not aware of the read orientation, should be the same between
        // F1R2 and F2R1 artifacts
        ArtifactState.getF1R2ArtifactStates().forEach(s -> alleleFractionPseudoCounts.put(s, new BetaDistributionShape(ALT_PSEUDOCOUNT, REF_PSEUDOCOUNT)));
        ArtifactState.getF2R1ArtifactStates().forEach(s -> alleleFractionPseudoCounts.put(s, new BetaDistributionShape(ALT_PSEUDOCOUNT, REF_PSEUDOCOUNT)));

        alleleFractionPseudoCounts.put(ArtifactState.HOM_REF, new BetaDistributionShape(PSEUDOCOUNT_OF_HOM_UNLIKELY, PSEUDOCOUNT_OF_HOM_LIKELY));
        alleleFractionPseudoCounts.put(ArtifactState.GERMLINE_HET, new BetaDistributionShape(BALANCED_HET_PSEUDOCOUNT, BALANCED_HET_PSEUDOCOUNT));
        alleleFractionPseudoCounts.put(ArtifactState.SOMATIC_HET, new BetaDistributionShape(PSEUDOCOUNT_OF_SOMATIC_ALT, PSEUDOCOUNT_OF_SOMATIC_REF));
        alleleFractionPseudoCounts.put(ArtifactState.HOM_VAR, new BetaDistributionShape(PSEUDOCOUNT_OF_HOM_LIKELY, PSEUDOCOUNT_OF_HOM_UNLIKELY));

        return alleleFractionPseudoCounts;
    }

    private static Map<ArtifactState, BetaDistributionShape> getPseudoCountsForAltF1R2Fraction(){
        final Map<ArtifactState, BetaDistributionShape> altF1R2FractionPseudoCounts = new HashMap<>(ArtifactState.values().length);

        ArtifactState.getF1R2ArtifactStates().forEach(z -> altF1R2FractionPseudoCounts.put(z, new BetaDistributionShape(PSEUDOCOUNT_OF_LIKELY_OUTCOME, PSEUDOCOUNT_OF_RARE_OUTCOME)));
        ArtifactState.getF2R1ArtifactStates().forEach(z -> altF1R2FractionPseudoCounts.put(z, new BetaDistributionShape(PSEUDOCOUNT_OF_RARE_OUTCOME, PSEUDOCOUNT_OF_LIKELY_OUTCOME)));
        ArtifactState.getNonArtifactStates().forEach(z -> altF1R2FractionPseudoCounts.put(z, new BetaDistributionShape(BALANCED_F1R2_PRIOR, BALANCED_F1R2_PRIOR)));

        return altF1R2FractionPseudoCounts;
    }

    @VisibleForTesting
    public double[] getRefResonsibilities(final int rowNum){
        return refResponsibilities.getRow(rowNum);
    }

    @VisibleForTesting
    public double[] getAltResonsibilities(final int rowNum){
        return altResponsibilities.getRow(rowNum);
    }

    @VisibleForTesting
    public double[] getAltDepth1Resonsibilities(final int rowNum){
        return null;
    }

    @VisibleForTesting
    public RealVector getEffectiveCounts(){
        return effectiveCounts;
    }

    @VisibleForTesting
    public double getEffectiveCounts(ArtifactState state){
        return effectiveCounts.getEntry(state.ordinal());
    }

    public static double[] getFlatPrior(final Nucleotide refAllele){
        // We skip the artifact states in which ref transitions to ref e.g. under the ref context AGT, F1R2_G and F2R1_G
        final List<ArtifactState> refToRefStates = ArtifactState.getRefToRefArtifacts(refAllele);

        double[] prior = new double[F1R2FilterConstants.NUM_STATES];

        Arrays.fill(prior, 1.0 / (F1R2FilterConstants.NUM_STATES - refToRefStates.size()));
        for (ArtifactState s : refToRefStates){
            prior[s.ordinal()] = 0;
        }

        return prior;
    }

    private Triple<Integer, Nucleotide, ReadOrientation> createKey(final int depth, final Nucleotide altAllele, final ReadOrientation orientation){
        return new ImmutableTriple<>(depth, altAllele, orientation);
    }
}