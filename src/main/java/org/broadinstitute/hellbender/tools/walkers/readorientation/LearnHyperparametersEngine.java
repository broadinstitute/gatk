package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Histogram;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.LongStream;

import static org.broadinstitute.hellbender.tools.walkers.readorientation.ArtifactType.F1R2;
import static org.broadinstitute.hellbender.tools.walkers.readorientation.ReadOrientationFilterConstants.*;


/**
 * Created by tsato on 7/26/17.
 */
public class LearnHyperparametersEngine {
    public static final int NUM_STATES = ArtifactState.values().length;

    // When the increase in likelihood falls under this value, we call the algorithm converged
    private final double convergenceThreshold;

    // If the EM does not converge in a few steps we should suspect that something went wrong
    private final int maxEMIterations;

    // Regularizer (?) TODO: think this through
    static final double EPSILON = 1e-3;


    final String referenceContext;

    final Nucleotide refAllele;

    final Histogram<Integer> refHistogram;

    final List<Histogram<Integer>> altHistograms;

    final List<AltSiteRecord> altDesignMatrix;

    // N by K matrix of posterior probabilities of latent variable z, where N is the number of alt sites,
    // evaluated at the current estimates of the mixture weights pi
    final double[][] altResponsibilities;


    // {@code MAX_COVERAGE + 1} by 2 by K matrix of responsibilities of a alt site with alt count = 1
    // [DEPTH][STATES][2], where 2 comes from {F1R2, F2R1}
    final double[][][] responsibilitiesOfAltDepth1Sites;

    // {@code MAX_COVERAGE + 1} by K matrix of a cache of responsibilities of a ref site (i.e. m = 0, x = 0)
    // for ref sites with coverage 1, 2,..., maxDepthForHistograms. To reiterate, the rows represent different coverages,
    // not samples (the count of samples in each coverage is stored in a separate histogram) and index 0 corresponds to
    // coverage 1; coverage = index + 1
    final double[][] refResponsibilities;

    final int numAltExamples;

    final int numRefExamples;

    final int numExamples;

    final List<ArtifactState> impossibleStates;

    // K-dimensional vector of effective sample counts for each class of z, weighted by the the altResponsibilities. For a fixed k,
    // we sum up the counts over all alleles. N_k in the docs.
    // TODO: should be final
    @VisibleForTesting
    double[] effectiveCounts = new double[NUM_STATES];

    /*** Hyperparameters of the model ***/

    // pi is the K-dimensional vector of probabilities for the categorical variable z. Adds up to 1.0
    double[] pi = new double[NUM_STATES];

    // K-dimensional vector of beta-binomial parameters to the conditional r.v. m given z
    public static List<Pair<Double, Double>> alleleFractionPseudoCounts = getPseudoCountsForAlleleFraction();

    // K-dimensional vector of beta-binomial parameters to the conditional r.v. m given z
    public static List<Pair<Double, Double>> altF1R2FractionPseudoCounts = getPseudoCountsForAltF1R2Fraction();

    private int numIterations = 0;

    final Logger logger;

    // statistics about the alt design matrix - use them to validate e.g. effective alt depth over all samples
    final long sumRefDepth;
    final int sumAltDepth;
    final int sumF1R2AltDepth;

    public LearnHyperparametersEngine(final Histogram<Integer> refHistogram, final List<Histogram<Integer>> altHistograms,
                                      final List<AltSiteRecord> altDesignMatrixForContext,
                                      final double convergenceThreshold, final int maxEMIterations, final Logger logger) {
        // TODO: count up the number of rows in altDataTable
        // TODO: restructure the code such that learn hyperparameters is a tool but LearnHyperparametersEngine is a separate class with a constructor,
        // and you initialize the object for each reference context
        Utils.nonNull(refHistogram);
        Utils.nonNull(altHistograms);
        Utils.nonNull(altDesignMatrixForContext);

        referenceContext = refHistogram.getValueLabel();
        Utils.validateArg(referenceContext.length() == REFERENCE_CONTEXT_SIZE,
                String.format("reference context must have length %d but got %s", REFERENCE_CONTEXT_SIZE, referenceContext));
        Utils.validateArg(referenceContext.matches(String.format("[ACGT]{%d}", REFERENCE_CONTEXT_SIZE)),
                String.format("reference context must consist of bases A,C,G,T but got %s", referenceContext));

        this.refHistogram = refHistogram;
        this.altHistograms = altHistograms;
        this.altDesignMatrix = altDesignMatrixForContext;
        numAltExamples = altDesignMatrix.size() + altHistograms.stream().mapToInt(h -> (int) h.getSumOfValues()).sum();
        numRefExamples = (int) refHistogram.getSumOfValues();
        numExamples = numAltExamples + numRefExamples;
        refResponsibilities = new double[maxDepthForHistograms][NUM_STATES];
        altResponsibilities = new double[altDesignMatrix.size()][NUM_STATES];
        responsibilitiesOfAltDepth1Sites = new double[maxDepthForHistograms][ArtifactState.values().length][NUM_STATES];
        refAllele = Nucleotide.valueOf(referenceContext.substring(MIDDLE_INDEX, MIDDLE_INDEX + 1));
        this.convergenceThreshold = convergenceThreshold;
        this.maxEMIterations = maxEMIterations;
        this.logger = logger;

        sumRefDepth = LongStream.range(1, maxDepthForHistograms).map(c -> c * (long) refHistogram.get((int) c).getValue()).sum();
        sumAltDepth = altDesignMatrixForContext.stream().mapToInt(alt -> alt.getBaseCounts()[alt.getAltAllele().ordinal()]).sum();
        sumF1R2AltDepth = altDesignMatrixForContext.stream().mapToInt(alt -> alt.getF1R2Counts()[alt.getAltAllele().ordinal()]).sum();

        // Initialize pi
        // ===========================
        // Some artifact states don't make sense for a given context and therefore should be given the probability of 0
        // e.g. under the ref context AGT, F1R2_G and F2R1_G states are impossible, because by definition a read orientation
        // artifact only applies to alt sites (REWORD). We would remove those entries from the vector but it simplifies the
        // implementation (i.e. we can share the same indices for the states across all contexts) if we just set some
        // probabilities to 0
        impossibleStates = ArtifactState.getImpossibleStates(refAllele);
        Arrays.fill(pi, 1.0 / (NUM_STATES - impossibleStates.size()));
        for (ArtifactState impossibleState : impossibleStates) {
            pi[impossibleState.ordinal()] = 0;
        }
    }

    public Hyperparameters runEMAlgorithm() {
        boolean converged = false;
        double[] oldPi = new double[NUM_STATES];
        // one may plot the changes in L2 distance of parameters to make sure that
        // EM is steadily moving towards the (local? global?) maximum
        double[] l2distancesOfParameters = new double[maxEMIterations];

        while (!converged && numIterations < maxEMIterations) {
            // TODO: stylistic problems here, there's too many side-effects
            takeEstep();
            // FIXME: perhaps M-step shoudl return a Hyperparameters data structure
            takeMstep(); // hyperparameters e.g. {@code pi} are updated via side-effect

            // assert newLikelihood >= oldLikelihood : "M step must increase the likelihood";
            final double l2Distance = MathArrays.distance(oldPi, pi);
            converged = l2Distance < convergenceThreshold;

            l2distancesOfParameters[numIterations] = l2Distance;

            oldPi = Arrays.copyOf(pi, NUM_STATES);


            numIterations++;
        }

        logger.info(String.format("Context %s: with %s ref and %s alt examples, EM converged in %d steps",
                referenceContext, numRefExamples, numAltExamples, numIterations));
        logger.info(String.format("The changes in L2 distance of pi between iterations %s",
                Hyperparameters.doubleArrayToString(l2distancesOfParameters)));
        return new Hyperparameters(referenceContext, pi, numExamples, numAltExamples);
    }

    // Given the current estimates of the parameters pi, compute the log10AltResponsibilities
    // gamma_nk = p(z_nk)
    private void takeEstep() {
        // We save some computation here by recognizing that ref sites with the same depth have the same alt depth and
        // alt F1R2 depth (i.e. m = x = 0). Thus responsibilities for ref sites are a function only of the depth (ref and
        // alt combined) and therefore we need only compute the responsibility once for unique depth 0, 1, ..., MAX_COVERAGE
        for (int depth = 1; depth <= maxDepthForHistograms; depth++) {
            final int r = depth; // another hack to use depth in a stream

            final double[] log10UnnormalizedResponsibilities = computeLog10Responsibilities(refAllele, refAllele, 0, 0, depth, pi);

            // TODO: computing responsibilities for ref will not be too expensive, but, as noted elsewhere,
            // It may make sense to peg z = hom ref at 1.0 at certain coverage
            refResponsibilities[r - 1] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);
            Utils.validate(Math.abs(MathUtils.sum(refResponsibilities[r - 1]) - 1.0) < EPSILON,
                    String.format("ref responsibility for depth = %d added up to %f", r, MathUtils.sum(refResponsibilities[r - 1])));
        }

        // Compute the altResponsibilities of each of n alt sites \gamma_{nk}
        for (int n = 0; n < altDesignMatrix.size(); n++) {
            final AltSiteRecord example = altDesignMatrix.get(n);

            final int depth = example.getDepth();

            // TODO: test that one alt allele approach. if it's successful, simplify the code (i.e. no array, use pileupsummary, etc.)
            final int[] baseCounts = example.getBaseCounts();
            final int[] f1r2Counts = example.getF1R2Counts();
            final Nucleotide altAllele = example.getAltAllele();

            final int depthOfMostLikelyAltAllele = baseCounts[altAllele.ordinal()];
            final int f1r2DepthOfMostLikelyAltAllele = f1r2Counts[altAllele.ordinal()];

            // K-dimensional array of one of the terms that comprises gamma*_{nk}
            double[] log10UnnormalizedResponsibilities = computeLog10Responsibilities(refAllele, altAllele, depthOfMostLikelyAltAllele, f1r2DepthOfMostLikelyAltAllele,
                    depth, pi);

            // we normalize responsibilities here because the M-step uses normalized responsibilities
            altResponsibilities[n] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);

            Utils.validate(Math.abs(MathUtils.sum(altResponsibilities[n]) - 1.0) < EPSILON,
                    String.format("responsibility for %dth example added up to %f", n, MathUtils.sumLog10(altResponsibilities[n])));
        }

        // Compute the responsibilities of sites of ref sites (alt depth = 0)
        // and alt sites with dept=1
        for (int d = 1; d <= maxDepthForHistograms; d++){
            // Streamify this?
            for (Nucleotide altAllele : REGULAR_BASES){
                for (ArtifactType artifactType : ArtifactType.values()){
                    if (altAllele.toString().equals(referenceContext.substring(MIDDLE_INDEX, MIDDLE_INDEX+1))){
                        continue;
                    }

                    final int[] baseCounts = new int[REGULAR_BASES.size()];
                    baseCounts[altAllele.ordinal()] = 1;

                    final int[] f1r2Counts = new int[REGULAR_BASES.size()];
                    f1r2Counts[altAllele.ordinal()] = artifactType == F1R2 ? 1 : 0;

                    final int f1r2Depth = artifactType == F1R2 ? 1 : 0;

                    double[] log10UnnormalizedResponsibilities = computeLog10Responsibilities(
                            refAllele, altAllele, 1, f1r2Depth, d, pi);

                    responsibilitiesOfAltDepth1Sites[d-1][artifactType.ordinal()] =
                            MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);
                }
            }
        }
    }

    // Given the current posterior distributions over z (aka repsonsibilities), compute the estimate of mixture weights
    // for the categorial variable z (pi) that maximize the lower bound on the marginal likelihood p(data).
    // We may achieve this by maximizing the expectation of the complete data likelihood
    // with respect to the posterior over z from the E-step
    private void takeMstep() {
        // First we compute the effective counts of each state, N_k in the docs. We do this separately over alt and ref sites
        final double[] effectiveAltCountsFromDesignMatrix = GATKProtectedMathUtils.sumArrayFunction(0, altDesignMatrix.size(), n -> altResponsibilities[n]);
        double[] effectiveAltCountsFromHistograms = new double[NUM_STATES];

        for (Histogram<Integer> h : altHistograms){
            final int artifactIndex = h.getValueLabel().endsWith(ArtifactType.F1R2.toString()) ? 0 : 1;
            effectiveAltCountsFromHistograms = MathArrays.ebeAdd(effectiveAltCountsFromHistograms,
                    GATKProtectedMathUtils.sumArrayFunction(0, maxDepthForHistograms, d ->
                            MathArrays.scale(h.get(d + 1).getValue(), responsibilitiesOfAltDepth1Sites[d][artifactIndex])));
        }
        final double[] effectiveAltCounts = MathArrays.ebeAdd(effectiveAltCountsFromDesignMatrix, effectiveAltCountsFromHistograms);

        Utils.validate(Math.abs(MathUtils.sum(effectiveAltCounts) - numAltExamples) < EPSILON,
                String.format("effective alt counts must add up to %d but got %f", numAltExamples, MathUtils.sum(effectiveAltCounts)));

        // TODO: at some depth, the responsibilities must be 1 for z = hom ref and 0 for everything else, we could probably save some time there
        // Over ref sites, we have a histogram of sites over different depths. At each depth we simply multiply the responsibilities by the number of sites,
        // and sum them over all of depths. Because we cut off the depth histogram at {@code MAX_COVERAGE}, we underestimate the ref effective counts by design
        final double[] effectiveRefCounts = GATKProtectedMathUtils.sumArrayFunction(0, maxDepthForHistograms,
                c -> MathArrays.scale(refHistogram.get(c + 1).getValue(), refResponsibilities[c]));
        Utils.validate(Math.abs(MathUtils.sum(effectiveRefCounts) - numRefExamples) < EPSILON,
                String.format("effective ref counts must add up to %d but got %f", numRefExamples, MathUtils.sum(effectiveRefCounts)));

        effectiveCounts = MathArrays.ebeAdd(effectiveAltCounts, effectiveRefCounts);
        Utils.validate(effectiveCounts.length == NUM_STATES, "effectiveCount must be a k-dimensional vector");
        Utils.validate(Math.abs(MathUtils.sum(effectiveCounts) - numExamples) < EPSILON,
                String.format("effective counts must add up to number of examples %d but got %f", numExamples, MathUtils.sum(effectiveCounts)));

        // Update pi
        // ===========================
        pi = MathArrays.scale(1.0 / numExamples, effectiveCounts);
        Utils.validate(Math.abs(MathUtils.sum(pi) - 1.0) < EPSILON, "pi must be normalized");

        return;
    }

    public static double[] computeLog10Responsibilities(final Nucleotide refAllele, final Nucleotide altAllele,
                                                      final int altDepth, final int f1r2AltCount, final int depth,
                                                      final double[] pi) {
        final double[] log10UnnormalizedResponsibilities = new double[LearnHyperparametersEngine.NUM_STATES];
        List<ArtifactState> impossibleStates = ArtifactState.getImpossibleStates(refAllele);

        for (ArtifactState z : ArtifactState.values()){
            final int k = z.ordinal();
            if (impossibleStates.contains(z)) {
                // this state is impossible e.g. F1R2_G under context AGT and should get the normalized probability of 0
                log10UnnormalizedResponsibilities[k] = Double.NEGATIVE_INFINITY;
                continue;
            }

            if (ArtifactState.artifactStates.contains(z) && ArtifactState.getAltAlleleOfArtifact(z) != altAllele) {
                // artifact states whose alt allele doesn't match up with the observed allele
                log10UnnormalizedResponsibilities[k] = Double.NEGATIVE_INFINITY;
                continue;
            }

            // If we get here, we have a non-artifact state i.e. { germline het, hom ref, hom var, somatic het }
            // or an artifact state whose transitions match the observed alt allele (e.g. alt allele = A, z = F1R2_A, F2R1_A)
            log10UnnormalizedResponsibilities[k] = computePosterior(z, altDepth, f1r2AltCount, depth, pi[k],
                    alleleFractionPseudoCounts.get(k), altF1R2FractionPseudoCounts.get(k));
        }

        return log10UnnormalizedResponsibilities;
    }


    // Compute the posterior probability of the state z given data. The caller is responsible for not calling
    // this method on inconsistent states e.g. z = F1R2_C where the reference context is ACT
    public static double computePosterior(final ArtifactState z, final int altDepth, final int altF1R2Depth, final int depth,
                                          final double pi, final Pair<Double, Double> afPseudoCounts,
                                          final Pair<Double, Double> f1r2PseudoCounts){
        Utils.validateArg(MathUtils.isAProbability(pi), String.format("pi must be a probability but got %f", pi));
        Utils.validateArg(afPseudoCounts.getFirst() > 0 && afPseudoCounts.getSecond() > 0,
                String.format("pseudocounts for allele fraction must be greater than 0 but got %f and %f",
                        afPseudoCounts.getFirst(), afPseudoCounts.getSecond()));
        Utils.validateArg(f1r2PseudoCounts.getFirst() > 0 && f1r2PseudoCounts.getSecond() > 0,
                String.format("pseudocounts for alt F1R2 fraction must be greater than 0 but got %f and %f",
                        f1r2PseudoCounts.getFirst(), f1r2PseudoCounts.getSecond()));

        return Math.log10(pi) +
                MathUtils.log10BetaBinomialProbability(altDepth, depth, afPseudoCounts.getFirst(), afPseudoCounts.getSecond()) +
                MathUtils.log10BetaBinomialProbability(altF1R2Depth, altDepth, f1r2PseudoCounts.getFirst(), f1r2PseudoCounts.getSecond());
    }

    // For each state z, define allele fraction distribution f_z and alt f1r2 fraction theta_z
    // They are both beta binomial distributions and are therefore parameterized by the pseudocounts alpha and beta
    // ===========================
    private static List<Pair<Double, Double>> getPseudoCountsForAlleleFraction(){
        final List<Pair<Double, Double>> alleleFractionPseudoCounts = new ArrayList<>(ArtifactState.values().length);
        final double pseudoCountOfAltUnderArtifact = 1.0; // give it a flat prior and see what happens
        final double pseudoCountOfRefUnderArtifact = 25.0; // alpha = 1, beta = 10 gives a nice distribution that gives high prob to very low allele fractions

        // for home ref sites assume Q35, but maintaining some width
        final double pseudoCountOfHomCorrect = 10000.0;
        final double pseudoCountOfHomError = 3;

        final double balancedPseudoCount = 50;

        // These variables define the distribution of allele fraction in a somatic variant. It should be learned from
        // the data in the future. In the meantime, one should tweak this by hand when e.g. applying the read orientation
        // filter on low allele fraction samples such as the blood biopsy
        final double pseudoCountOfSomaticAlt = 2;
        final double pseudoCountOfSomaticRef = 5;

        // The allele fraction distribution, which is not aware of the read orientation, should be the same between
        // F1R2 and F2R1 artifacts
        for (ArtifactState z : ArtifactState.getF1R2States()){
            alleleFractionPseudoCounts.add(z.ordinal(), new Pair<>(pseudoCountOfAltUnderArtifact, pseudoCountOfRefUnderArtifact));
        }

        for (ArtifactState z : ArtifactState.getF2R1States()){
            alleleFractionPseudoCounts.add(z.ordinal(), new Pair<>(pseudoCountOfAltUnderArtifact, pseudoCountOfRefUnderArtifact));
        }

        for (ArtifactState z : ArtifactState.getNonArtifactStates()){
            switch (z) {
                case HOM_REF:
                    alleleFractionPseudoCounts.add(z.ordinal(), new Pair<>(pseudoCountOfHomError, pseudoCountOfHomCorrect));
                    break;
                case GERMLINE_HET:
                    alleleFractionPseudoCounts.add(z.ordinal(), new Pair<>(balancedPseudoCount, balancedPseudoCount));
                    break;
                case SOMATIC_HET:
                    alleleFractionPseudoCounts.add(z.ordinal(), new Pair<>(pseudoCountOfSomaticAlt, pseudoCountOfSomaticRef));
                    break;
                case HOM_VAR:
                    alleleFractionPseudoCounts.add(z.ordinal(), new Pair<>(pseudoCountOfHomCorrect, pseudoCountOfHomError));
                    break;
            }
        }
        return alleleFractionPseudoCounts;
    }

    private static List<Pair<Double, Double>> getPseudoCountsForAltF1R2Fraction(){
        final List<Pair<Double, Double>> altF1R2FractionPseudoCounts = new ArrayList<>(ArtifactState.values().length);

        final double pseudoCountOfLikelyOutcome = 500.0;
        final double pseudoCountOfRareOutcome = 1.0;
        final double balancedPseudoCount = 50;


        for (ArtifactState z : ArtifactState.getF1R2States()){
            altF1R2FractionPseudoCounts.add(z.ordinal(), new Pair<>(pseudoCountOfLikelyOutcome, pseudoCountOfRareOutcome));
        }

        for (ArtifactState z : ArtifactState.getF2R1States()){
            altF1R2FractionPseudoCounts.add(z.ordinal(), new Pair<>(pseudoCountOfRareOutcome, pseudoCountOfLikelyOutcome));
        }

        for (ArtifactState z : ArtifactState.getNonArtifactStates()){
            altF1R2FractionPseudoCounts.add(z.ordinal(), new Pair<>(balancedPseudoCount, balancedPseudoCount));

        }
        return altF1R2FractionPseudoCounts;
    }

    /**
     * This enum encapsulates the domain of the discrete latent random variable z
     */
    public enum ArtifactState {
        // F1R2 artifact to a particular alt base. The F1R2_{ref} will be ignored (e.g. under the ref context AGT,
        // we ignore F1R2_G
        F1R2_A,
        F1R2_C,
        F1R2_G,
        F1R2_T,
        // F2R1 artifact states. We ignore F2R1_{ref}
        F2R1_A,
        F2R1_C,
        F2R1_G,
        F2R1_T,
        // What follows below are states in which no read orientation artifact is detected
        HOM_REF,
        GERMLINE_HET,
        SOMATIC_HET,
        HOM_VAR;

        public static List<ArtifactState> getStates(){
            return Arrays.stream(ArtifactState.values()).collect(Collectors.toList());
        }

        static ArtifactState[] getF1R2States(){
            return new ArtifactState[]{F1R2_A, F1R2_C, F1R2_G, F1R2_T};
        }

        static ArtifactState[] getF2R1States(){
            return new ArtifactState[]{F2R1_A, F2R1_C, F2R1_G, F2R1_T};
        }

        public static List<ArtifactState> getNonArtifactStates(){
            return Arrays.asList(HOM_REF, GERMLINE_HET, SOMATIC_HET, HOM_VAR);
        }

        public static List<ArtifactState> getImpossibleStates(final Nucleotide refAllele){
            switch (refAllele){
                case A : return Arrays.asList( F1R2_A, F2R1_A );
                case C : return Arrays.asList( F1R2_C, F2R1_C );
                case G : return Arrays.asList( F1R2_G, F2R1_G );
                case T : return Arrays.asList( F1R2_T, F2R1_T );
                default: throw new UserException(String.format("Invalid nucleotide given: %s", refAllele));
            }
        }

        // Given a state z, return the alt allele of the artifact that the state encodes
        public static Nucleotide getAltAlleleOfArtifact(final ArtifactState z){
            Utils.validateArg(Arrays.asList(ArtifactState.F1R2_A, ArtifactState.F1R2_C, ArtifactState.F1R2_G, ArtifactState.F1R2_T,
                    ArtifactState.F2R1_A, ArtifactState.F2R1_C, ArtifactState.F2R1_G, ArtifactState.F2R1_T).contains(z),
                    String.format("ArtifactState must be F1R2_a or F2R1_a but got %s", z));
            switch (z){
                case F1R2_A : return Nucleotide.A;
                case F1R2_C : return Nucleotide.C;
                case F1R2_G : return Nucleotide.G;
                case F1R2_T : return Nucleotide.T;
                case F2R1_A : return Nucleotide.A;
                case F2R1_C : return Nucleotide.C;
                case F2R1_G : return Nucleotide.G;
                case F2R1_T : return Nucleotide.T;
                default: throw new UserException(String.format("Invalid state: %s", z));
            }
        }

        static List<ArtifactState> artifactStates = Arrays.asList(F1R2_A, F1R2_C, F1R2_G, F1R2_T, F2R1_A, F2R1_C, F2R1_G, F2R1_T);


        public static ArtifactState getF1R2StateOfInterest(final Nucleotide altAllele) {
            switch (altAllele) {
                case A:
                    return ArtifactState.F1R2_A;
                case C:
                    return ArtifactState.F1R2_C;
                case G:
                    return ArtifactState.F1R2_G;
                case T:
                    return ArtifactState.F1R2_T;
                default:
                    throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
            }
        }

        public static ArtifactState getF2R1StateOfInterest(final Nucleotide altAllele) {
            switch (altAllele) {
                case A:
                    return ArtifactState.F2R1_A;
                case C:
                    return ArtifactState.F2R1_C;
                case G:
                    return ArtifactState.F2R1_G;
                case T:
                    return ArtifactState.F2R1_T;
                default:
                    throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
            }
        }

        public static ArtifactState getRevCompState(final ArtifactState state){
            switch (state) {
                case F1R2_A:
                    return F2R1_T;
                case F1R2_C:
                    return F2R1_G;
                case F1R2_G:
                    return F2R1_C;
                case F1R2_T:
                    return F2R1_A;
                case F2R1_A:
                    return F1R2_T;
                case F2R1_C:
                    return F1R2_G;
                case F2R1_G:
                    return F1R2_C;
                case F2R1_T:
                    return F1R2_A;
                default:
                    return state;
            }
        }
    }

}
