package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.apache.commons.math3.util.MathArrays;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.spark_project.guava.annotations.VisibleForTesting;

import java.util.*;
import java.util.stream.IntStream;


/**
 * Created by tsato on 7/26/17.
 */
public class LearnHyperparametersEngine {
    public static final int NUM_STATES = State.values().length;

    // When the increase in likelihood falls under this value, we call the algorithm converged
    private final double convergenceThreshold;

    // If the EM does not converge in a few steps we should suspect that something went wrong
    private final int maxEMIterations;

    // Regularizer (?) TODO: think this through
    static final double EPSILON = 1e-4;


    final String referenceContext;

    final Nucleotide refAllele;

    final RefSiteHistogram refHistogram;

    final List<AltSiteRecord> altDesignMatrix;

    // N by K matrix of posterior probabilities of latent variable z, where N is the number of alt sites,
    // evaluated at the current estimates of the hyperparameters pi, f, and theta
    final double[][] altResponsibilities;

    // {@code MAX_COVERAGE} by K matrix of a cache of responsibilities of a ref site (i.e. m = 0, x = 0)
    // for ref sites with coverage 0, 1, ..., MAX_VALUE - 1. To reiterate, the rows represent different coverages,
    // not samples (the count of samples in each coverage is stored in a separate histogram)
    final double[][] refResponsibilities = new double[RefSiteHistogram.MAX_DEPTH][NUM_STATES];

    final int numAltExamples;

    final int numExamples;

    final List<State> impossibleStates;

    // K-dimensional vector of effective sample counts for each class of z, weighted by the the altResponsibilities. For a fixed k,
    // we sum up the counts over all alleles. N_k in the docs.
    // TODO: should be final
    @VisibleForTesting
    double[] effectiveCounts = new double[NUM_STATES];

    /*** Hyperparameters of the model ***/

    // pi is the K-dimensional vector of probabilities for the categorical variable z. Adds up to 1.0
    double[] pi = new double[NUM_STATES];

    // K-dimensional vector of conditional probabilities for the binomial r.v. m given z.
    // m|z represents the number of alt reads at site n given state z
    double[] f = new double[NUM_STATES];

    // In case we observe no data - assume that the allele fraction given artifact is this value
    // TODO: we also use this as the initial value for f - perhaps rename
    static final double DEFAULT_ARTIFACT_ALLELE_FRACTION = 0.3;

    // K-dimensional vector of probabilities for the binomial x given z, which represents the number of F1R2 alt reads at site n
    double[] theta = new double[NUM_STATES];

    // parameters to the beta prior for the allele fraction f given z = somatic het
    public static final int alpha = 3;
    public static final int beta = 7;

    private int numIterations = 0;

    public LearnHyperparametersEngine(final RefSiteHistogram refHistogram, final List<AltSiteRecord> altDesignMatrixForContext,
                                      final double convergenceThreshold, final int maxEMIterations){
        // TODO: count up the number of rows in altDataTable
        // TODO: restructure the code such that learn hyperparameters is a tool but LearnHyperparametersEngine is a separate class with a constructor,
        // and you initialize the object for each reference context

        this.refHistogram = refHistogram;
        this.altDesignMatrix = altDesignMatrixForContext;
        numAltExamples = altDesignMatrix.size();
        numExamples = numAltExamples + IntStream.of(refHistogram.getCounts()).sum();
        altResponsibilities = new double[numAltExamples][NUM_STATES];
        referenceContext = refHistogram.getReferenceContext();
        refAllele = Nucleotide.valueOf(referenceContext.substring(1,2));
        this.convergenceThreshold = convergenceThreshold;
        this.maxEMIterations = maxEMIterations;

        // initialize pi
        // ===========================
        // Some artifact states don't make sense for a given context and therefore should be given the porbability of 0
        // e.g. under the ref context AGT, F1R2_G and F2R1_G states are impossible, because by definition a read orientation
        // artifact only applies to alt sites (REWORD). We would remove those entries from the vector but it simplifies the
        // implementation (i.e. we can share the same indices for the states across all contexts) if we just set some
        // probabilities to 0
        impossibleStates = State.getImpossibleStates(refAllele);
        Arrays.fill(pi, 1.0/(NUM_STATES - impossibleStates.size()));
        for (State impossibleState : impossibleStates){
            pi[impossibleState.ordinal()] = 0;
        }

        // initialize f and theta
        // ===========================
        // we fix some of the parameters to entice the model to assign particular states to the indices into z
        // for instance, we fix the allele fraction parameter f for z = Balanced Het to be 0.5.
        for (State z : State.artifactStates){
            f[z.ordinal()] = DEFAULT_ARTIFACT_ALLELE_FRACTION;
        }
        f[State.GERMLINE_HET.ordinal()] = 0.5;
        f[State.HOM_REF.ordinal()] = EPSILON;
        f[State.HOM_VAR.ordinal()] = 1 - EPSILON;
        f[State.SOMATIC_HET.ordinal()] = 0; // m is parameterized by alpha and beta when z = somatic het

        // similarly, we may fix some of theta_z
        // TODO: or should we learn them?
        for (State state : State.getF1R2States()){
            theta[state.ordinal()] = 1 - EPSILON;
        }

        for (State state : State.getF2R1States()){
            theta[state.ordinal()] = EPSILON;
        }

        for (State state : State.getNonArtifactStates()){
            theta[state.ordinal()] = 0.5;
        }

        // initialize responsibilities TODO: no need to initialize responsibilities
        // ===========================
        // FIXME: revisit after the test - we only need to initialize either resonsibilities (M-step first) or
        // hyperparameters (E-step first). Probably makes sense to initialize hyperparameters.
        final double initialResponsibility = 1.0/NUM_STATES;
        for (int n = 0; n < numAltExamples; n++ ) {
            Arrays.fill(altResponsibilities[n], initialResponsibility);
        }

        /**
         * TODO: explain how we optimize for ref sites
         */
        for (int n = 0; n < RefSiteHistogram.MAX_DEPTH; n++ ) {
            Arrays.fill(refResponsibilities[n], initialResponsibility);
        }
    }

    public Hyperparameters runEMAlgorithm(final Logger logger){
        boolean converged = false;
        double[] oldPi = new double[NUM_STATES];
        // one may plot the changes in L2 distance of parameters to make sure that
        // EM is steadily moving towards the (local? global?) maximum
        double[] l2distancesOfParameters = new double[maxEMIterations];

        while (!converged && numIterations < maxEMIterations){
            // TODO: stylistic problems here, there's too many side-effects
            takeEstep();
            // FIXME: perhaps M-step shoudl return a Hyperparameters data structure
            takeMstep(); // hyperparameters e.g. {@code pi} are updated via side-effect

            // assert newLikelihood >= oldLikelihood : "M step must increase the likelihood";
            // TODO: just pi? what about f and theta?
            final double l2Distance = MathArrays.distance(oldPi, pi);
            converged = l2Distance < convergenceThreshold;

            l2distancesOfParameters[numIterations] = l2Distance;

            oldPi = Arrays.copyOf(pi, NUM_STATES);


            numIterations++;
        }

        logger.info(String.format("Context %s, EM converged in %d steps", referenceContext, numIterations));
        logger.info(String.format("The changes in L2 distance of pi between iterations %s",
                Hyperparameters.doubleArrayToString(l2distancesOfParameters)));
        return new Hyperparameters(referenceContext, pi, f, theta, numExamples, numAltExamples);
    }

    // Given the current estimates of the parameters pi, f, and theta, compute the log10AltResponsibilities
    // gamma_nk = p(z_nk)
    private void takeEstep(){
        // we save some computation here by recognizing that ref sites with the same depth have the same alt depth and
        // alt F1R2 depth (i.e. m = x = 0). Thus responsibilities for ref sites are a function only of the depth (ref and
        // alt combined) and therefore we need only compute the responsibility once for unique depth 0, 1, ..., MAX_COVERAGE
        for (int depth = 0; depth < RefSiteHistogram.MAX_DEPTH; depth++){
            final int r = depth; // another hack to use depth in a stream

            final double[] log10UnnormalizedResponsibilities = new double[NUM_STATES];
            for (State z : State.values()){
                final int k = z.ordinal();
                if (z == State.SOMATIC_HET){
                    // note log(Binom(k=0|n=0, p)) = log(1) = 0, so we may ignore the term coming from the likelihood
                    // of the number of alt F1R2 reads
                    log10UnnormalizedResponsibilities[k] = Math.log10(pi[k]) +
                            MathUtils.log10BetaBinomialDensity(0, r, alpha, beta);
                } else if (impossibleStates.contains(z)){
                    // give 0 probabilities to impossible states
                    log10UnnormalizedResponsibilities[k] = Double.NEGATIVE_INFINITY;
                    continue;
                } else {
                    log10UnnormalizedResponsibilities[k] = Math.log10(pi[k]) +
                            MathUtils.log10BinomialProbability(r, 0, Math.log10(f[k]));
                }
            }

            // TODO: computing responsibilities for ref will not be too expensive, but, as noted elsewhere, it may make sense to peg z = hom ref at 1.0 at certain coverage
            refResponsibilities[r] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);
            assert Math.abs(MathUtils.sum(refResponsibilities[r]) - 1.0) < EPSILON :
                    String.format("ref responsibility for depth = %d added up to %f", r, MathUtils.sum(refResponsibilities[r]));
        }

        // compute the altResponsibilities of each of n alt sites \gamma_{nk}
        for (int n = 0; n < numAltExamples; n++){
            final AltSiteRecord example = altDesignMatrix.get(n);

            final int depth = example.getDepth();
            final int[] baseCounts = example.getBaseCounts();
            final int[] f1r2Counts = example.getF1R2Counts();
            final Nucleotide altAllele = example.getAltAllele();

            final int depthOfMostLikelyAltAllele = baseCounts[altAllele.ordinal()];
            final int f1r2DepthOfMostLikelyAltAllele = f1r2Counts[altAllele.ordinal()];

            // K-dimensional array of one of the terms that comprises gamma*_{nk}
            double[] log10UnnormalizedResponsibilities = new double[NUM_STATES];

            for (State z : State.values()){
                final int k = z.ordinal();
                if (impossibleStates.contains(z)){
                    // this state is impossible e.g. F1R2_G under context AGT and should get the normalized probability of 0
                    log10UnnormalizedResponsibilities[k] = Double.NEGATIVE_INFINITY;
                    continue;
                }

                // breaking Java's variable naming convention to make clear how the variables relate to the docs
                final int m_nk;
                final int x_nk;

                if (State.artifactStates.contains(z)){
                    // for an artifact state, use the corresponding alt base count of the transition in question
                    // e.g. for z = F1R2_c use the count of C's
                    final Nucleotide altAlleleOfTransition = State.getAltAlleleOfTransition(z);
                    m_nk = baseCounts[altAlleleOfTransition.ordinal()];
                    x_nk = f1r2Counts[altAlleleOfTransition.ordinal()];
                } else {
                    // we have a non-artifact state i.e. { germline het, hom ref, hom var, somatic het }
                    m_nk = depthOfMostLikelyAltAllele;
                    x_nk = f1r2DepthOfMostLikelyAltAllele;
                }

                log10UnnormalizedResponsibilities[k] = computePosterior(z, m_nk, x_nk, depth, pi[k], f[k], theta[k]);

            }


            // TODO: maybe we don't need to keep the matrix of log10 alt responsibilities
            // we normalize responsibilities here because the M-step uses normalized responsibilities
            altResponsibilities[n] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);

            assert Math.abs(MathUtils.sum(altResponsibilities[n]) - 1.0) < EPSILON :
                    String.format("responsibility for %dth example added up to %f", n,  MathUtils.sumLog10(altResponsibilities[n]));
        }
    }

    // given the current posterior distributions over z (aka repsonsibilities), compute the estimate for
    // the probabilities for the categorial variable z (pi), allele fractions (f), and alt F1R2 fraction (theta) that maximize the lower bound
    // on the marginal likelihood p(data). We may achieve this by maximizing the expectation of the complete data likelihood
    // with respect to the posterior over z evaluated at the old parameter estimates
    private void takeMstep(){
        /*** compute responsibility-based statistics based on the current log10AltResponsibilities ***/

        // First we compute the effective counts of each state, N_k in the docs. We do so separately over alt and ref sites
        final double[] effectiveAltCounts = GATKProtectedMathUtils.sumArrayFunction(0, numAltExamples, n -> altResponsibilities[n]);

        // TODO: at some depth, the responsibilities must be 1 for z = hom ref and 0 for everything else, we could probably save some time there
        // Over ref sites, we have a histogram of sites over different depths. At each depth we simply multiply the responsibilities by the number of sites,
        // and sum them over all of depths. Because we cut off the depth histogram at {@code MAX_COVERAGE}, we underestimate the ref effective counts by design
        final double[] effectiveRefCounts = GATKProtectedMathUtils.sumArrayFunction(0, RefSiteHistogram.MAX_DEPTH, c ->
                        MathArrays.scale((double)refHistogram.getCounts()[c], refResponsibilities[c]));
        effectiveCounts = MathArrays.ebeAdd(effectiveAltCounts, effectiveRefCounts);

        assert effectiveCounts.length == NUM_STATES : "effectiveCount must be a k-dimensional vector";
        assert Math.abs(MathUtils.sum(effectiveCounts) - numExamples) < EPSILON :
                String.format("effective counts must add up to number of examples %d but got %f", numExamples, MathUtils.sum(effectiveCounts));

        // K-dimensional vector of mean alt depth weighted by the responsibilities, times the effective count. N_k \bar{m} in the docs
        // ref sites do not contribute to the sum as m_nk = 0 for ref sites
        final double[] weightedEffectiveAltDepth = GATKProtectedMathUtils.sumArrayFunction(0, numAltExamples,
                n -> MathArrays.ebeMultiply(getAltCountsArrangedByState(n, ObservedData.ALT_DEPTH), altResponsibilities[n]));
        assert weightedEffectiveAltDepth.length == NUM_STATES : "weightedAvgAltDepth should have length K";

        // K-dimensional vector of mean alt F1R2 depth weighted by the responsibilityes, times the effective count. N_k \bar{x} in the docs
        // ref sites do not contribute to the sum as x_n = 0 for ref sites
        final double[] weightedEffectiveAltF1R2Depth = GATKProtectedMathUtils.sumArrayFunction(0, numAltExamples,
                n -> MathArrays.ebeMultiply(getAltCountsArrangedByState(n, ObservedData.ALT_F1R2_DEPTH), altResponsibilities[n]));
        assert weightedEffectiveAltF1R2Depth.length == NUM_STATES : "weightedAvgAltF1R2Depth should have length K";

        final double[] weightedEffectiveDepthOverAltSites = GATKProtectedMathUtils.sumArrayFunction(0, numAltExamples,
                n -> MathArrays.scale((double) altDesignMatrix.get(n).getDepth(), altResponsibilities[n]));
        final double[] weightedAvgDepthOverRefSites = GATKProtectedMathUtils.sumArrayFunction(0, RefSiteHistogram.MAX_DEPTH, c ->
                {
                    final double numSitesAtThisCoverage = (double) refHistogram.getCounts()[c];
                    // TODO: coverage should really start at one - revisit
                    final double[] coverageWeightedByResponsibilities = MathArrays.scale((double) c, refResponsibilities[c]);
                    return MathArrays.scale(numSitesAtThisCoverage, coverageWeightedByResponsibilities);
                });
        // K-dimensional vector of avg depth weighted by the responsibilities, times the effective count. N_k \bar{r} in the docs
        final double[] weightedEffectiveDepth = MathArrays.ebeAdd(weightedEffectiveDepthOverAltSites, weightedAvgDepthOverRefSites);
        assert weightedEffectiveDepth.length == NUM_STATES : "weightedAvgAltDepth should have length K";

        // FIXME: placeholder for learning theta

        // update the vector of allele fractions f
        // ===========================
        // final double[] updatedAlleleFractions = MathArrays.ebeDivide(weightedAvgAltDepth, weightedAvgDepth);
        for (State z : State.values()){
            final int k = z.ordinal();
            if (impossibleStates.contains(z)){
                f[k] = 0;
            } else {
                // if the effective alt count is essentially 0 it could still give us a false but believable allele fraction
                // like 0.2 if the denominator is small too. In such a case we just fix the allele fraction to 0
                f[k] = weightedEffectiveAltDepth[k] < EPSILON ? 0.0 : weightedEffectiveAltDepth[k] / weightedEffectiveDepth[k];
            }

        }


        // FIXME: this part of the code needs some attention
        // For some states we fix f e.g. f[GERMLINE_HET] = 0.5. Also, f[SOMATIC_HET] should be ignored
        // FIXME: Alternatively, for some states we do not update when the value exceeds certain limit e.g. If f[HOM_VAR] deviates
        // an epsilon = 0.05 away from what we expect it to be --- but for now let's just forge ahead.
        final double delta = 0.05;
        // For now, let's not fix anything. TODO: see how the results turn out. Will fix as needed.

        // If anything becomes 0 i.e. we didn't observe a particular alt allele under this context, we just set the value
        // back to the default


        // update pi
        // ===========================
        pi = MathArrays.scale(1.0 / numExamples, effectiveCounts);
        assert  Math.abs(MathUtils.sum(pi) - 1.0) < EPSILON : "pi must be normalized";

        return;
    }

    // Compute the posterior probability of the state z given data. The caller is responsible for not calling
    // this method on inconsistent states e.g. z = F1R2_C where the reference context is ACT
    public static double computePosterior(final State z, final int altDepth, final int altF1R2Depth, final int depth,
                                          final double pi, final double f, final double theta){
        Utils.validate(MathUtils.isAProbability(f), String.format("f must be a probability but got %f", f));
        Utils.validate(MathUtils.isAProbability(pi), String.format("pi must be a probability but got %f", pi));
        Utils.validate(MathUtils.isAProbability(theta), String.format("theta must be a probability but got %f", theta));

        if (z == State.SOMATIC_HET){
            return Math.log10(pi) +
                    MathUtils.log10BetaBinomialDensity(altDepth, depth, alpha, beta) +
                    MathUtils.log10BinomialProbability(altDepth, altF1R2Depth, Math.log10(theta));
        } else {
            return Math.log10(pi) +
                    MathUtils.log10BinomialProbability(depth, altDepth, Math.log10(f)) +
                    MathUtils.log10BinomialProbability(altDepth, altF1R2Depth, Math.log10(theta));
        }
    }

    private enum ObservedData {
        ALT_F1R2_DEPTH, ALT_DEPTH
    }

    // Return K-dimensional array of alt counts m_nk. Each state k uses a different value of alt depth
    // e.g. Z = F1R2_A uses the base count of A's, and Z = F1R2_T uses the count of T's.
    // Returns a double[], not int[], because the downstream API (MathArrays.ebeMultiply()) requires an array of doubles
    private double[] getAltCountsArrangedByState(final int n, final ObservedData observedData){
        final AltSiteRecord altSiteRecord = altDesignMatrix.get(n);
        final Nucleotide refAllele = Nucleotide.valueOf(referenceContext.substring(1, 2));
        final Nucleotide altAllele = altSiteRecord.getAltAllele();
        final double[] arrangedByState = new double[State.values().length];

        int[] counts;
        switch (observedData) {
            case ALT_DEPTH: counts = altSiteRecord.getBaseCounts(); break;
            case ALT_F1R2_DEPTH: counts = altSiteRecord.getF1R2Counts(); break;
            default : throw new UserException("Should be impossible to get here");
        }

        final double countOfAltAllele = counts[altAllele.ordinal()];
        final List<State> impossibleStates = State.getImpossibleStates(refAllele);

        // for the non-artifact states, set m_nk equal to the count of the most likely alt allele
        for (State z : State.getNonArtifactStates()){
            arrangedByState[z.ordinal()] = countOfAltAllele;
        }

        for (State z : State.artifactStates){
            if (impossibleStates.contains(z)){
                continue;
            }

            // for artifact states, get the count of the allele that the artifact e.g. the count of A's for Z=F1R2_A
            arrangedByState[z.ordinal()] = counts[State.getAltAlleleOfTransition(z).ordinal()];
        }
        return arrangedByState;
    }

    private boolean checkLikelihoodHasConverged(final double oldLikelihood, final double newLikelihood){
        return Math.abs(newLikelihood - oldLikelihood) < convergenceThreshold;
    }

    private static boolean checkLikelihoodHasConverged(final int numIterations){
        return numIterations > 10;
    }

    /**
     * This enum encapsulates the domain of the discrete latent random variable z
     */
    public enum State {
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

        static State[] getF1R2States(){
            return new State[]{F1R2_A, F1R2_C, F1R2_G, F1R2_T};
        }

        static State[] getF2R1States(){
            return new State[]{F2R1_A, F2R1_C, F2R1_G, F2R1_T};
        }

        public static List<State> getNonArtifactStates(){
            return Arrays.asList(HOM_REF, GERMLINE_HET, SOMATIC_HET, HOM_VAR);
        }

        public static List<State> getImpossibleStates(final Nucleotide refAllele){
            switch (refAllele){
                case A : return Arrays.asList( F1R2_A, F2R1_A );
                case C : return Arrays.asList( F1R2_C, F2R1_C );
                case G : return Arrays.asList( F1R2_G, F2R1_G );
                case T : return Arrays.asList( F1R2_T, F2R1_T );
                default: throw new UserException(String.format("Invalid nucleotide given: %s", refAllele));
            }
        }

        // Given a state z, return the alt allele of the artifact that the state encodes
        public static Nucleotide getAltAlleleOfTransition(final State z){
            Utils.validateArg(Arrays.asList(State.F1R2_A, State.F1R2_C, State.F1R2_G, State.F1R2_T,
                    State.F2R1_A, State.F2R1_C, State.F2R1_G, State.F2R1_T).contains(z),
                    String.format("State must be F1R2_a or F2R1_a but got %s", z));
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

        static List<State> artifactStates = Arrays.asList(F1R2_A, F1R2_C, F1R2_G, F1R2_T, F2R1_A, F2R1_C, F2R1_G, F2R1_T);


        /** Contract: this method always returns the F1R2 state in the 0th index, F2R1 state in the 1st index **/
        public static State[] getArtifactStatesOfInterest(final Nucleotide altAllele) {
            final State[] relevantArtifactStates;
            switch (altAllele) {
                case A:
                    relevantArtifactStates = new State[]{State.F1R2_A, State.F2R1_A};
                    break;
                case C:
                    relevantArtifactStates = new State[]{State.F1R2_C, State.F2R1_C};
                    break;
                case G:
                    relevantArtifactStates = new State[]{State.F1R2_G, State.F2R1_G};
                    break;
                case T:
                    relevantArtifactStates = new State[]{State.F1R2_T, State.F2R1_T};
                    break;
                default:
                    throw new UserException(String.format("Alt allele must be in {A, C, G, T} but got %s", altAllele));
            }

            Utils.validate(relevantArtifactStates.length == 2, "there should only be 2 artifact states of interest for each alt allele");
            return relevantArtifactStates;
        }
    }



}
