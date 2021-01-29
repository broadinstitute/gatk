package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;

public class VariationalContamination {
    private static int numIterations = 10;
    private double errorRate;

    BetaDistributionShape contamination = BetaDistributionShape.FLAT_BETA;
    private static final double ALPHA_PSEUDOCOUNT = 1.0;
    private static final double BETA_PSEUDOCOUNT = 1.0;

    public VariationalContamination(final double errorRate){
        this.errorRate = errorRate;
    }

    public Pair<Double, Double> calculateContaminationFromHoms(List<PileupSummary> sites){
        int allBases = sites.stream().mapToInt(s -> s.getTotalCount()).sum();
        if (errorRate == 0.0){
            errorRate = 1.0/allBases; // If the reported sequencing error rate is 0
        }

        // Variational E-step
        double effectiveContaminationBaseCount = 0.0;
        double effectiveGoodBaseCount = 0.0;

        // to check for convergence
        double c = -1.0;
        double std = -1.0;
        final double CONVERGENCE_THRESHOLD = 1e-4;

        for (int i = 0; i < numIterations; i++){
            // Reset
            effectiveContaminationBaseCount = 0.0;
            effectiveGoodBaseCount = 0.0;

            // Variational E-step
            for (PileupSummary site : sites){
                final double f = site.getAlleleFrequency();
                final int refCount = site.getRefCount();
                final int altCount = site.getAltCount();
                final int otherAltCount = site.getOtherAltCount();

                // I believe we don't need to worry about losing precision since we normalize for each base
                double effectiveContaminationBaseCountForSite = 0.0;
                double effectiveGoodBaseCountForSite = 0.0;


                final Pair<Double, Double> refEffectiveCounts = getEffectiveCounts(refCount, ReadAllele.REF, f);
                effectiveContaminationBaseCountForSite += refEffectiveCounts.getLeft();
                effectiveGoodBaseCountForSite += refEffectiveCounts.getRight();

                final Pair<Double, Double> altEffectiveCounts = getEffectiveCounts(altCount, ReadAllele.ALT, f);
                effectiveContaminationBaseCountForSite += altEffectiveCounts.getLeft();
                effectiveGoodBaseCountForSite += altEffectiveCounts.getRight();

                final Pair<Double, Double> otherAltEffectiveCounts = getEffectiveCounts(otherAltCount, ReadAllele.OTHER_ALT, f);
                effectiveContaminationBaseCountForSite += otherAltEffectiveCounts.getLeft();
                effectiveGoodBaseCountForSite += otherAltEffectiveCounts.getRight();

                effectiveContaminationBaseCount += effectiveContaminationBaseCountForSite;
                effectiveGoodBaseCount += effectiveGoodBaseCountForSite;
                if (Math.abs(effectiveContaminationBaseCountForSite + effectiveGoodBaseCountForSite - site.getTotalCount()) > 1e-3){
                    throw new UserException("Effective Counts don't add up.");
                }
            }

            // Variaitonal M-step
            contamination = new BetaDistributionShape(ALPHA_PSEUDOCOUNT + effectiveContaminationBaseCount, BETA_PSEUDOCOUNT + effectiveGoodBaseCount);

            // Sanity check
            final double totalEffectiveCount = effectiveContaminationBaseCount + effectiveGoodBaseCount;
            final double diff = totalEffectiveCount - allBases;
            if (Math.abs(diff) > 1.0){
                throw new UserException("Effective Counts don't add up.");
            }

            // Convergence check
            if (Math.abs(c - contamination.getMean()) < CONVERGENCE_THRESHOLD){
                break;
            }

            c = contamination.getMean();
            std = Math.sqrt(contamination.getVariance());
        }

        return new ImmutablePair<>(contamination.getMean(), Math.sqrt(contamination.getVariance()));
    }

    public enum State {
        CONTAMINATION(0), GOOD(1);

        int index = -1;
        State(int index){
            this.index = index;
        }
    }

    public enum ReadAllele {
        REF, ALT, OTHER_ALT
    }

    public Pair<Double, Double> getEffectiveCounts(int count, ReadAllele allele, double f){
        final int CONTAMINTION_INDEX = 0;
        final int GOOD_INDEX = 1;
        double logAlleleLikelihoodGivenContamination = getConditionalLogProbability(allele, State.CONTAMINATION, f);
        double logAlleleLikelihoodGivenGood = getConditionalLogProbability(allele, State.GOOD, f);
        double[] responsibility = MathUtils.normalizeFromLog10ToLinearSpace(new double[]{
                contamination.getExpectationOfLog() + logAlleleLikelihoodGivenContamination, contamination.getExpectationOfLog1MinusP() + logAlleleLikelihoodGivenGood });
        double effectiveContamCount = count * responsibility[CONTAMINTION_INDEX];
        double effectiveGoodCount = count * responsibility[GOOD_INDEX];
        if (Math.abs(effectiveContamCount + effectiveGoodCount - count) > 1e-3){
            throw new UserException("Effective counts don't add up");
        }
        return new ImmutablePair<>(effectiveContamCount, effectiveGoodCount);
    }


    /**
     *
     * Computes log p(x=allle|z=state)
     * @param f allele frequency at the site
     */
    public double getConditionalLogProbability(final ReadAllele readAllele, final State state, final double f){
        if (readAllele == ReadAllele.REF){
            return state == State.CONTAMINATION ? Math.log((1.0-f)*(1-errorRate) + f*errorRate/3.0) : Math.log(errorRate/3.0);
        } else if (readAllele == ReadAllele.ALT) {
            return state == State.CONTAMINATION ? Math.log(f*(1-errorRate) + (1-f)*errorRate/3.0) : Math.log(1.0-errorRate);
        } else {
            return Math.log(2.0 * errorRate / 3.0);
        }
    }
}
