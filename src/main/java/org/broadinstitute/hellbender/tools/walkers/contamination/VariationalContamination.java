package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;

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
            errorRate = 1.0/allBases; // If the reported error rate is none
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
                double refGivenContam = getConditionalLogProbability(ReadAllele.REF, State.CONTAMINATION, f);
                double refGivenGood = getConditionalLogProbability(ReadAllele.REF, State.GOOD, f);
                double unnormalizedResponsibilityContamRef = Math.exp(contamination.getLogMean() + refGivenContam);
                double unnormalizedResponsibilityGoodRef = Math.exp(contamination.getExpectationLog1MinusP() + refGivenGood);
                double responsibilityRefBase = unnormalizedResponsibilityContamRef / (unnormalizedResponsibilityContamRef + unnormalizedResponsibilityGoodRef);
                double effectiveContamRefBaseCountForSite = refCount * responsibilityRefBase;
                double effectiveGoodRefBaseCountForSite = refCount * (1-responsibilityRefBase);
                if (Math.abs(effectiveContamRefBaseCountForSite + effectiveGoodRefBaseCountForSite - refCount) > 1e-3){
                    throw new UserException("Ref effective counts don't add up");
                }
                effectiveContaminationBaseCount += effectiveContamRefBaseCountForSite;
                effectiveGoodBaseCount += effectiveGoodRefBaseCountForSite;

                double altGivenContam = getConditionalLogProbability(ReadAllele.ALT, State.CONTAMINATION, f);
                double altGivenGood = getConditionalLogProbability(ReadAllele.ALT, State.GOOD, f);
                double unnormalizedResponsibilityContamAlt = Math.exp(contamination.getLogMean() + altGivenContam);
                double unnormalizedResponsibilityGoodAlt = Math.exp(contamination.getExpectationLog1MinusP() + altGivenGood);
                double responsibilityAltBase = unnormalizedResponsibilityContamAlt / (unnormalizedResponsibilityContamAlt + unnormalizedResponsibilityGoodAlt);
                double effectiveContamAltBaseCountForSite = altCount * responsibilityAltBase;
                double effectiveGoodAltBaseCountForSite = altCount * (1-responsibilityAltBase);
                if (Math.abs(effectiveContamAltBaseCountForSite + effectiveGoodAltBaseCountForSite - altCount) > 1e-3){
                    throw new UserException("Alt effective counts don't add up");
                }
                effectiveContaminationBaseCount += effectiveContamAltBaseCountForSite;
                effectiveGoodBaseCount += effectiveGoodAltBaseCountForSite;

                double otherAltGivenContam = getConditionalLogProbability(ReadAllele.OTHER_ALT, State.CONTAMINATION, f);
                double otherAltGivenGood = getConditionalLogProbability(ReadAllele.OTHER_ALT, State.GOOD, f);
                double unnormalizedResponsibilityContamOther = Math.exp(contamination.getLogMean() + otherAltGivenContam);
                double unnormalizedResponsibilityGoodOther = Math.exp(contamination.getLogMean() + otherAltGivenGood);
                double responsibilityOtherAltBase = unnormalizedResponsibilityContamOther / (unnormalizedResponsibilityContamOther + unnormalizedResponsibilityGoodOther);
                double effectiveContamOtherAltBaseCountForSite = otherAltCount * responsibilityOtherAltBase;
                double effectiveGoodOtherAltBaseCountForSite = otherAltCount * (1.0 - responsibilityOtherAltBase);
                if (Math.abs(effectiveContamOtherAltBaseCountForSite + effectiveGoodOtherAltBaseCountForSite - otherAltCount) > 1e-3){
                    throw new UserException("Other alt effective counts don't add up");
                }

                effectiveContaminationBaseCount += effectiveContamOtherAltBaseCountForSite;
                effectiveGoodBaseCount += effectiveGoodOtherAltBaseCountForSite;
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

    private enum State {
        CONTAMINATION, GOOD
    }

    private enum ReadAllele {
        REF, ALT, OTHER_ALT
    }


    public double getConditionalLogProbability(final ReadAllele readAllele, final State state, final double af){
        if (readAllele == ReadAllele.REF){
            return state == State.CONTAMINATION ? Math.log((1.0-af)*(1-errorRate) + af*errorRate/3.0) : Math.log(errorRate/3.0);
        } else if (readAllele == ReadAllele.ALT) {
            return state == State.CONTAMINATION ? Math.log(af*(1-errorRate) + (1-af)*errorRate/3.0) : Math.log(1.0-errorRate);
        } else {
            return Math.log(2.0 * errorRate / 3.0);
        }
    }
}
