package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

public class SVGraphEdgePrior {

    double[] logPrior;
    private final double BREAKPOINT_EVIDENCE_COUNT_FACTOR = 0.5;
    private final double VARIANT_CALL_EVIDENCE_COUNT_FACTOR = 0.5;
    private final double EVIDENCE_TARGET_LINK_READ_PAIR_FACTOR = 0.5;
    private final double EVIDENCE_TARGET_LINK_SPLIT_READ_FACTOR = 1.0 - EVIDENCE_TARGET_LINK_READ_PAIR_FACTOR;
    private final double EVIDENCE_COUNT_STD_FACTOR = 0.5;

    public SVGraphEdgePrior(final SVGraphEdgeEvidence evidence, final double meanDepth, final int maxVisits) {

        logPrior = new double[maxVisits + 1];
        if (evidence.getEvidenceTargetLink() == null && evidence.getCalledVariant() == null && evidence.getBreakpointPair() == null) {
            for (int i = 0; i < logPrior.length; i++) {
                logPrior[i] = 0;
            }
        } else {
            final int evidenceCount;
            if (evidence.getBreakpointPair() != null) {
                evidenceCount = (int) Math.round(meanDepth * BREAKPOINT_EVIDENCE_COUNT_FACTOR);
            } else if (evidence.getCalledVariant() != null) {
                evidenceCount = (int) Math.round(meanDepth * VARIANT_CALL_EVIDENCE_COUNT_FACTOR);
            } else if (evidence.getEvidenceTargetLink() != null) {
                evidenceCount = Math.max(evidence.getEvidenceTargetLink().getReadPairs(), evidence.getEvidenceTargetLink().getSplitReads());
            } else {
                evidenceCount = 0;
            }
            for (int i = 0; i <= maxVisits; i++) {
                final double std = EVIDENCE_COUNT_STD_FACTOR * meanDepth;
                final double mean = i * meanDepth;
                logPrior[i] = 1; //TODO (new NormalDistribution(mean,std)).logDensity(evidenceCount);
                //final double lambda = i > 0 ? i * meanDepth / 2.0 : 1;
                //logPrior[i] = (new PoissonDistribution(lambda)).logProbability(evidenceCount);
            }
        }
    }

    public double getLogPrior(int n) {
        return logPrior[n];
    }
}
