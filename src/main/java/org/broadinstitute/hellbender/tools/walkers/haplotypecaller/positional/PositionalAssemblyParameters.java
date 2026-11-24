package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;

/**
 * Parameters of the positional assembler.
 *
 * @param k kmer size
 * @param delta how far (in bases) a repeated kmer may join a node from where the read expects it
 * @param maxJump how far a kmer unique in the graph may join a node from where the read expects it ...
 * @param lookahead ... and, beyond delta, only if the read then matches the reference for this many more bases
 * @param minSupport reads needed for an event, fragments needed for an allele read pattern or a new haplotype
 * @param adaptivePruning whether an event also needs the log odds against sequencing error that adaptive pruning
 *                        requires of a graph chain
 * @param initialErrorRate the error rate with which adaptive pruning first finds likely errors, from which it estimates
 *                         the region's error rate
 * @param pruningLogOddsThreshold the natural log odds that adaptive pruning requires
 * @param minBaseQuality bases of lower quality are not used in kmers and don't count as mismatches
 * @param minDanglingBranchLength a read's dangling end with fewer kmers than this, like a de Bruijn graph's dangling
 *                                branch, can't make an event by itself; it only supports events that other reads make
 * @param maxSiteMismatches a read with more mismatches than this against every allele of a site is uninformative there
 * @param maxTailMismatches a dangling end completed against the reference may carry at most this many events
 * @param eventScores alignment scores for decomposing read sequence into events
 */
record PositionalAssemblyParameters(int k, double delta, int maxJump, int lookahead, int minSupport, boolean adaptivePruning,
                                    double initialErrorRate, double pruningLogOddsThreshold, int minBaseQuality,
                                    int minDanglingBranchLength, int maxSiteMismatches, int maxTailMismatches, SWParameters eventScores) {
    static final double DEFAULT_DELTA = 40;
    static final int DEFAULT_MAX_JUMP = 300;
    static final int DEFAULT_LOOKAHEAD = 10;
    static final int DEFAULT_MAX_SITE_MISMATCHES = 2;
    static final int DEFAULT_MAX_TAIL_MISMATCHES = 2;

    PositionalAssemblyParameters {
        if (k < 5 || k > 31) {
            throw new IllegalArgumentException("positional assembly kmer size must be between 5 and 31, not " + k);
        }
        if (minSupport < 1) {
            throw new IllegalArgumentException("positional assembly needs a minimum support (--min-pruning) of at least 1, not " + minSupport);
        }
        if (adaptivePruning && !(initialErrorRate > 0 && initialErrorRate < 1)) {
            throw new IllegalArgumentException("adaptive pruning needs an initial error rate between 0 and 1, not " + initialErrorRate);
        }
        if (minDanglingBranchLength < 0) {
            throw new IllegalArgumentException("the minimum dangling branch length must be non-negative, not " + minDanglingBranchLength);
        }
    }

    static PositionalAssemblyParameters defaults(final int k, final int minSupport, final boolean adaptivePruning, final double initialErrorRate,
                                                 final double pruningLogOddsThreshold, final int minBaseQuality, final int minDanglingBranchLength,
                                                 final SWParameters eventScores) {
        return new PositionalAssemblyParameters(k, DEFAULT_DELTA, DEFAULT_MAX_JUMP, DEFAULT_LOOKAHEAD, minSupport, adaptivePruning,
                initialErrorRate, pruningLogOddsThreshold, minBaseQuality, minDanglingBranchLength, DEFAULT_MAX_SITE_MISMATCHES,
                DEFAULT_MAX_TAIL_MISMATCHES, eventScores);
    }
}
