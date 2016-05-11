package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Holds the core Smith-Waterman alignment parameters of
 * match value, and mismatch, gap open and gap extension penalties.
 */
public final class SmithWatermanParameters {
    public final int w_match;
    public final int w_mismatch;
    public final int w_open;
    public final int w_extend;
    public final OverhangStrategy overhangStrategy;

    /**
     * Create a new set of SW parameters
     *
     * @param w_match    the match score
     * @param w_mismatch the mismatch penalty
     * @param w_open     the gap open penalty
     * @param w_extend   the gap extension penalty
     * @param overhangStrategy What strategy should we use when the best path does not start/end at the corners of the matrix
     */
    public SmithWatermanParameters(final int w_match, final int w_mismatch, final int w_open, final int w_extend, final OverhangStrategy overhangStrategy) {
        if (w_mismatch > 0) throw new IllegalArgumentException("w_mismatch must be <= 0 but got " + w_mismatch);
        if (w_open > 0) throw new IllegalArgumentException("w_open must be <= 0 but got " + w_open);
        if (w_extend > 0) throw new IllegalArgumentException("w_extend must be <= 0 but got " + w_extend);
        Utils.nonNull(overhangStrategy);

        this.w_match = w_match;
        this.w_mismatch = w_mismatch;
        this.w_open = w_open;
        this.w_extend = w_extend;
        this.overhangStrategy = overhangStrategy;
    }

    public int[][] toSubstitutionMatrix() {
        int[][] score = new int[128][128];
        for (int i = 0; i < 128; i++) {
            for (int j = 0; j < 128; j++) {
                score[i][j] = i == j ? w_match : w_mismatch;
            }
        }
        return score;
    }
}
