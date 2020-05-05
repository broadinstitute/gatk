package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.lang.math.IntRange;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.EstimateDragstrParameters;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.stream.IntStream;

public class DragstrParametersEstimator {

    private static final Logger logger = LogManager.getLogger(DragstrParametersEstimator.class);

    private final double[] phred_gp_range;
    private final double[] phred_api_range;
    private final double[] log10_gp_range;
    private final double[] log10_api_range;
    private final double het_hom_ratio;
    private final double log10_het_hom_ratio;
    private final int min_loci_count;
    private static final double LOG_10_OF_2 = Math.log10(2);
    private final double api_mono_threshold;
    private final double min_gop;
    private final double max_gop;
    private final int[] min_gp_idx_by_period;
    private final double[][][] log10_p_err_by_len;
    private final double[][][] log10_p_no_err_by_len;
    private final int maxPeriod;
    private final int maxRepeats;

    public DragstrParametersEstimator(final DragstrCasesSamplerArgumentCollection argumentCollection) {
        maxPeriod = argumentCollection.maxPeriod;
        maxRepeats = argumentCollection.maxRepeats;
        phred_gp_range = argumentCollection.phredGpValues.toDoubleArray();
        phred_api_range = argumentCollection.phredApiValues.toDoubleArray();
        log10_gp_range = MathUtils.applyToArray(phred_gp_range, d -> -.1 * d);
        log10_api_range = MathUtils.applyToArray(phred_api_range, d -> -.1 * d);
        het_hom_ratio = argumentCollection.hetToHomRatio;
        log10_het_hom_ratio = Math.log10(het_hom_ratio);
        min_loci_count = argumentCollection.minLociCount;
        api_mono_threshold = argumentCollection.apiMonothresh;
        min_gop = argumentCollection.minGOP;
        max_gop = argumentCollection.maxGOP;
        log10_p_err_by_len = new double[phred_gp_range.length][maxPeriod][maxRepeats];
        log10_p_no_err_by_len = new double[phred_gp_range.length][maxPeriod][maxRepeats];
        min_gp_idx_by_period = new int[maxPeriod];
        for (int i = 0; i < phred_gp_range.length; i++) {
            for (int k = 0; k < maxPeriod; k++) {
                final int period = k + 1;
                final double log10_p_no_err_per_pos = MathUtils.log10OneMinusPow10(LOG_10_OF_2 + log10_gp_range[i]);
                for (int j = 0; j < maxRepeats; j++) {
                    final int repeats = j + 1;
                    final int lengthInBases = repeats * period;
                    log10_p_no_err_by_len[i][k][j] = lengthInBases * log10_p_no_err_per_pos;
                    log10_p_err_by_len[i][k][j] = MathUtils.log10OneMinusPow10(log10_p_no_err_by_len[i][k][j]);
                }
            }
        }
        for (int i = 0; i < maxPeriod; i++) {
            final int period = i + 1;
            final double gp_min = Math.ceil(-10 * Math.log10((1 - Math.pow(0.5, (1.0/ (20.0 * period)) / 2.0))));
            final int index = Arrays.binarySearch(phred_gp_range, gp_min);
            // since we are looking for a double, we have to be a bit tolerant in terms of differences,
            // so if no exact was found we look the position before the insertion position in case that value
            // is close enough (less than 0.001 way).
            min_gp_idx_by_period[i] = index >= 0 ? index :
                    (index < -1 && Math.abs(gp_min - phred_api_range[-index - 2]) < 0.001) ? -index - 2 : -index - 1;
        }
    }


    public DragstrParams estimateSequencial(final EstimateDragstrParameters.StratifiedDragstrLocusCases cases) {
        final DragstrParamsBuilder builder = new DragstrParamsBuilder(maxPeriod, maxRepeats);
        for (int period = 1; period <= maxPeriod; period++) {
            estimatePeriod(period, builder, cases);
        }
        return builder.make(min_gop, max_gop, .25);
    }

    public DragstrParams estimateParallel(final EstimateDragstrParameters.StratifiedDragstrLocusCases cases, final int threads) {
        final DragstrParamsBuilder builder = new DragstrParamsBuilder(maxPeriod, maxRepeats);
        final int actualThreads = Math.min(threads, maxPeriod + 1);
        Utils.runInParallel(actualThreads, () -> IntStream.range(1, maxPeriod + 1)
                .parallel()
                .forEach(period  -> estimatePeriod(period, builder, cases) ));
        return builder.make(min_gop, max_gop, .25);
    }

    private void estimatePeriod(final int period, final DragstrParamsBuilder builder, final EstimateDragstrParameters.StratifiedDragstrLocusCases cases) {
        // We calculate the "flanks" repeat counts that have too little data
        // to be analyzed separately:
        //cases.removeLowDepth(minimum_depth);
        int accum = 0;
        int leftFlank;
        for (leftFlank = 0; leftFlank < maxRepeats; leftFlank++) {
            if ((accum += cases.get(period, leftFlank + 1).size()) >= min_loci_count) {
                break;
            }
        }
        int rightFlank; accum = 0;
        for (rightFlank = maxRepeats - 1; rightFlank >= 0; rightFlank--) {
            if ((accum += cases.get(period, rightFlank + 1).size()) >= min_loci_count) {
                break;
            }
        }
        if (rightFlank <= leftFlank) {
            throw new UserException.BadInput("not enough cases per period " + period);
        }
        //TODO objection: Taken as done in DRAGstr however...
        //                this code wont combine repeat counts that have less than min_loci_count
        //                between flanks. Inf act leftFlank + 1 may well have less than min_loci_count.
        //                would simply make more sense to keep combining repeat counts to be analyzed
        //                together?
        final ArrayDeque<IntRange> pending = new ArrayDeque<>(maxRepeats);
        final ArrayDeque<IntRange> done = new ArrayDeque<>(maxRepeats);

        pending.add(new IntRange(1, leftFlank + 1));
        for (int i = leftFlank + 1; i < rightFlank; i++) {
            pending.add(new IntRange(i + 1));
        }
        pending.add(new IntRange(rightFlank + 1, maxRepeats));
        IntRange last = null;
        do {
            final IntRange next = pending.pop();
            if (next.getMinimumInteger() != next.getMaximumInteger()) {
                logger.debug("calculating parameters for cases with repeat count " + next.getMinimumInteger() + " to " + next.getMaximumNumber());
            } else {
                logger.debug("calculating parameters for cases with repeat count " + next.getMinimumInteger());
            }
            estimatePeriodRepeatInterval(period, next, builder, cases);
            final double gp1 = builder.gp[period - 1][next.getMinimumInteger() - 1];
            final double api1 = builder.api[period - 1][next.getMinimumInteger() - 1];
            if (done.isEmpty() || (builder.gp[period - 1][last.getMaximumInteger() - 1] >= gp1 &&
                                   builder.api[period - 1][last.getMaximumInteger() - 1] + api_mono_threshold >= api1)) {
                done.addLast(last = next);
            } else {
                pending.push(new IntRange(last.getMinimumNumber(), next.getMaximumNumber()));
                done.removeLast();
                last = !done.isEmpty() ? done.getLast() : null;
            }
        } while (!pending.isEmpty());
    }

    // Given a observed het/hom ratio the total number of reads, the variable ones and the length of the STR in bases.
    // calculate the optimal gp and api that maximizes the likelihood.
    private void estimatePeriodRepeatInterval(final int period, final IntRange repeatRange,
                                              final DragstrParamsBuilder builder, final EstimateDragstrParameters.StratifiedDragstrLocusCases cases) {
        int min_i = -1;
        int min_j = -1;
        final int minRepeat = repeatRange.getMinimumInteger();
        final int maxRepeat = repeatRange.getMaximumInteger();
        final double maxLog10PHet = Math.log10(het_hom_ratio) - Math.log10(1 + het_hom_ratio);
        double min_sum = Double.NEGATIVE_INFINITY;
// gp_min = -10*log10((1 - 0.5.^(1/(params.L_max*period))) / 2);
        for (int i = 0; i < phred_api_range.length; i++) {
            final double log10_p_het = Math.min(log10_api_range[i], maxLog10PHet);
            final double log10_p_hom = log10_p_het - log10_het_hom_ratio;
            final double log10_p_ref = Math.log10(1 - Math.pow(10, log10_p_het) - Math.pow(10, log10_p_hom));
            gp_for: for (int j = min_gp_idx_by_period[period - 1]; j < phred_gp_range.length; j++) {
                double sum = 0;
                for (int r = minRepeat; r <= maxRepeat; r++) {
                    final EstimateDragstrParameters.DragstrLocusCases repeatCases = cases.get(period, r);
                    final int size = repeatCases.size();
                    final int[] n = repeatCases.depth.toIntArray(), k = repeatCases.nonref.toIntArray();
                    final double log10_p_err = log10_p_err_by_len[j][period - 1][r - 1];
                    final double log10_p_no_err = log10_p_no_err_by_len[j][period - 1][r - 1];
                    for (int s = 0; s < size; s++) {
                        if ((sum += log10ProbFunc(n[s], k[s], log10_p_err, log10_p_no_err, log10_p_ref, log10_p_het, log10_p_hom)) < min_sum) {
                            continue gp_for;
                        };
                    }
                }
                if (sum > min_sum) {
                    min_i = i;
                    min_j = j;
                    min_sum = sum;
                }
            }
        }
        builder.set(period, repeatRange, phred_gp_range[min_j], 10.0 / period, phred_api_range[min_i]);

    }




    //%% prob_func runs into numerical problems for n > ~1000, so we must do the calculation
    //%% in the log domain.
    //log10sum = @(log10a, log10b) max(max(log10a, log10b), log10(10.^(log10a) + 10.^(log10b)));
    //
    //log10prob_func = @(n, k, p_err, p_het, p_hom) ...
    //log10sum(k .* log10(p_err) + (n-k) .* log10(1 - p_err) + log10(1 - p_het - p_hom), ...
    //log10sum(n .* log10(0.5) + log10(p_het), log10(p_hom * (n==k))) );
    //
    // We use a more accurate way to calculate log10-sum-log10 and log10(1+p).
    // In pratice should do very little to change things.
    private double log10ProbFunc(final int n, final int k, final double log10_p_err, final double log10_p_no_err, final double log10_p_ref, final double log10_p_het, final double log10_p_hom) {
        return MathUtils.log10SumLog10(
                log10_p_ref + k * log10_p_err + (n - k) * log10_p_no_err,
                log10_p_het + n * MathUtils.LOG10_ONE_HALF,
                   n == k ? log10_p_hom : Double.NEGATIVE_INFINITY); // this last hom_ref lk part is not accurate as is possible that an error that reversts to the reference does occur.
    }
}
