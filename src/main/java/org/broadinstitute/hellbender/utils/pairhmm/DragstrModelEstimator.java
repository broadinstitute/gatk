package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.lang.math.IntRange;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;


public class DragstrModelEstimator {

    private static double[][] DEFAULT_GOP = {
            {45.00, 45.00, 45.00, 45.00, 45.00, 45.00, 40.50, 33.50, 28.00, 24.00, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75},
            {39.50, 39.50, 39.50, 39.50, 36.00, 30.00, 27.25, 25.00, 24.25, 24.75, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.75},
            {38.50, 41.00, 41.00, 41.00, 41.00, 37.50, 35.25, 34.75, 34.75, 33.25, 33.25, 33.25, 32.50, 30.75, 28.50, 29.00, 29.00, 29.00, 29.00, 29.00},
            {37.50, 39.00, 39.00, 37.75, 34.00, 34.00, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 31.75, 31.75, 31.75, 31.75, 31.75},
            {37.00, 40.00, 40.00, 40.00, 36.00, 35.00, 24.50, 24.50, 24.50, 24.50, 22.50, 22.50, 22.50, 23.50, 23.50, 23.50, 23.50, 23.50, 23.50, 23.50},
            {36.25, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00},
            {36.00, 40.50, 40.50, 40.50, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75},
            {36.25, 39.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75}};

    private static double[][] DEFAULT_API = {
            {39.00, 39.00, 37.00, 35.00, 32.00, 26.00, 20.00, 16.00, 12.00, 10.00, 8.00, 7.00, 7.00, 6.00, 6.00, 5.00, 5.00, 4.00, 4.00, 4.00},
            {30.00, 30.00, 29.00, 22.00, 17.00, 14.00, 11.00, 8.00, 6.00, 5.00, 4.00, 4.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 2.00, 2.00},
            {27.00, 27.00, 25.00, 18.00, 14.00, 12.00, 9.00, 7.00, 5.00, 4.00, 3.00, 3.00, 3.00, 3.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00},
            {27.00, 27.00, 18.00, 9.00, 9.00, 9.00, 9.00, 3.00, 3.00, 3.00, 3.00, 3.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00},
            {29.00, 29.00, 18.00, 8.00, 8.00, 8.00, 4.00, 3.00, 3.00, 3.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00},
            {25.00, 25.00, 10.00, 10.00, 10.00, 4.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00},
            {21.00, 21.00, 11.00, 11.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00},
            {18.00, 18.00, 10.00, 6.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00}};

    private static final Logger logger = LogManager.getLogger(DragstrModelEstimator.class);

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
    private final boolean dont_adjust_gop;
    private final int minimum_depth = 10;

    public DragstrModelEstimator(final DragstrModelEstimatorArgumentCollection argumentCollection) {
        phred_gp_range = argumentCollection.phredGpValues.toDoubleArray();
        phred_api_range = argumentCollection.phredApiValues.toDoubleArray();
        log10_gp_range = MathUtils.applyToArray(phred_gp_range, d -> -.1 * d);
        log10_api_range = MathUtils.applyToArray(phred_api_range, d -> -0.1 * d);
        het_hom_ratio = argumentCollection.hetToHomRatio;
        log10_het_hom_ratio = Math.log10(het_hom_ratio);
        min_loci_count = argumentCollection.minLociCount;
        api_mono_threshold = argumentCollection.apiMonothresh;
        min_gop = argumentCollection.minGOP;
        max_gop = argumentCollection.maxGOP;
        log10_p_err_by_len = new double[phred_gp_range.length][8][20];
        log10_p_no_err_by_len = new double[phred_gp_range.length][8][20];
        min_gp_idx_by_period = new int[8];
        for (int i = 0; i < phred_gp_range.length; i++) {
            for (int k = 0; k < 8; k++) {
                final int period = k + 1;
                final double log10_p_no_err_per_pos = MathUtils.log10OneMinusPow10(LOG_10_OF_2 + log10_gp_range[i]);
                for (int j = 0; j < 20; j++) {
                    final int repeats = j + 1;
                    final int lengthInBases = repeats * period;
                    log10_p_no_err_by_len[i][k][j] = lengthInBases * log10_p_no_err_per_pos;
                    log10_p_err_by_len[i][k][j] = MathUtils.log10OneMinusPow10(log10_p_no_err_by_len[i][k][j]);
                }
            }
        }
        for (int i = 0; i < min_gp_idx_by_period.length; i++) {
            final int period = i + 1;
            final double gp_min = Math.ceil(-10 * Math.log10((1 - Math.pow(0.5, (1.0/ (20.0 * period)) / 2.0))));
            final int index = Arrays.binarySearch(phred_gp_range, gp_min);
            // since we are looking for a double, we have to be a bit tolerant in terms of differences,
            // so if no exact was found we look the position before the insertion position in case that value
            // is close enough (less than 0.001 way).
            min_gp_idx_by_period[i] = index >= 0 ? index :
                    (index < -1 && Math.abs(gp_min - phred_api_range[-index - 2]) < 0.001) ? -index - 2 : -index - 1;
        }
        dont_adjust_gop = argumentCollection.dontPostAdjustmentOfGOP;
    }

    public Estimate createEstimate(final int maxPeriod, final int maxRepeats) {
        return new Estimate(maxPeriod, maxRepeats);
    }

    public PeriodCases createPeriodCases(final int period, final int maxRepeats, final int totalCases) {
        return new PeriodCases(period, maxRepeats, totalCases);
    }

    public static class PeriodCases {
        public final int period;
        public final RepeatCases[] repeatCases;

        public PeriodCases(final int period, final int maxRepeats, final int casesPerRepeatCapacity) {
            this.period = period;
            this.repeatCases = new RepeatCases[maxRepeats];
            for (int i = 0; i < maxRepeats; i++) {
                this.repeatCases[i] = new RepeatCases(period,i + 1, casesPerRepeatCapacity);
            }
        }

        /**
         * Returns the {@link RepeatCases} instance that holds the cases for this period an a particular number of
         * repeat units. Notice that repeats counts over the maximum are collapsed into that maximum.
         * @param repeats
         * @return never {@code null}.
         * @throws IllegalArgumentException if {@code repeats} is 0 or less.
         */
        public RepeatCases getRepeatCases(final int repeats) {
            if (repeats > repeatCases.length) {
                return repeatCases[repeatCases.length - 1];
            } else {
                return repeatCases[repeats - 1];
            }
        }

        /**
         * The maximum number of repeats str.
         * @return 1 or greater.
         */
        public int getMaxRepeats() {
            return repeatCases.length;
        }

        public void removeLowDepth(final int minimum_depth) {
            for (final RepeatCases rc : repeatCases) {
                rc.removeLowDepth(minimum_depth);
            }
        }
    }

    public static class RepeatCases {
        final int period;
        final int repeat; // number of repeat units in this case.
        int[] n; // total number of reads/fragments considered.
        int[] k; // total number of reads/fragments that have a non-zero total indel sum within the str.
        int size;


        public int removeLowDepth(final int minDepth) {
            int i, newSize;
            for (i = 0, newSize = 0; i < size; newSize++) {
                if (n[i++] < minDepth) {
                    break;
                }
            }
            for (; i < size; i++) {
                if (n[i] >= minDepth) {
                   n[newSize] = n[i];
                   k[newSize++] = k[i];
                }
            }
            return size = newSize;
        }

        public int effectiveSampleSize(final int minDepth) {
            int i, result;
            for (i = 0, result = 0; i < size; result++) {
                if (n[i++] < minDepth) {
                    break;
                }
            }
            for (; i < size; i++) {
                if (n[i] >= minDepth) {
                    result++;
                }
            }
            return result;
        }

        public RepeatCases(final int period, final int repeat, final int capacity) {
            this.period = period;
            this.repeat = repeat;
            final int initialCapacity = capacity < 10 ? 10 : capacity;
            n = new int[initialCapacity];
            k = new int[initialCapacity];
            size = 0;
        }

        public int getPeriod() {
            return period;
        }

        public int getRepeats() {
            return repeat;
        }

        public int size() {
            return size;
        }

        public void add(final int n, final int k) {
            if (this.n.length == size) {
                this.n = Arrays.copyOf(this.n, size << 1);
                this.k = Arrays.copyOf(this.k, this.n.length);
            }
            this.n[size] = n;
            this.k[size++] = k;
        }
    }

    public class Estimate {
        public final double[][] gp;
        public final double[][] ph_het_variant;
        private final double gop[][];
        private final double[] gcpByPeriod;

        private Estimate(final int maxPeriod, final int maxRepeats) {
            this.gp = new double[maxPeriod][maxRepeats];
            this.ph_het_variant = new double[maxPeriod][maxRepeats];
            for (int i = 0; i < maxPeriod; i++) {
                final double[] defaultAPIValues = DEFAULT_API.length > i ? DEFAULT_API[i] : DEFAULT_API[DEFAULT_API.length - 1];
                final double[] defaultGOPValues = DEFAULT_GOP.length > i ? DEFAULT_GOP[i] : DEFAULT_GOP[DEFAULT_GOP.length - 1];
                for (int j = 0; j < maxRepeats; j++) {
                    ph_het_variant[i][j] = defaultAPIValues.length > i ? defaultAPIValues[j] : defaultAPIValues[defaultAPIValues.length - 1];
                    gp[i][j] = defaultAPIValues.length > i ? defaultGOPValues[j] : defaultGOPValues[defaultGOPValues.length - 1];
                }
            }
            this.gop = gp.clone();
            this.gcpByPeriod = new double[maxPeriod];
            for (int i = 0; i < maxPeriod; i++) {
                gcpByPeriod[i] = 10.0 / (i + 1);
            }
        }

        void set(final int period, final int repeat, final double gp, final double ph_het_variant) {
            this.gp[period - 1][repeat - 1] = gp;
            this.gop[period - 1] = null; // make sure it is recalculated.
            this.ph_het_variant[period - 1][repeat - 1] = ph_het_variant;
        }

        void set(final int period, final int minRepeat, final int maxRepeat, final double gp, final double ph_het_variant) {
            for (int repeat = minRepeat; repeat <= maxRepeat; repeat++) {
                set(period, repeat, gp, ph_het_variant);
            }
        }

        public double gcp(final int period, @SuppressWarnings("unused") final int repeat) {
            return gcpByPeriod[period - 1];
        }

        public double api(final int period, final int repeat) {
            return ph_het_variant[period - 1][repeat >= ph_het_variant[0].length ? ph_het_variant[0].length - 1 : repeat - 1];
        }

        public double gop(final int period, final int repeat) {
            double[] periodGop = gop[period - 1];
            if (periodGop == null) {
                gop[period -1] = periodGop = new double[gp[period - 1].length];
                for (int i = 0; i < periodGop.length; i++) {
                    periodGop[i] = Math.max(min_gop, gopCalculation(gp[period - 1][i], gcpByPeriod[period - 1], period, .25));
                }
            }
            return periodGop[repeat >= periodGop.length ? periodGop.length - 1 : repeat - 1];
        }
    }

    private double gopCalculation(final double gp, final double gcp, final int period, final double stepwise) {

        double gpProb = Math.pow(10, -.1 * gp);
        double bestGop = -1;
        final double c = Math.pow(10, -0.1 * gcp);
        double bestCost = Double.POSITIVE_INFINITY;
        for (double gop = 0.0; max_gop - gop > -0.001; gop += stepwise) {
            final double g = Math.pow(10, -0.1 * gop);
            final double p_gap = g * Math.pow(c, period - 1) * (1.0 - c);
            final double p_no_gap = Math.pow(1 - 2 * g, period + 1);
            final double cost = Math.abs(p_gap / p_no_gap - gpProb);
            if (cost < bestCost) {
                bestGop = gop;
                bestCost = cost;
            }
        }
        return bestGop;
    }

    public void estimate(final Estimate result, final PeriodCases cases) {
        final int period = cases.period;
        // We calculate the "flanks" repeat counts that have too little data
        // to be analyzed separately:
        //cases.removeLowDepth(minimum_depth);
        int accum = 0;
        int leftFlank;
        for (leftFlank = 0; leftFlank < cases.repeatCases.length; leftFlank++) {
            if ((accum += cases.repeatCases[leftFlank].effectiveSampleSize(minimum_depth)) > min_loci_count) {
                break;
            }
        }
        int rightFlank; accum = 0;
        for (rightFlank = cases.repeatCases.length - 1; rightFlank >= 0; rightFlank--) {
            if ((accum += cases.repeatCases[rightFlank].effectiveSampleSize(minimum_depth)) > min_loci_count) {
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
        final ArrayDeque<IntRange> pending = new ArrayDeque<>(cases.getMaxRepeats());
        final ArrayDeque<IntRange> done = new ArrayDeque<>(cases.getMaxRepeats());

        pending.add(new IntRange(1, leftFlank + 1));
        for (int i = leftFlank + 1; i < rightFlank; i++) {
            pending.add(new IntRange(i + 1));
        }
        pending.add(new IntRange(rightFlank + 1, cases.getMaxRepeats()));
        IntRange last = null;
        do {
            final IntRange next = pending.pop();
            if (next.getMinimumInteger() != next.getMaximumInteger()) {
                logger.debug("calculating parameters for cases with repeat count " + next.getMinimumInteger() + " to " + next.getMaximumNumber());
            } else {
                logger.debug("calculating parameters for cases with repeat count " + next.getMinimumInteger());
            }
            estimate(result, period, next.getMinimumInteger(), next.getMaximumInteger(), cases);
            final double gp1 = result.gp[period - 1][next.getMinimumInteger() - 1];
            final double api1 = result.ph_het_variant[period - 1][next.getMinimumInteger() - 1];
            if (done.isEmpty() || (result.gp[period - 1][last.getMaximumInteger() - 1] >= gp1 &&
                                   result.ph_het_variant[period - 1][last.getMaximumInteger() - 1] + api_mono_threshold >= api1)) {
                done.addLast(last = next);
            } else {
                pending.push(new IntRange(last.getMinimumNumber(), next.getMaximumNumber()));
                done.removeLast();
                last = !done.isEmpty() ? done.getLast() : null;
            }
        } while (!pending.isEmpty());
    }

    // This seem to accom
    private void smoothParamValuesOverNumberOfRepeats(final Estimate result, final int period) {
        final int index = period - 1;
        final double[] gp = result.gp[index];
        final double[] api = result.ph_het_variant[index];
        double minApi = api[0];
        for (int i = 1; i < gp.length; i++) {
            if (gp[i] > gp[i - 1] || api[i] > minApi + api_mono_threshold) {
                gp[i] = gp[i - 1];
                api[i] = api[i - 1];
            } else if (api[i] < minApi) {
                minApi = api[i];
            }
        }
    }


    // Given a observed het/hom ratio the total number of reads, the variable ones and the length of the STR in bases.
    // calculate the optimal gp and api that maximizes the likelihood.
    private void estimate(final Estimate result, final int period, final int minRepeat, final int maxRepeat, final PeriodCases cases) {
        //final double[][] y = new double[phred_api_range.length][phred_gp_range.length];
        int min_i = -1;
        int min_j = -1;
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
                    final RepeatCases repeatCases = cases.getRepeatCases(r);
                    final int size = repeatCases.size;
                    final int[] n = repeatCases.n, k = repeatCases.k;
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
        result.set(period, minRepeat, maxRepeat, phred_gp_range[min_j], phred_api_range[min_i]);
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
