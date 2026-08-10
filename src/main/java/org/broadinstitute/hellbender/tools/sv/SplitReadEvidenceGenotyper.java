package org.broadinstitute.hellbender.tools.sv;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.aggregation.EvidenceStatUtils;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.*;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class SplitReadEvidenceGenotyper {

    private final Integer minSize;
    private final int trainingCountCutoff;
    private final double targetCoverage;
    private final int maxQuality;
    private final Map<String, Double> sampleCoverageMap;
    private final Map<String, FirstPassResult> firstPassCounts = new HashMap<>();
    private Double hetMedian = null;
    private Double hetMad = null;
    private Double hetCutoff = null;
    private final List<Double> hetCounts = new ArrayList<>();
    private final List<Double> homCounts = new ArrayList<>();
    private boolean firstPassMade = false;
    private boolean secondPassMade = false;
    private boolean thirdPassMade = false;

    private final int rareMin;
    private final int rareMax;
    private final int commonMin;
    private final int commonMax;

    // Recovery histograms for SR frequency cutoff optimization, matching v1.1's union semantics.
    // In v1.1 (SR_genotype.opt_part1.sh + optimalsrcutoff.sh), ALL non-ref entries go into
    // recover.single.txt, and variants with bothside support ALSO go into recover.bothsides.txt.
    // The optimalsrcutoff.sh script combines them with `sort -u` on VID@Sample (set union).
    // To replicate this efficiently, we use three histograms and inclusion-exclusion:
    //   union(i,j) = single(i) + both(j) - overlap(i,j)
    //
    // singleHistogram: ALL non-ref sample entries, indexed by single-sided frac
    //   Dimensions: [freq (rare=0, common=1)][pass=0, fail=1][single_frac_bin 0..10]
    // bothHistogram: non-ref sample entries from bothside variants, indexed by both-sided frac
    //   Dimensions: [freq][pass][both_frac_bin 0..10]
    // overlapHistogram: non-ref sample entries from bothside variants, indexed by BOTH fracs
    //   Dimensions: [freq][pass][single_frac_bin 0..10][both_frac_bin 0..10]
    private static final int NUM_FRAC_BINS = 11;
    private final long[][][] singleHistogram = new long[2][2][NUM_FRAC_BINS];
    private final long[][][] bothHistogram = new long[2][2][NUM_FRAC_BINS];
    private final long[][][][] overlapHistogram = new long[2][2][NUM_FRAC_BINS][NUM_FRAC_BINS];

    private static final Median MEDIAN = new Median();
    private static final NormalDistribution Z_DISTRIBUTION = new NormalDistribution();

    private static final Set<Integer> HET_COPY_STATES = Set.of(1, 3);

    // ------------------------------------------------------------------------
    // Diagnostics. Purely observational: nothing below participates in
    // genotyping or cutoff selection. Rendered by cutoffDiagnosticsReport().
    //
    // Motivation: when the frac histograms concentrate in a single bin the
    // 11x11 cutoff grid becomes non-discriminative and every cell scores alike.
    // validateCutoffGrid() now rejects that outright, but knowing WHY it
    // collapsed still requires the histograms and the scored grid, so we
    // record them and report them alongside the rejection.
    // ------------------------------------------------------------------------

    /** Buckets for the integer sample-count inputs to the two frac metrics. */
    private static final String[] COUNT_BIN_LABELS = {
            "0", "1", "2", "3", "4", "5", "6-7", "8-15", "16-31", "32-63",
            "64-127", "128-255", "256-511", "512+"
    };
    private static final int NUM_COUNT_BINS = COUNT_BIN_LABELS.length;

    // Variant-level filter tallies from the histogram accumulation passes
    private long diagVariantsSeen = 0;
    private long diagVariantsNotCnv = 0;
    private long diagVariantsTooSmall = 0;
    private long diagVariantsNoSrEvidence = 0;
    private long diagVariantsNoNonRef = 0;
    private long diagVariantsPassFlag = 0;
    private long diagVariantsNoSamplesOverOne = 0;
    private long diagVariantsNoBothsideNonZero = 0;
    private long diagVariantsNoTwoSidedPass = 0;
    private long diagVariantsContributedSingle = 0;
    private long diagVariantsContributedBoth = 0;
    private long diagVariantsOutsideBothFreqBins = 0;

    // Frac saturation detail. The 11-bin histograms cannot distinguish a frac of
    // exactly 1.0 from one exceeding 1.0, but the two have different causes: the
    // numerator threshold (sr_count) being at or below the denominator threshold
    // permits frac > 1. Tracked at variant level.
    private long diagSingleFracGtOne = 0;
    private long diagSingleFracEqOne = 0;
    private long diagSingleFracNearOne = 0; // [0.95, 1.0)
    private long diagBothFracGtOne = 0;
    private long diagBothFracEqOne = 0;
    private long diagBothFracNearOne = 0;
    private double diagSingleFracSum = 0;
    private long diagSingleFracN = 0;
    private double diagBothFracSum = 0;
    private long diagBothFracN = 0;

    // Distributions of the raw counts feeding the fracs
    private final long[] diagNonRefCountHist = new long[NUM_COUNT_BINS];
    private final long[] diagSamplesOverOneHist = new long[NUM_COUNT_BINS];
    private final long[] diagTwoSidedPassHist = new long[NUM_COUNT_BINS];
    private final long[] diagBothsideNonZeroHist = new long[NUM_COUNT_BINS];

    // Training-pass statistics, captured before the count lists are released
    private long diagFirstPassVariants = 0;
    private long diagFirstPassHetN = 0;
    private long diagSecondPassHetN = 0;
    private long diagSecondPassHomN = 0;
    private Double diagSecondPassHetMedian = null;
    private Double diagSecondPassHomMedian = null;
    private Double diagSecondPassHetMad = null;

    // Scored cutoff grids and their selection outcomes, retained from finalizeThirdPass
    private CutoffGrid diagRareGrid = null;
    private CutoffGrid diagCommonGrid = null;
    private SelectionOutcome diagRareOutcome = null;
    private SelectionOutcome diagCommonOutcome = null;

    /** {@link CutoffGrid#selectedIndex()} when the objective could not distinguish a cell. */
    public static final int NO_CUTOFF_SELECTED = -1;

    /**
     * A fully scored cutoff grid: every (fracSingle, fracBoth) cell, its score, the
     * baseline it was scored against, and the selected index, which is
     * {@link #NO_CUTOFF_SELECTED} if no cell scored a finite value.
     */
    public record CutoffGrid(List<CutoffResult> cells, double[] scores, double baseline,
                             int selectedIndex, int freqMin, int freqMax) {}

    public SplitReadEvidenceGenotyper(final Map<String,Double> sampleCoverageMap, final int numSamples, final double qualityCutoff, final Integer minSize, final double targetCoverage, final int maxQuality) {
        this.sampleCoverageMap = Utils.nonNull(sampleCoverageMap);
        this.trainingCountCutoff = computeCountCutoff(qualityCutoff);
        Utils.validateArg(maxQuality > 0, "Maximum quality must be greater than zero");
        this.maxQuality = maxQuality;
        this.minSize = minSize;
        this.targetCoverage = targetCoverage;
        this.rareMin = 0;
        this.rareMax = Math.max(numSamples / 100, 2);
        this.commonMin = rareMax;
        this.commonMax = numSamples;
    }

    private static int computeCountCutoff(final double qualityCutoff) {
        int i = 1;
        while (true) {
            final PoissonDistribution dist = new PoissonDistribution(i);
            final double p = dist.cumulativeProbability(0);
            if (p == 0) {
                throw new IllegalArgumentException("Precision error - quality cutoff " + qualityCutoff + " is too high");
            }
            final double qual = QualityUtils.errorProbToQual(p);
            if (qual > qualityCutoff) {
                return Math.max(i - 1, 1);
            }
            i++;
        }
    }

    public void addFirstPass(final SVCallRecord record, final List<SplitReadEvidence> startEvidence,
                             final List<SplitReadEvidence> endEvidence, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype, final List<String> depthGenotypeSamples) {
        Utils.nonNull(startEvidence);
        Utils.nonNull(endEvidence);
        Utils.validate(!firstPassMade, "First pass already made");
        // TODO check all samples are in coverage map
        final int minCount = Math.max(trainingCountCutoff / 2, 1);
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final Set<String> bothSideSupportSamples = hasBothSideSupport(startCounts, endCounts, minCount);
        if (!bothSideSupportSamples.isEmpty()) {
            firstPassCounts.put(record.getId(), new FirstPassResult(bothSideSupportSamples, startCounts, endCounts, depthGenotype, depthGenotypeSamples));
        }
    }

    // Assumes samples with 0 counts are not present in the key set
    // Uses strict > to match v1.1's awk '$NF>(sr_count/2)' and stay consistent with
    // countBothSideSupport(), which also uses >. Both methods port the same v1.1 pattern.
    // Note: v1.1 uses awk float division (sr_count/2=3.5 for sr_count=7) while v1.2 uses
    // Java integer division (trainingCountCutoff/2=3 for 7). The threshold value differs
    // but the comparison operator must be consistent across hasBothSideSupport and
    // countBothSideSupport to avoid a sample being counted as two-sided in Phase 1 but
    // not in Phase 2/3.
    private static Set<String> hasBothSideSupport(final Map<String, Double> startCounts,
                                              final Map<String, Double> endCounts, final double threshold) {
        final Set<String> result = new HashSet<>();
        for (final Map.Entry<String, Double> entry : startCounts.entrySet()) {
            if (entry.getValue() > threshold) {
                if (endCounts.getOrDefault(entry.getKey(), 0.0) > threshold) {
                    result.add(entry.getKey());
                }
            }
        }
        return result;
    }

    // Assumes samples with 0 counts are not present in the key set
    private static int countBothSideSupport(final Map<String, Double> startCounts,
                                            final Map<String, Double> endCounts, final double threshold) {
        int count = 0;
        for (final Map.Entry<String, Double> entry : startCounts.entrySet()) {
            if (entry.getValue() > threshold) {
                if (endCounts.getOrDefault(entry.getKey(), 0.0) > threshold) {
                    count++;
                }
            }
        }
        return count;
    }

    /**
     * Count samples whose normalized summed split read evidence exceeds the threshold.
     * This matches the legacy SR training path, which learns single-sided background ratios
     * from {@code SR_sum} after coverage normalization rather than raw start/end counts.
     */
    private static int countNormalizedSummedSupport(final Map<String, Double> startCounts,
                                                    final Map<String, Double> endCounts,
                                                    final double threshold) {
        final Set<String> samples = new HashSet<>(startCounts.keySet());
        samples.addAll(endCounts.keySet());
        int count = 0;
        for (final String sample : samples) {
            final double sum = startCounts.getOrDefault(sample, 0.0) + endCounts.getOrDefault(sample, 0.0);
            if (sum > threshold) {
                count++;
            }
        }
        return count;
    }

    private Map<String, Double> normalizeCounts(final List<SplitReadEvidence> evidence) {
        final Map<String, Double> counts = new HashMap<>();
        for (final SplitReadEvidence e : evidence) {
            final double sampleCoverage = sampleCoverageMap.getOrDefault(e.getSample(), 0.);
            final double normalizedCount = EvidenceStatUtils.getNormalizedCount(e.getCount(), sampleCoverage, targetCoverage);
            if (normalizedCount > 0) {
                counts.put(e.getSample(), normalizedCount);
            }
        }
        return counts;
    }

    public void finalizeFirstPass() {
        Utils.validate(!firstPassCounts.isEmpty(), "No split read counts after first pass");
        Utils.validate(!firstPassMade, "First pass has already been made");
        final double[] hetCounts = firstPassCounts.values().stream().map(c -> c.getCounts(HET_COPY_STATES)).flatMapToDouble(DoubleStream::of).toArray();
        hetMedian = MEDIAN.evaluate(hetCounts);
        final double[] deviations = DoubleStream.of(hetCounts).map(d -> Math.abs(d - hetMedian)).toArray();
        hetMad = MEDIAN.evaluate(deviations);
        hetCutoff = hetMedian + 1.645 * hetMad;
        diagFirstPassVariants = firstPassCounts.size();
        diagFirstPassHetN = hetCounts.length;
        firstPassMade = true;
    }

    public record CutoffResult(double fracSingle, double fracBoth, long countPass, long countFail, int freqMin, int freqMax) {}

    /**
     * Compute suffix sums from a histogram: result[k] = number of entries with frac >= cutoffs[k].
     * The histogram bin j holds the count of entries with frac in [j*0.1, (j+1)*0.1); bin 10 holds frac == 1.0.
     */
    private static long[] suffixSumsFromHistogram(final long[] histogram) {
        final long[] suffixSums = new long[NUM_FRAC_BINS];
        suffixSums[NUM_FRAC_BINS - 1] = histogram[NUM_FRAC_BINS - 1];
        for (int k = NUM_FRAC_BINS - 2; k >= 0; k--) {
            suffixSums[k] = suffixSums[k + 1] + histogram[k];
        }
        return suffixSums;
    }

    /**
     * Map a frac value in [0, 1] to a histogram bin index in [0, 10].
     * Bin k covers [k*0.1, (k+1)*0.1); bin 10 covers exactly 1.0.
     */
    private static int fracToBin(final double frac) {
        final int bin = (int) (frac * 10);
        return Math.min(bin, NUM_FRAC_BINS - 1);
    }

    /**
     * Increment the single-sided recovery histogram. Called for ALL non-ref samples.
     * @param count the variant-level count (nonRefCount) used for frequency bin selection
     * @param frac the single-sided ratio (nonRefCount / samplesOverOneCount)
     * @param pass whether this entry passes training criteria
     */
    private void addToSingleHistogram(final int count, final double frac, final boolean pass) {
        final int fracBin = fracToBin(frac);
        final int passIndex = pass ? 0 : 1;
        if (count > rareMin && count <= rareMax) {
            singleHistogram[0][passIndex][fracBin]++;
        }
        if (count > commonMin && count <= commonMax) {
            singleHistogram[1][passIndex][fracBin]++;
        }
    }

    /**
     * Increment the both-sided recovery histogram. Called for non-ref samples of variants
     * with bothside support.
     * @param count the variant-level count (twoSidedPassCount) used for frequency bin selection
     * @param frac the both-sided ratio (twoSidedPassCount / bothsideNonZeroCount)
     * @param pass whether this entry passes training criteria
     */
    private void addToBothHistogram(final int count, final double frac, final boolean pass) {
        final int fracBin = fracToBin(frac);
        final int passIndex = pass ? 0 : 1;
        if (count > rareMin && count <= rareMax) {
            bothHistogram[0][passIndex][fracBin]++;
        }
        if (count > commonMin && count <= commonMax) {
            bothHistogram[1][passIndex][fracBin]++;
        }
    }

    /**
     * Increment the overlap histogram (entries that appear in BOTH single and both histograms).
     * Used for inclusion-exclusion: union(i,j) = single(i) + both(j) - overlap(i,j).
     * An entry only appears in the overlap for a given frequency bin if it passes
     * BOTH the single-sided and both-sided frequency range tests for that bin.
     */
    private void addToOverlapHistogram(final int singleCount, final double singleFrac,
                                       final int bothCount, final double bothFrac,
                                       final boolean pass) {
        final int singleFracBin = fracToBin(singleFrac);
        final int bothFracBin = fracToBin(bothFrac);
        final int passIndex = pass ? 0 : 1;
        // Rare bin: must pass both single (nonRefCount) and both (twoSidedPassCount) freq checks
        if (singleCount > rareMin && singleCount <= rareMax
                && bothCount > rareMin && bothCount <= rareMax) {
            overlapHistogram[0][passIndex][singleFracBin][bothFracBin]++;
        }
        // Common bin: must pass both freq checks
        if (singleCount > commonMin && singleCount <= commonMax
                && bothCount > commonMin && bothCount <= commonMax) {
            overlapHistogram[1][passIndex][singleFracBin][bothFracBin]++;
        }
    }

    public SplitReadGenotypeFrequencyCutoffs finalizeThirdPass() {
        Utils.validate(!thirdPassMade, "Third pass has already been made");
        final double[] cutoffs = IntStream.rangeClosed(0, 10).mapToDouble(i -> i * 0.1).toArray();
        diagRareGrid = cutoffOptimizationFromHistograms(cutoffs, 0);
        diagCommonGrid = cutoffOptimizationFromHistograms(cutoffs, 1);
        diagRareOutcome = classifyCutoffGrid(diagRareGrid, "rare");
        diagCommonOutcome = classifyCutoffGrid(diagCommonGrid, "common");
        thirdPassMade = true;
        // A rejected grid falls back to the (0.0, 0.0) cell, matching prior behavior. The
        // rejection is recorded in the diagnostics report and surfaced by
        // cutoffSelectionOutcomes() so the caller can log it and a downstream task can act on
        // it; see classifyCutoffGrid for why this does not throw.
        final CutoffResult rareResult = diagRareGrid.cells().get(Math.max(diagRareGrid.selectedIndex(), 0));
        final CutoffResult commonResult = diagCommonGrid.cells().get(Math.max(diagCommonGrid.selectedIndex(), 0));
        return new SplitReadGenotypeFrequencyCutoffs(rareResult, commonResult);
    }

    /**
     * Cutoff optimization using inclusion-exclusion on histogram suffix sums to replicate
     * v1.1's {@code sort -u} union semantics across single-sided and both-sided entries.
     *
     * <p>In v1.1, a variant with bothside support appears in BOTH recover.single.txt and
     * recover.both.txt. The optimalsrcutoff.sh script concatenates them and deduplicates
     * with {@code sort -u} on VID@Sample. To replicate this efficiently:</p>
     * <pre>
     *   union_pass(i, j) = singlePass(i) + bothPass(j) - overlapPass(i, j)
     * </pre>
     * where {@code overlapPass(i, j)} is the count of entries that pass BOTH thresholds.
     *
     * <p>Returns the whole scored grid rather than only the winning cell so that the
     * selection can be audited; {@link CutoffGrid#selectedIndex()} is the argmax that
     * callers should use, leaving selection behavior unchanged.</p>
     */
    private CutoffGrid cutoffOptimizationFromHistograms(final double[] cutoffs, final int freqBinIndex) {
        final long[] singlePassSuffix = suffixSumsFromHistogram(singleHistogram[freqBinIndex][0]);
        final long[] singleFailSuffix = suffixSumsFromHistogram(singleHistogram[freqBinIndex][1]);
        final long[] bothPassSuffix = suffixSumsFromHistogram(bothHistogram[freqBinIndex][0]);
        final long[] bothFailSuffix = suffixSumsFromHistogram(bothHistogram[freqBinIndex][1]);

        // 2D suffix sums for overlap: overlapPassSuffix2D[s][b] = sum of overlap[s'][b'] for s'>=s AND b'>=b
        final long[][] overlapPassSuffix2D = suffix2DFromHistogram(overlapHistogram[freqBinIndex][0]);
        final long[][] overlapFailSuffix2D = suffix2DFromHistogram(overlapHistogram[freqBinIndex][1]);

        final int freqMin = freqBinIndex == 0 ? rareMin : commonMin;
        final int freqMax = freqBinIndex == 0 ? rareMax : commonMax;

        final List<CutoffResult> combine = new ArrayList<>(cutoffs.length * cutoffs.length);
        for (int s = 0; s < cutoffs.length; s++) {
            for (int b = 0; b < cutoffs.length; b++) {
                // Inclusion-exclusion: union = single + both - overlap
                final long passCount = singlePassSuffix[s] + bothPassSuffix[b] - overlapPassSuffix2D[s][b];
                final long failCount = singleFailSuffix[s] + bothFailSuffix[b] - overlapFailSuffix2D[s][b];
                combine.add(new CutoffResult(cutoffs[s], cutoffs[b], passCount, failCount, freqMin, freqMax));
            }
        }
        final double baseline = computeBaseline(combine);
        final double[] scores = combine.stream().mapToDouble(c -> computeCutoffScore(c, baseline)).toArray();
        return new CutoffGrid(combine, scores, baseline, selectCutoffIndex(scores), freqMin, freqMax);
    }

    /**
     * Compute 2D suffix sums from a 2D histogram: result[s][b] = sum of hist[s'][b'] for s'>=s AND b'>=b.
     */
    private static long[][] suffix2DFromHistogram(final long[][] histogram) {
        final long[][] suffix = new long[NUM_FRAC_BINS][NUM_FRAC_BINS];
        // Start from bottom-right corner and work backwards
        for (int s = NUM_FRAC_BINS - 1; s >= 0; s--) {
            for (int b = NUM_FRAC_BINS - 1; b >= 0; b--) {
                suffix[s][b] = histogram[s][b]
                        + (s + 1 < NUM_FRAC_BINS ? suffix[s + 1][b] : 0)
                        + (b + 1 < NUM_FRAC_BINS ? suffix[s][b + 1] : 0)
                        - (s + 1 < NUM_FRAC_BINS && b + 1 < NUM_FRAC_BINS ? suffix[s + 1][b + 1] : 0);
            }
        }
        return suffix;
    }

    private static double computeBaseline(final List<CutoffResult> list) {
        for (final CutoffResult result : list) {
            if (result.fracSingle == 0 && result.fracBoth == 0) {
                return result.countPass;
            }
        }
        throw new IllegalArgumentException("List did not contain 0-fraction entry");
    }

    final double computeCutoffScore(final CutoffResult cutoffResult, final double baseline) {
        final double a = cutoffResult.countFail / (double) (cutoffResult.countFail + cutoffResult.countPass);
        final double b = (cutoffResult.countPass / baseline) - 1;
        return -((a * a) + (b * b));
    }

    /**
     * Pick the winning grid cell, or {@link #NO_CUTOFF_SELECTED} if the objective cannot
     * distinguish one.
     *
     * <p>Ties keep the first, and therefore lowest, cutoff pair, which is both the historical
     * behavior and the right convention: when every threshold from some point upward separates
     * the training data equally well, the lowest is the least extrapolation beyond what the
     * data shows. Raising it further would reject ratios that were never observed to be bad.</p>
     *
     * <p>Two differences from {@code MathUtils.maxElementIndex}, which this replaces. It skips
     * NaN scores rather than letting a NaN at index 0 suppress every later comparison, and it
     * reports {@link #NO_CUTOFF_SELECTED} instead of falling back to index 0 when nothing
     * scored finite. Index 0 is the (0.0, 0.0) cell, and zero cutoffs are not a harmless
     * default: at application time {@code rare.freqMax == common.freqMin}, so a zero cutoff
     * makes the frequency predicate a tautology and disables SR background filtering
     * outright. A grid with no usable signal has to be reported, not quietly resolved.</p>
     */
    private static int selectCutoffIndex(final double[] scores) {
        int best = NO_CUTOFF_SELECTED;
        for (int i = 0; i < scores.length; i++) {
            if (Double.isNaN(scores[i])) {
                continue;
            }
            if (best == NO_CUTOFF_SELECTED || scores[i] > scores[best]) {
                best = i;
            }
        }
        return best;
    }

    /** Machine-readable outcome of cutoff selection for one frequency bin. */
    public enum SelectionStatus {
        OK,
        /** Baseline pass count was zero: no positive training examples, objective undefined. */
        REJECTED_NO_PASSING_ENTRIES,
        /** Every cell scored identically: the recovery fractions carry no signal. */
        REJECTED_DEGENERATE_GRID
    }

    public record SelectionOutcome(SelectionStatus status, String detail) {
        public boolean rejected() {
            return status != SelectionStatus.OK;
        }
    }

    /**
     * Classify a cutoff grid, reporting whether it carries usable signal.
     *
     * <p>Deliberately returns an outcome rather than throwing. Cromwell only delocalizes task
     * outputs when the command succeeds, so failing here would strand the diagnostics report on
     * the worker, in exactly the case where it is needed. Detection therefore stays in this
     * tool and enforcement belongs downstream of it, where the report has already been copied
     * out. See the ValidateSRCutoffs task in GenotypeBatch.wdl.</p>
     *
     * <p>A rejected grid still resolves to the (0.0, 0.0) cell for the emitted parameter table,
     * preserving prior behavior, but that choice is now recorded rather than silent. Zero
     * cutoffs are not harmless: at application time {@code rare.freqMax == common.freqMin}, so a
     * zero cutoff makes the frequency predicate a tautology and disables SR background
     * filtering.</p>
     */
    private SelectionOutcome classifyCutoffGrid(final CutoffGrid grid, final String freqBinName) {
        final String context = String.format(
                "%s frequency bin (freq_min=%d, freq_max=%d): pass=%d fail=%d at zero cutoffs, "
                        + "sr_count=%d, two-sided threshold=%d. See the .sr_cutoff_diagnostics.txt "
                        + "output for the frac histograms and the full scored grid.",
                freqBinName, grid.freqMin(), grid.freqMax(),
                grid.cells().get(0).countPass(), grid.cells().get(0).countFail(),
                trainingCountCutoff, Math.max(trainingCountCutoff / 2, 1));

        if (grid.selectedIndex() == NO_CUTOFF_SELECTED) {
            return new SelectionOutcome(SelectionStatus.REJECTED_NO_PASSING_ENTRIES,
                    "No split read training entry satisfied the pass criteria, so every candidate cutoff "
                            + "scored undefined, for the " + context);
        }

        int tiedAtMax = 0;
        final double max = grid.scores()[grid.selectedIndex()];
        for (final double score : grid.scores()) {
            if (score == max) {
                tiedAtMax++;
            }
        }
        if (tiedAtMax == grid.scores().length) {
            return new SelectionOutcome(SelectionStatus.REJECTED_DEGENERATE_GRID,
                    "All " + grid.scores().length + " candidate cutoffs scored identically, so the split read "
                            + "recovery fractions carry no signal to optimize against. This happens when the "
                            + "fractions collapse into one histogram bin, most often because the SR quality "
                            + "cutoff (--sr-quality) yields an sr_count low enough that background samples are "
                            + "genotyped non-ref. Affects the " + context);
        }
        return new SelectionOutcome(SelectionStatus.OK, "");
    }

    public void addSecondPass(final SVCallRecord record, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype,
                              final List<String> depthGenotypeSamples) {
        Utils.nonNull(record);
        Utils.nonNull(depthGenotype);
        Utils.nonNull(depthGenotypeSamples);
        Utils.validate(firstPassMade, "First pass must be made");
        Utils.validate(!secondPassMade, "First pass has already been made");
        if (!firstPassCounts.containsKey(record.getId())) {
            // Skip if this variant wasn't added in the first pass
            return;
        }
        final FirstPassResult firstPassResult = firstPassCounts.get(record.getId());
        if (firstPassResult == null) {
            throw new IllegalArgumentException("Record " + record.getId() + " wasn't added in the first pass");
        }
        final int[] copyStates = depthGenotype.copyStates();
        Utils.validateArg(copyStates.length == depthGenotypeSamples.size(), "Copy states array and sample list must be the same size");
        for (int i = 0; i < copyStates.length; i++) {
            final String sample = depthGenotypeSamples.get(i);
            final Double count = firstPassResult.getCount(sample);
            if (count != null && !(count < hetCutoff && (copyStates[i] == 0 || copyStates[i] >= 4))) {
                if (copyStates[i] == 0 || copyStates[i] == 4) {
                    homCounts.add(count);
                } else if (copyStates[i] == 1 || copyStates[i] == 3) {
                    hetCounts.add(count);
                }
            }
        }
    }

    public SplitReadGenotypeParameters finalizeSecondPass() {
        Utils.validate(!homCounts.isEmpty(), "No split read counts after second pass");
        Utils.validate(!hetCounts.isEmpty(), "No split read counts after second pass");
        Utils.validate(!secondPassMade, "Second pass has already been made");
        final double[] homArr = homCounts.stream().mapToDouble(Double::doubleValue).toArray();
        final double[] hetArr = hetCounts.stream().mapToDouble(Double::doubleValue).toArray();
        final double homMedian = MEDIAN.evaluate(homArr);
        final double hetMedian = MEDIAN.evaluate(hetArr);
        final double hetMadValue = MEDIAN.evaluate(DoubleStream.of(hetArr).map(d -> Math.abs(d - hetMedian)).toArray());
        final double sdHet = 1.645 * 1.4826 * hetMadValue;
        diagSecondPassHetN = hetArr.length;
        diagSecondPassHomN = homArr.length;
        diagSecondPassHetMedian = hetMedian;
        diagSecondPassHomMedian = homMedian;
        diagSecondPassHetMad = hetMadValue;
        secondPassMade = true;
        // Free training accumulation data that is no longer needed
        firstPassCounts.clear();
        hetCounts.clear();
        homCounts.clear();
        return new SplitReadGenotypeParameters(trainingCountCutoff, homMedian, sdHet);
    }
    public SplitReadGenotypeResult genotype(final SVCallRecord record,
                                                    final List<SplitReadEvidence> startEvidence,
                                                    final List<SplitReadEvidence> endEvidence,
                                                    final SplitReadGenotypeMetrics metrics,
                                                    final int medianHomIns,
                                                    final double medianHomCutoffMultiplier,
                                                    final List<String> samples) {
        final SplitReadGenotypeParameters parameters = metrics.parameters;
        final SplitReadGenotypeFrequencyCutoffs frequencyCutoffs = metrics.cutoffs;
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final int[] genotypes = new int[samples.size()];
        final double[] countSum = new double[samples.size()];
        final int[] genotypeQuals = new int[samples.size()];
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype = record.getType();
        final double medianHom = svtype == GATKSVVCFConstants.StructuralVariantAnnotationType.INS ? medianHomIns : parameters.medianHom;
        final double bothsideMedianHet = 0.5 * medianHom;
        final double bothsideHomCutoff = medianHomCutoffMultiplier * bothsideMedianHet;
        final double bothsideSdHet = parameters.sdHet;
        final double onesideMedianHet = 0.5 * bothsideMedianHet;
        final double onesideHomCutoff = medianHomCutoffMultiplier * onesideMedianHet;
        final double onesideSdHet = 0.5 * parameters.sdHet;
        final double twoSidedThreshold = parameters.minCount / 2.;
        int nonRefCount = 0;
        final boolean[] twoSided = new boolean[samples.size()];
        final List<Double> bothsideNonRefCounts = new ArrayList<>();
        final List<Double> onesideNonRefCounts = new ArrayList<>();
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double startCount = startCounts.getOrDefault(sample, 0.);
            final double endCount = endCounts.getOrDefault(sample, 0.);
            countSum[i] = startCount + endCount;
            twoSided[i] = startCount > twoSidedThreshold && endCount > twoSidedThreshold;
            final double medianHet;
            final double homCutoff;
            if (twoSided[i]) {
                medianHet = bothsideMedianHet;
                homCutoff = bothsideHomCutoff;
            } else {
                medianHet = onesideMedianHet;
                homCutoff = onesideHomCutoff;
            }
            if (countSum[i] < parameters.minCount) {
                genotypes[i] = 0;
            } else if (countSum[i] <= homCutoff) {
                genotypes[i] = 1;
            } else {
                genotypes[i] = Math.max((int) ((countSum[i] / medianHet) + 0.5), 2);
            }
            if (genotypes[i] != 0) {
                ++nonRefCount;
                if (twoSided[i]) {
                    bothsideNonRefCounts.add(countSum[i]);
                } else {
                    onesideNonRefCounts.add(countSum[i]);
                }
            }
        }
        final boolean hasSplitReadEvidence = record.getEvidence().contains(GATKSVVCFConstants.EvidenceTypes.SR);
        final int twosidedMinCount = Math.max((int) parameters.minCount()/ 2, 1);
        final int twoSidedPassCount = countBothSideSupport(startCounts, endCounts, twosidedMinCount);
        final int bothsideNonZeroCount = countBothSideSupport(startCounts, endCounts, 0);
        final int samplesOverOneCount = countNormalizedSummedSupport(startCounts, endCounts, 1.0);
        boolean bothsidePass = false;
        boolean onesidePass = false;
        if (bothsideNonZeroCount > 0) {
            final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
            // v1.1 uses twoSidedPassCount ($2 in recover.bothsides.txt) for frequency binning
            bothsidePass = (backgroundRatio >= frequencyCutoffs.rare.fracBoth && twoSidedPassCount <= frequencyCutoffs.rare.freqMax) ||
                    (backgroundRatio >= frequencyCutoffs.common.fracBoth && twoSidedPassCount >= frequencyCutoffs.common.freqMin);
        }
        if (samplesOverOneCount > 0) {
            final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
            onesidePass = (genotypeRatio >= frequencyCutoffs.rare.fracSingle && nonRefCount <= frequencyCutoffs.rare.freqMax) ||
                    (genotypeRatio >= frequencyCutoffs.common.fracSingle && nonRefCount >= frequencyCutoffs.common.freqMin);
        }
        final boolean backgroundFail = !(onesidePass || bothsidePass) && hasSplitReadEvidence && nonRefCount > 0;
        final double normalization = maxQuality / Z_DISTRIBUTION.cumulativeProbability(0);
        for (int i = 0; i < samples.size(); i++) {
            final double medianHet;
            final double sdHet;
            if (twoSided[i]) {
                medianHet = bothsideMedianHet;
                sdHet = bothsideSdHet;
            } else {
                medianHet = onesideMedianHet;
                sdHet = onesideSdHet;
            }
            if (genotypes[i] == 0) {
                if (countSum[i] == 0) {
                    genotypeQuals[i] = maxQuality;
                } else {
                    final PoissonDistribution dist = new PoissonDistribution(countSum[i]);
                    genotypeQuals[i] = Math.min((int) Math.round((normalization * (1. - dist.cumulativeProbability(0)))), maxQuality);

                }
            } else if (genotypes[i] == 1) {
                genotypeQuals[i] = Math.min((int) Math.round(normalization * Z_DISTRIBUTION.cumulativeProbability((countSum[i] - medianHet) / sdHet)), maxQuality);
            } else {
                genotypeQuals[i] = Math.min((int) Math.round(normalization * Z_DISTRIBUTION.cumulativeProbability((countSum[i] - 2. * medianHet) / sdHet)), maxQuality);
            }
        }
        final int onesideVariantQuality = computeVariantQuality(onesideNonRefCounts, onesideMedianHet);
        final int bothsideVariantQuality = computeVariantQuality(bothsideNonRefCounts, bothsideMedianHet);
        int variantQuality = 0;
        if (bothsideNonZeroCount > 0) {
            final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
            variantQuality = (int) Math.round(backgroundRatio * bothsideVariantQuality + (1. - backgroundRatio) * onesideVariantQuality);
        } else {
            variantQuality = onesideVariantQuality;
        }
        return new SplitReadGenotypeResult(genotypes, genotypeQuals, variantQuality, bothsidePass, backgroundFail);
    }

    private int computeVariantQuality(final List<Double> nonRefCounts, final double medianHet) {
        final PoissonDistribution poissonDistribution = new PoissonDistribution(medianHet);
        final double normalizationVariant = maxQuality / (-10. * Math.log10(poissonDistribution.cumulativeProbability(0)));
        if (!nonRefCounts.isEmpty()) {
            final double median = MEDIAN.evaluate(nonRefCounts.stream().mapToDouble(Double::doubleValue).toArray());
            return (int) Math.round(Math.min(-10. * normalizationVariant * Math.log10(new PoissonDistribution(median).cumulativeProbability(0)), maxQuality));
        } else {
            return 0;
        }
    }

    /**
     * Computes SR genotypes for a variant and accumulates recovery histogram data.
     * This is the full computation used when writing the training VCF.
     *
     * <p>Note: The {@code pass} flag for histogram binning simplifies in practice to
     * {@code isCNV && nonRefCount > 0 && largeEnough && hasSplitReadEvidence} because
     * PE GQ is clamped to min 1 (making nonRefDiscordantPair always true) and nonRefDepth
     * is subsumed by the depthGenotype null check. The full expression is preserved here
     * for clarity and future-proofing. See {@link #accumulateHistogramOnly} for a
     * streamlined version used when only histogram accumulation is needed.</p>
     */
    public SplitReadGenotypeResult genotypeTraining(final SVCallRecord record,
                                                    final List<SplitReadEvidence> startEvidence,
                                                    final List<SplitReadEvidence> endEvidence,
                                                    final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype,
                                                    final DiscordantPairEvidenceGenotyper.DiscordantPairGenotypeResult discordantPairGenotype,
                                                    final SplitReadGenotypeParameters parameters,
                                                    final List<String> samples) {
        // TODO parameters.minCount and this.countCutoff are redundant ?
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final int[] genotypes = new int[samples.size()];
        final double[] countSum = new double[samples.size()];
        int nonRefCount = 0;
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double startCount = startCounts.getOrDefault(sample, 0.);
            final double endCount = endCounts.getOrDefault(sample, 0.);
            countSum[i] = startCount + endCount;
            if (countSum[i] < parameters.minCount) {
                genotypes[i] = 0;
            } else if (countSum[i] <= parameters.medianHom - parameters.sdHet) {
                genotypes[i] = 1;
            } else {
                genotypes[i] = (int) ((countSum[i] / (parameters.medianHom * 0.5)) + 0.5);
            }
            if (genotypes[i] != 0) {
                ++nonRefCount;
            }
        }
        final boolean largeEnough = record.getLength() != null && record.getLength() >= minSize;
        final boolean hasSplitReadEvidence = record.getEvidence().contains(GATKSVVCFConstants.EvidenceTypes.SR);
        final int minCount = Math.max(trainingCountCutoff / 2, 1);
        final int twoSidedPassCount = countBothSideSupport(startCounts, endCounts, minCount);
        final int bothsideNonZeroCount = countBothSideSupport(startCounts, endCounts, 0);
        final int samplesOverOneCount = countNormalizedSummedSupport(startCounts, endCounts, 1.0);
        // Diagnostics only. In this path the pass flag is per-sample, so record the
        // variant-level approximation that the Phase 2b path uses (depthGenotype != null
        // stands in for isCNV); see accumulateHistogramOnly for why they coincide.
        recordVariantDiagnostics(depthGenotype != null, largeEnough, hasSplitReadEvidence,
                depthGenotype != null && nonRefCount > 0 && largeEnough && hasSplitReadEvidence,
                nonRefCount, samplesOverOneCount, twoSidedPassCount, bothsideNonZeroCount);
        for (int i = 0; i < samples.size(); i++) {
            // Note: nonRefDiscordantPair checks GQ > 0, which is always true since GQ is clamped to min 1.
            // This matches v1.1 behavior (SR_genotype.opt_part1.sh checks $NF>0 on pe.geno.withquality.txt.gz,
            // where PE GQ is also clamped to min 1), making pass/fail effectively per-variant not per-sample.
            final boolean nonRefDiscordantPair = discordantPairGenotype != null && discordantPairGenotype.genotypeQuals()[i] > 0;
            final boolean nonRefDepth = depthGenotype != null && depthGenotype.copyStates()[i] != 2;
            final boolean pass = depthGenotype != null && nonRefCount > 0 && largeEnough && hasSplitReadEvidence && (nonRefDiscordantPair || nonRefDepth);
            if (genotypes[i] > 0) {
                // ALL non-ref samples go into single-sided histogram (matching v1.1's recover.txt
                // which is built from sr.geno.final.oneside.txt.gz covering all variants)
                if (samplesOverOneCount > 0 && nonRefCount > 0) {
                    final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
                    addToSingleHistogram(nonRefCount, genotypeRatio, pass);
                }
                // Variants with bothside support ALSO go into both-sided histogram
                // (matching v1.1's recover.bothsides.txt) and the overlap histogram
                if (bothsideNonZeroCount > 0 && twoSidedPassCount > 0) {
                    final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
                    addToBothHistogram(twoSidedPassCount, backgroundRatio, pass);
                    if (samplesOverOneCount > 0 && nonRefCount > 0) {
                        final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
                        addToOverlapHistogram(nonRefCount, genotypeRatio, twoSidedPassCount, backgroundRatio, pass);
                    }
                }
            }
        }
        return new SplitReadGenotypeResult(genotypes, null, maxQuality, null, null);
    }

    /**
     * Accumulates SR recovery histogram data without computing full genotypes or querying
     * depth/PE evidence. Used by the fast (non-VCF) Phase 2 path in TrainSVGenotyping.
     *
     * <p>The {@code pass} flag is derived from {@code isCNV} (equivalent to
     * {@code depthGenotype != null} since the depth genotyper returns null for non-CNV types)
     * combined with SR-derived metrics. This is equivalent to the full {@link #genotypeTraining}
     * pass flag because PE GQ is clamped to min 1 (making nonRefDiscordantPair always true)
     * and nonRefDepth is subsumed by the depthGenotype null check.</p>
     */
    public void accumulateHistogramOnly(final SVCallRecord record,
                                        final List<SplitReadEvidence> startEvidence,
                                        final List<SplitReadEvidence> endEvidence,
                                        final boolean isCNV,
                                        final SplitReadGenotypeParameters parameters,
                                        final List<String> samples) {
        final Map<String, Double> startCounts = normalizeCounts(startEvidence);
        final Map<String, Double> endCounts = normalizeCounts(endEvidence);
        final int[] genotypes = new int[samples.size()];
        int nonRefCount = 0;
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final double startCount = startCounts.getOrDefault(sample, 0.);
            final double endCount = endCounts.getOrDefault(sample, 0.);
            final double countSum = startCount + endCount;
            if (countSum < parameters.minCount) {
                genotypes[i] = 0;
            } else if (countSum <= parameters.medianHom - parameters.sdHet) {
                genotypes[i] = 1;
            } else {
                genotypes[i] = (int) ((countSum / (parameters.medianHom * 0.5)) + 0.5);
            }
            if (genotypes[i] != 0) {
                ++nonRefCount;
            }
        }
        final boolean largeEnough = record.getLength() != null && record.getLength() >= minSize;
        final boolean hasSplitReadEvidence = record.getEvidence().contains(GATKSVVCFConstants.EvidenceTypes.SR);
        final int minCount = Math.max(trainingCountCutoff / 2, 1);
        final int twoSidedPassCount = countBothSideSupport(startCounts, endCounts, minCount);
        final int bothsideNonZeroCount = countBothSideSupport(startCounts, endCounts, 0);
        final int samplesOverOneCount = countNormalizedSummedSupport(startCounts, endCounts, 1.0);
        final boolean pass = isCNV && nonRefCount > 0 && largeEnough && hasSplitReadEvidence;
        recordVariantDiagnostics(isCNV, largeEnough, hasSplitReadEvidence, pass, nonRefCount,
                samplesOverOneCount, twoSidedPassCount, bothsideNonZeroCount);
        for (int i = 0; i < samples.size(); i++) {
            if (genotypes[i] > 0) {
                // ALL non-ref samples go into single-sided histogram
                if (samplesOverOneCount > 0 && nonRefCount > 0) {
                    final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
                    addToSingleHistogram(nonRefCount, genotypeRatio, pass);
                }
                // Variants with bothside support ALSO go into both-sided + overlap histograms
                if (bothsideNonZeroCount > 0 && twoSidedPassCount > 0) {
                    final double backgroundRatio = twoSidedPassCount / (double) bothsideNonZeroCount;
                    addToBothHistogram(twoSidedPassCount, backgroundRatio, pass);
                    if (samplesOverOneCount > 0 && nonRefCount > 0) {
                        final double genotypeRatio = nonRefCount / (double) samplesOverOneCount;
                        addToOverlapHistogram(nonRefCount, genotypeRatio, twoSidedPassCount, backgroundRatio, pass);
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------------
    // Diagnostics recording and reporting
    // ------------------------------------------------------------------------

    private static int countBin(final int count) {
        if (count <= 5) {
            return Math.max(count, 0);
        } else if (count <= 7) {
            return 6;
        } else if (count <= 15) {
            return 7;
        } else if (count <= 31) {
            return 8;
        } else if (count <= 63) {
            return 9;
        } else if (count <= 127) {
            return 10;
        } else if (count <= 255) {
            return 11;
        } else if (count <= 511) {
            return 12;
        }
        return 13;
    }

    /**
     * Tally one variant's contribution to the cutoff histograms. Observational only.
     *
     * <p>Counts are variant-level, whereas the frac histograms are incremented once per
     * non-ref sample, so histogram totals will exceed these tallies. The purpose here is
     * to attribute where variants are lost and to expose whether the frac values are
     * saturating at or above 1.0, which the 11-bin histograms alone cannot show.</p>
     */
    private void recordVariantDiagnostics(final boolean isCNV, final boolean largeEnough,
                                          final boolean hasSplitReadEvidence, final boolean pass,
                                          final int nonRefCount, final int samplesOverOneCount,
                                          final int twoSidedPassCount, final int bothsideNonZeroCount) {
        diagVariantsSeen++;
        if (!isCNV) {
            diagVariantsNotCnv++;
        }
        if (!largeEnough) {
            diagVariantsTooSmall++;
        }
        if (!hasSplitReadEvidence) {
            diagVariantsNoSrEvidence++;
        }
        if (pass) {
            diagVariantsPassFlag++;
        }
        diagNonRefCountHist[countBin(nonRefCount)]++;
        diagSamplesOverOneHist[countBin(samplesOverOneCount)]++;
        diagTwoSidedPassHist[countBin(twoSidedPassCount)]++;
        diagBothsideNonZeroHist[countBin(bothsideNonZeroCount)]++;
        if (nonRefCount == 0) {
            diagVariantsNoNonRef++;
            return;
        }
        if (samplesOverOneCount == 0) {
            diagVariantsNoSamplesOverOne++;
        } else {
            diagVariantsContributedSingle++;
            final double frac = nonRefCount / (double) samplesOverOneCount;
            diagSingleFracSum += frac;
            diagSingleFracN++;
            if (frac > 1.0) {
                diagSingleFracGtOne++;
            } else if (frac == 1.0) {
                diagSingleFracEqOne++;
            } else if (frac >= 0.95) {
                diagSingleFracNearOne++;
            }
        }
        if (bothsideNonZeroCount == 0) {
            diagVariantsNoBothsideNonZero++;
        } else if (twoSidedPassCount == 0) {
            diagVariantsNoTwoSidedPass++;
        } else {
            diagVariantsContributedBoth++;
            final double frac = twoSidedPassCount / (double) bothsideNonZeroCount;
            diagBothFracSum += frac;
            diagBothFracN++;
            if (frac > 1.0) {
                diagBothFracGtOne++;
            } else if (frac == 1.0) {
                diagBothFracEqOne++;
            } else if (frac >= 0.95) {
                diagBothFracNearOne++;
            }
        }
        // A non-ref variant whose count falls in neither the rare nor the common frequency
        // range contributes nothing to either histogram, regardless of its frac.
        final boolean inRare = nonRefCount > rareMin && nonRefCount <= rareMax;
        final boolean inCommon = nonRefCount > commonMin && nonRefCount <= commonMax;
        if (!inRare && !inCommon) {
            diagVariantsOutsideBothFreqBins++;
        }
    }

    private static void appendCountHistogram(final StringBuilder sb, final String name, final long[] histogram) {
        for (int i = 0; i < histogram.length; i++) {
            sb.append(name).append('\t').append(COUNT_BIN_LABELS[i]).append('\t').append(histogram[i]).append('\n');
        }
    }

    private static void appendFracHistogram(final StringBuilder sb, final String name, final long[][][] histogram) {
        for (int f = 0; f < 2; f++) {
            for (int p = 0; p < 2; p++) {
                sb.append(name).append('\t').append(f == 0 ? "rare" : "common").append('\t')
                        .append(p == 0 ? "pass" : "fail");
                for (int b = 0; b < NUM_FRAC_BINS; b++) {
                    sb.append('\t').append(histogram[f][p][b]);
                }
                sb.append('\n');
            }
        }
    }

    private static void appendGrid(final StringBuilder sb, final String name, final CutoffGrid grid) {
        if (grid == null) {
            sb.append(name).append("\tNOT_COMPUTED\n");
            return;
        }
        sb.append("# ").append(name).append(": baseline=").append(grid.baseline())
                .append(" freq_min=").append(grid.freqMin()).append(" freq_max=").append(grid.freqMax())
                .append(" selected_index=").append(grid.selectedIndex()).append('\n');
        sb.append(name).append("\tindex\tfrac_single\tfrac_both\tpass\tfail\tfail_frac_a\tpass_ratio_b\tscore\tselected\n");
        for (int i = 0; i < grid.cells().size(); i++) {
            final CutoffResult c = grid.cells().get(i);
            final double a = c.countFail() / (double) (c.countFail() + c.countPass());
            final double b = (c.countPass() / grid.baseline()) - 1;
            sb.append(name).append('\t').append(i)
                    .append('\t').append(String.format("%.1f", c.fracSingle()))
                    .append('\t').append(String.format("%.1f", c.fracBoth()))
                    .append('\t').append(c.countPass())
                    .append('\t').append(c.countFail())
                    .append('\t').append(a)
                    .append('\t').append(b)
                    .append('\t').append(grid.scores()[i])
                    .append('\t').append(i == grid.selectedIndex() ? "*" : "")
                    .append('\n');
        }
    }

    /**
     * Summarize how the argmax arrived at its answer: how many cells tied at the maximum,
     * how many scored NaN, and how the selected cell compares to the (0.0, 0.0) cell.
     *
     * <p>These are the quantities that distinguish a real optimum from a degenerate grid.
     * An all-NaN or all-tied score array means the objective carried no signal;
     * {@link #validateCutoffGrid} rejects those rather than resolving them to zero cutoffs.</p>
     */
    private static void appendSelectionSummary(final StringBuilder sb, final String name, final CutoffGrid grid) {
        if (grid == null) {
            sb.append(name).append("\tNOT_COMPUTED\n");
            return;
        }
        final double[] scores = grid.scores();
        int nanCount = 0;
        double max = Double.NEGATIVE_INFINITY;
        for (final double s : scores) {
            if (Double.isNaN(s)) {
                nanCount++;
            } else if (s > max) {
                max = s;
            }
        }
        int tiedAtMax = 0;
        for (final double s : scores) {
            if (s == max) {
                tiedAtMax++;
            }
        }
        // Best cell that is not the all-zero cell, for comparison against what was chosen
        int bestNonZeroIndex = -1;
        for (int i = 0; i < scores.length; i++) {
            final CutoffResult c = grid.cells().get(i);
            if (c.fracSingle() == 0 && c.fracBoth() == 0) {
                continue;
            }
            if (!Double.isNaN(scores[i]) && (bestNonZeroIndex < 0 || scores[i] > scores[bestNonZeroIndex])) {
                bestNonZeroIndex = i;
            }
        }
        final boolean selectionMade = grid.selectedIndex() != NO_CUTOFF_SELECTED;
        sb.append(name).append("\tselection_made\t").append(selectionMade).append('\n');
        sb.append(name).append("\tnum_cells\t").append(scores.length).append('\n');
        sb.append(name).append("\tnum_nan_scores\t").append(nanCount).append('\n');
        // Either of the next two being true means the argmax had nothing to discriminate on
        // and fell through to index 0, which is the (0.0, 0.0) cell.
        sb.append(name).append("\tall_scores_nan\t").append(nanCount == scores.length).append('\n');
        sb.append(name).append("\tall_scores_tied\t").append(tiedAtMax == scores.length).append('\n');
        sb.append(name).append("\tmax_score\t").append(max).append('\n');
        sb.append(name).append("\tnum_cells_tied_at_max\t").append(tiedAtMax).append('\n');
        sb.append(name).append("\tselected_index\t").append(grid.selectedIndex()).append('\n');
        if (selectionMade) {
            final CutoffResult selected = grid.cells().get(grid.selectedIndex());
            sb.append(name).append("\tselected_frac_single\t").append(selected.fracSingle()).append('\n');
            sb.append(name).append("\tselected_frac_both\t").append(selected.fracBoth()).append('\n');
            sb.append(name).append("\tselected_is_all_zero_cell\t")
                    .append(selected.fracSingle() == 0 && selected.fracBoth() == 0).append('\n');
        }
        sb.append(name).append("\tscore_at_zero_cell\t").append(scores[0]).append('\n');
        if (bestNonZeroIndex >= 0) {
            final CutoffResult alt = grid.cells().get(bestNonZeroIndex);
            sb.append(name).append("\tbest_nonzero_index\t").append(bestNonZeroIndex).append('\n');
            sb.append(name).append("\tbest_nonzero_frac_single\t").append(alt.fracSingle()).append('\n');
            sb.append(name).append("\tbest_nonzero_frac_both\t").append(alt.fracBoth()).append('\n');
            sb.append(name).append("\tbest_nonzero_score\t").append(scores[bestNonZeroIndex]).append('\n');
        }
    }

    /**
     * Render the full SR cutoff diagnostic report. Intended to be written to a file; the
     * report is section-delimited with {@code ##} headers and tab-separated bodies.
     */
    public String cutoffDiagnosticsReport() {
        final StringBuilder sb = new StringBuilder(1 << 16);

        sb.append("## SR_MODEL_PARAMETERS\n");
        sb.append("training_count_cutoff_sr_count\t").append(trainingCountCutoff).append('\n');
        sb.append("two_sided_training_threshold\t").append(Math.max(trainingCountCutoff / 2, 1)).append('\n');
        sb.append("min_size\t").append(minSize).append('\n');
        sb.append("target_coverage\t").append(targetCoverage).append('\n');
        sb.append("max_quality\t").append(maxQuality).append('\n');
        sb.append("num_samples_in_coverage_map\t").append(sampleCoverageMap.size()).append('\n');
        sb.append("rare_min\t").append(rareMin).append('\n');
        sb.append("rare_max\t").append(rareMax).append('\n');
        sb.append("common_min\t").append(commonMin).append('\n');
        sb.append("common_max\t").append(commonMax).append('\n');

        sb.append("## SR_TRAINING_PASSES\n");
        sb.append("first_pass_variants\t").append(diagFirstPassVariants).append('\n');
        sb.append("first_pass_het_observations\t").append(diagFirstPassHetN).append('\n');
        sb.append("first_pass_het_median\t").append(hetMedian).append('\n');
        sb.append("first_pass_het_mad_raw\t").append(hetMad).append('\n');
        sb.append("first_pass_het_cutoff\t").append(hetCutoff).append('\n');
        sb.append("second_pass_het_observations\t").append(diagSecondPassHetN).append('\n');
        sb.append("second_pass_hom_observations\t").append(diagSecondPassHomN).append('\n');
        sb.append("second_pass_het_median\t").append(diagSecondPassHetMedian).append('\n');
        sb.append("second_pass_hom_median\t").append(diagSecondPassHomMedian).append('\n');
        sb.append("second_pass_het_mad_raw\t").append(diagSecondPassHetMad).append('\n');

        sb.append("## SR_VARIANT_TALLIES\n");
        sb.append("variants_seen\t").append(diagVariantsSeen).append('\n');
        sb.append("variants_not_cnv\t").append(diagVariantsNotCnv).append('\n');
        sb.append("variants_below_min_size\t").append(diagVariantsTooSmall).append('\n');
        sb.append("variants_without_sr_evidence\t").append(diagVariantsNoSrEvidence).append('\n');
        sb.append("variants_with_pass_flag\t").append(diagVariantsPassFlag).append('\n');
        sb.append("variants_no_nonref_samples\t").append(diagVariantsNoNonRef).append('\n');
        sb.append("variants_no_samples_over_one\t").append(diagVariantsNoSamplesOverOne).append('\n');
        sb.append("variants_no_bothside_nonzero\t").append(diagVariantsNoBothsideNonZero).append('\n');
        sb.append("variants_no_two_sided_pass\t").append(diagVariantsNoTwoSidedPass).append('\n');
        sb.append("variants_contributed_single\t").append(diagVariantsContributedSingle).append('\n');
        sb.append("variants_contributed_both\t").append(diagVariantsContributedBoth).append('\n');
        sb.append("variants_outside_both_freq_bins\t").append(diagVariantsOutsideBothFreqBins).append('\n');

        sb.append("## SR_FRAC_SATURATION\n");
        sb.append("# frac > 1 is possible when the numerator threshold is at or below the\n");
        sb.append("# denominator threshold; the 11-bin histograms cannot distinguish it from frac == 1.\n");
        sb.append("single_frac_variants\t").append(diagSingleFracN).append('\n');
        sb.append("single_frac_mean\t").append(diagSingleFracN == 0 ? "NA" : diagSingleFracSum / diagSingleFracN).append('\n');
        sb.append("single_frac_gt_one\t").append(diagSingleFracGtOne).append('\n');
        sb.append("single_frac_eq_one\t").append(diagSingleFracEqOne).append('\n');
        sb.append("single_frac_in_0.95_to_1\t").append(diagSingleFracNearOne).append('\n');
        sb.append("both_frac_variants\t").append(diagBothFracN).append('\n');
        sb.append("both_frac_mean\t").append(diagBothFracN == 0 ? "NA" : diagBothFracSum / diagBothFracN).append('\n');
        sb.append("both_frac_gt_one\t").append(diagBothFracGtOne).append('\n');
        sb.append("both_frac_eq_one\t").append(diagBothFracEqOne).append('\n');
        sb.append("both_frac_in_0.95_to_1\t").append(diagBothFracNearOne).append('\n');

        sb.append("## SR_COUNT_DISTRIBUTIONS\n");
        sb.append("# metric\tcount_bin\tnum_variants\n");
        appendCountHistogram(sb, "non_ref_count", diagNonRefCountHist);
        appendCountHistogram(sb, "samples_over_one_count", diagSamplesOverOneHist);
        appendCountHistogram(sb, "two_sided_pass_count", diagTwoSidedPassHist);
        appendCountHistogram(sb, "bothside_non_zero_count", diagBothsideNonZeroHist);

        sb.append("## SR_FRAC_HISTOGRAMS\n");
        sb.append("# Columns are frac bins: [0,0.1) [0.1,0.2) ... [0.9,1.0) and >=1.0 (bin 10).\n");
        sb.append("# Entries are sample-level, so totals exceed the variant tallies above.\n");
        sb.append("# histogram\tfreq_bin\tclass\tbin0..bin10\n");
        appendFracHistogram(sb, "single", singleHistogram);
        appendFracHistogram(sb, "both", bothHistogram);

        sb.append("## SR_OVERLAP_HISTOGRAM\n");
        sb.append("# freq_bin\tclass\tsingle_bin\tbin0..bin10 over both_bin\n");
        for (int f = 0; f < 2; f++) {
            for (int p = 0; p < 2; p++) {
                for (int s = 0; s < NUM_FRAC_BINS; s++) {
                    sb.append("overlap\t").append(f == 0 ? "rare" : "common").append('\t')
                            .append(p == 0 ? "pass" : "fail").append('\t').append(s);
                    for (int b = 0; b < NUM_FRAC_BINS; b++) {
                        sb.append('\t').append(overlapHistogram[f][p][s][b]);
                    }
                    sb.append('\n');
                }
            }
        }

        sb.append("## SR_CUTOFF_GRID\n");
        appendGrid(sb, "rare", diagRareGrid);
        appendGrid(sb, "common", diagCommonGrid);

        sb.append("## SR_SELECTION_SUMMARY\n");
        appendSelectionSummary(sb, "rare", diagRareGrid);
        appendSelectionSummary(sb, "common", diagCommonGrid);

        // Machine-readable status, parsed by the ValidateSRCutoffs task. Keep the key names and
        // the SelectionStatus spellings stable; that task greps for them.
        sb.append("## SR_SELECTION_STATUS\n");
        appendSelectionStatus(sb, "rare", diagRareOutcome);
        appendSelectionStatus(sb, "common", diagCommonOutcome);

        return sb.toString();
    }

    private static void appendSelectionStatus(final StringBuilder sb, final String name,
                                              final SelectionOutcome outcome) {
        final SelectionStatus status = outcome == null ? SelectionStatus.OK : outcome.status();
        sb.append(name).append("_selection_status\t").append(status).append('\n');
        if (outcome != null && outcome.rejected()) {
            sb.append(name).append("_selection_rejection_reason\t").append(outcome.detail()).append('\n');
        }
    }

    /**
     * Selection outcomes for the rare and common frequency bins, in that order. Empty until
     * {@link #finalizeThirdPass()} has run. Callers should log any rejection; enforcement
     * happens downstream, once the diagnostics report has been delocalized.
     */
    public List<SelectionOutcome> cutoffSelectionOutcomes() {
        if (diagRareOutcome == null || diagCommonOutcome == null) {
            return Collections.emptyList();
        }
        return List.of(diagRareOutcome, diagCommonOutcome);
    }

    public boolean trainableRecord(final SVCallRecord record,
                                   final boolean discordantPairEligible,
                                   final SVStratificationEngine exclusionEngine) {
        if (!discordantPairEligible) {
            return false;
        }
        if (minSize != null && (record.getLength() == null || record.getLength() < minSize)) {
            return false;
        }
        if (exclusionEngine != null && !exclusionEngine.getMatches(record, 0, 1, 1).isEmpty()) {
            return false;
        }
        if (!record.getEvidence().contains(GATKSVVCFConstants.EvidenceTypes.SR)) {
            return false;
        }
        return true;
    }

    public boolean trainableRecord(final SVCallRecord record,
                                   final DiscordantPairEvidenceGenotyper discordantPairGenotyper,
                                   final SVStratificationEngine exclusionEngine) {
        return trainableRecord(record, discordantPairGenotyper.isTrainingRecord(record), exclusionEngine);
    }

    private static final class FirstPassResult {
        private final Map<String, Double> counts;
        private final Map<String, Integer> depthGenotypes;
        public FirstPassResult(final Set<String> passingSamples, final Map<String, Double> startCounts, final Map<String, Double> endCounts, final DepthEvidenceGenotyper.DepthGenotypeResult depthGenotype, final List<String> depthGenotypeSamples) {
            this.counts = new HashMap<>();
            final Set<String> keys = new HashSet<>(startCounts.keySet());
            keys.addAll(endCounts.keySet());
            for (final String key : keys) {
                if (passingSamples.contains(key)) {
                    final double startCount = startCounts.getOrDefault(key, 0.);
                    final double endCount = endCounts.getOrDefault(key, 0.);
                    counts.put(key, startCount + endCount);
                }
            }
            this.depthGenotypes = new HashMap<>();
            for (int i = 0; i < depthGenotypeSamples.size(); i++) {
                final String sample = depthGenotypeSamples.get(i);
                if (counts.containsKey(sample)) {
                    depthGenotypes.put(sample, depthGenotype.copyStates()[i]);
                }
            }
        }

        public Double getCount(final String sample) {
            return counts.get(sample);
        }

        public double[] getCounts(final Set<Integer> validDepthStates) {
            final List<Double> countsList = new ArrayList<>(counts.size());
            for (final Map.Entry<String, Double> entry : counts.entrySet()) {
                if (depthGenotypes.containsKey(entry.getKey())) {
                    final int copyState = depthGenotypes.get(entry.getKey());
                    if (validDepthStates.contains(copyState)) {
                        countsList.add(entry.getValue());
                    }
                } else {
                    throw new IllegalArgumentException("Sample '" + entry.getKey() + "' does not exist in depth genotype");
                }
            }
            return countsList.stream().mapToDouble(Double::doubleValue).toArray();
        }
    }

    public record SplitReadGenotypeResult(int[] genotypes, int[] genotypeQuals, double variantQual, Boolean bothsidePass, Boolean backgroundFail) {}
    public record SplitReadGenotypeParameters(double minCount, double medianHom, double sdHet) {}
    public record SplitReadGenotypeFrequencyCutoffs(CutoffResult rare, CutoffResult common) {}

    public static final class SplitReadTableParser {

        private static final String MIN_COUNT_COLUMN = "sr_count";
        private static final String MEDIAN_HOM_COLUMN = "median_hom";
        private static final String SD_HET_COLUMN = "sd_het";
        private static final String RARE_MIN_COLUMN = "rare_min";
        private static final String RARE_MAX_COLUMN = "rare_max";
        private static final String COMMON_MIN_COLUMN = "common_min";
        private static final String COMMON_MAX_COLUMN = "common_max";
        private static final String RARE_PASS_COLUMN = "rare_pass";
        private static final String RARE_FAIL_COLUMN = "rare_fail";
        private static final String COMMON_PASS_COLUMN = "common_pass";
        private static final String COMMON_FAIL_COLUMN = "common_fail";
        private static final String RARE_SINGLE_COLUMN = "rare_single";
        private static final String RARE_BOTH_COLUMN = "rare_both";
        private static final String COMMON_SINGLE_COLUMN = "common_single";
        private static final String COMMON_BOTH_COLUMN = "common_both";
        public static final TableColumnCollection CUTOFFS_COLUMNS = new TableColumnCollection(Arrays.asList(
                MIN_COUNT_COLUMN, MEDIAN_HOM_COLUMN, SD_HET_COLUMN, RARE_MIN_COLUMN, RARE_MAX_COLUMN, COMMON_MIN_COLUMN,
                COMMON_MAX_COLUMN, RARE_PASS_COLUMN, RARE_FAIL_COLUMN, COMMON_PASS_COLUMN, COMMON_FAIL_COLUMN,
                RARE_SINGLE_COLUMN, RARE_BOTH_COLUMN, COMMON_SINGLE_COLUMN, COMMON_BOTH_COLUMN));

        public void composeCutoffsLine(final SplitReadGenotypeMetrics stats,
                                        final DataLine dataLine) {
            dataLine.append(stats.parameters().minCount());
            dataLine.append(stats.parameters().medianHom());
            dataLine.append(stats.parameters().sdHet());
            dataLine.append(stats.cutoffs().rare().freqMin());
            dataLine.append(stats.cutoffs().rare().freqMax());
            dataLine.append(stats.cutoffs().common().freqMin());
            dataLine.append(stats.cutoffs().common().freqMax());
            dataLine.append(stats.cutoffs().rare().countPass());
            dataLine.append(stats.cutoffs().rare().countFail());
            dataLine.append(stats.cutoffs().common().countPass());
            dataLine.append(stats.cutoffs().common().countFail());
            dataLine.append(stats.cutoffs().rare().fracSingle());
            dataLine.append(stats.cutoffs().rare().fracBoth());
            dataLine.append(stats.cutoffs().common().fracSingle());
            dataLine.append(stats.cutoffs().common().fracBoth());
        }

        public Function<DataLine, SplitReadGenotypeMetrics> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
            // Check for expected columns
            for (final String column : CUTOFFS_COLUMNS.names()) {
                if (!columns.contains(column)) {
                    throw exceptionFactory.apply("Missing column " + column);
                }
            }
            // Check there are no extra columns
            if (columns.columnCount() != CUTOFFS_COLUMNS.columnCount()) {
                throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
            }
            return this::parseTableLine;
        }

        public SplitReadGenotypeMetrics parseTableLine(final DataLine dataLine) {
            final double minCount = Double.parseDouble(dataLine.get(0));
            final double medianHom = Double.parseDouble(dataLine.get(1));
            final double sdHet = Double.parseDouble(dataLine.get(2));
            final int rareMin = Integer.parseInt(dataLine.get(3));
            final int rareMax = Integer.parseInt(dataLine.get(4));
            final int commonMin = Integer.parseInt(dataLine.get(5));
            final int commonMax = Integer.parseInt(dataLine.get(6));
            final long rarePass = Long.parseLong(dataLine.get(7));
            final long rareFail = Long.parseLong(dataLine.get(8));
            final long commonPass = Long.parseLong(dataLine.get(9));
            final long commonFail = Long.parseLong(dataLine.get(10));
            final double rareSingle = Double.parseDouble(dataLine.get(11));
            final double rareBoth = Double.parseDouble(dataLine.get(12));
            final double commonSingle = Double.parseDouble(dataLine.get(13));
            final double commonBoth = Double.parseDouble(dataLine.get(14));
            return new SplitReadGenotypeMetrics(
                    new SplitReadEvidenceGenotyper.SplitReadGenotypeParameters(minCount, medianHom, sdHet),
                    new SplitReadEvidenceGenotyper.SplitReadGenotypeFrequencyCutoffs(
                            new SplitReadEvidenceGenotyper.CutoffResult(rareSingle, rareBoth, rarePass, rareFail, rareMin, rareMax),
                            new SplitReadEvidenceGenotyper.CutoffResult(commonSingle, commonBoth, commonPass, commonFail, commonMin, commonMax)
                    )
            );
        }
    }

    public record SplitReadGenotypeMetrics(SplitReadGenotypeParameters parameters,
                                           SplitReadGenotypeFrequencyCutoffs cutoffs) {}
}
