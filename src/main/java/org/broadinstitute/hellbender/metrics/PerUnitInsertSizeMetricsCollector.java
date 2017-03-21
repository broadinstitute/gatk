package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.*;

/**
 * A Collector for individual InsertSizeMetrics for a given SAMPLE or SAMPLE/LIBRARY or
 * SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels)
 */
public final class PerUnitInsertSizeMetricsCollector
        implements PerUnitMetricCollector<InsertSizeMetrics, Integer, InsertSizeMetricsCollectorArgs>,
                    Serializable
{
    private static final long serialVersionUID = 1L;

    // kryo can't serialize an EnumMap so use a LinkedHashMap (also so we can maintain order
    // the insertion order of entries since the combining code below counts on it
    private final Map<SamPairUtil.PairOrientation, Histogram<Integer>> histograms = new LinkedHashMap<>();

    private final String sample;
    private final String library;
    private final String readGroup;

    // When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this
    // percentage of overall reads. (Range: 0 to 1)
    private final double minimumPct;

    // Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION.
    // This is done because insert size data typically includes enough anomolous values from chimeras and other
    // artifacts to make the mean and sd grossly misleading regarding the real distribution.
    private final double deviations;

    //Explicitly sets the Histogram width if set, overriding automatic truncation of Histogram tail.
    //Also, when calculating mean and stdev, only bins <= HISTOGRAM_WIDTH will be included.
    private Integer histogramWidth;

    private double totalInserts = 0;

    public PerUnitInsertSizeMetricsCollector(
            final String sample,
            final String library,
            final String readGroup,
            final double minimumPct,
            final double deviations,
            final Integer histogramWidth) {

        this.sample = sample;
        this.library = library;
        this.readGroup = readGroup;
        this.minimumPct = minimumPct;
        this.histogramWidth = histogramWidth;
        this.deviations = deviations;

        String prefix = createHistogramValuePrefix();

        histograms.put(SamPairUtil.PairOrientation.FR,     new Histogram<>("insert_size", prefix + "fr_count"));
        histograms.put(SamPairUtil.PairOrientation.TANDEM, new Histogram<>("insert_size", prefix + "tandem_count"));
        histograms.put(SamPairUtil.PairOrientation.RF,     new Histogram<>("insert_size", prefix + "rf_count"));
    }

    @Override
    public void acceptRecord(final InsertSizeMetricsCollectorArgs args) {
        histograms.get(args.getPairOrientation()).increment(args.getInsertSize());
    }

    @Override
    public void finish() { }

    public double getTotalInserts() {
        return totalInserts;
    }

    @Override
    public void addMetricsToFile(final MetricsFile<InsertSizeMetrics, Integer> file) {

        for(final Map.Entry<SamPairUtil.PairOrientation, Histogram<Integer>> entry : histograms.entrySet()) {
            totalInserts += entry.getValue().getCount();
        }

        for(final Map.Entry<SamPairUtil.PairOrientation, Histogram<Integer>> entry : histograms.entrySet()) {

            final SamPairUtil.PairOrientation pairOrientation = entry.getKey();
            final Histogram<Integer> Histogram = entry.getValue();

            final double total = Histogram.getCount();

            // Only include a category if it has a sufficient percentage of the data in it
            if( total > totalInserts * minimumPct ) {
                final InsertSizeMetrics metrics = new InsertSizeMetrics();
                metrics.SAMPLE             = this.sample;
                metrics.LIBRARY            = this.library;
                metrics.READ_GROUP         = this.readGroup;
                metrics.PAIR_ORIENTATION   = pairOrientation;
                metrics.READ_PAIRS         = (long) total;
                metrics.MAX_INSERT_SIZE    = (int) Histogram.getMax();
                metrics.MIN_INSERT_SIZE    = (int) Histogram.getMin();
                metrics.MEDIAN_INSERT_SIZE = Histogram.getMedian();
                metrics.MEDIAN_ABSOLUTE_DEVIATION = Histogram.getMedianAbsoluteDeviation();

                final double median  = Histogram.getMedian();
                double covered = 0;
                double low  = median;
                double high = median;

                while (low >= Histogram.getMin() || high <= Histogram.getMax()) {
                    final Histogram.Bin<Integer> lowBin = Histogram.get((int) low);
                    if (lowBin != null) covered += lowBin.getValue();

                    if (low != high) {
                        final Histogram.Bin<Integer> highBin = Histogram.get((int) high);
                        if (highBin != null) covered += highBin.getValue();
                    }

                    final double percentCovered = covered / total;
                    final int distance = (int) (high - low) + 1;
                    if (percentCovered >= 0.1  && metrics.WIDTH_OF_10_PERCENT == 0) metrics.WIDTH_OF_10_PERCENT = distance;
                    if (percentCovered >= 0.2  && metrics.WIDTH_OF_20_PERCENT == 0) metrics.WIDTH_OF_20_PERCENT = distance;
                    if (percentCovered >= 0.3  && metrics.WIDTH_OF_30_PERCENT == 0) metrics.WIDTH_OF_30_PERCENT = distance;
                    if (percentCovered >= 0.4  && metrics.WIDTH_OF_40_PERCENT == 0) metrics.WIDTH_OF_40_PERCENT = distance;
                    if (percentCovered >= 0.5  && metrics.WIDTH_OF_50_PERCENT == 0) metrics.WIDTH_OF_50_PERCENT = distance;
                    if (percentCovered >= 0.6  && metrics.WIDTH_OF_60_PERCENT == 0) metrics.WIDTH_OF_60_PERCENT = distance;
                    if (percentCovered >= 0.7  && metrics.WIDTH_OF_70_PERCENT == 0) metrics.WIDTH_OF_70_PERCENT = distance;
                    if (percentCovered >= 0.8  && metrics.WIDTH_OF_80_PERCENT == 0) metrics.WIDTH_OF_80_PERCENT = distance;
                    if (percentCovered >= 0.9  && metrics.WIDTH_OF_90_PERCENT == 0) metrics.WIDTH_OF_90_PERCENT = distance;
                    if (percentCovered >= 0.99 && metrics.WIDTH_OF_99_PERCENT == 0) metrics.WIDTH_OF_99_PERCENT = distance;

                    --low;
                    ++high;
                }

                // Trim the Histogram down to get rid of outliers that would make the chart useless.
                final Histogram<Integer> trimmedHisto = Histogram; //alias it
                if (histogramWidth == null) {
                    histogramWidth = (int) (metrics.MEDIAN_INSERT_SIZE + (deviations * metrics.MEDIAN_ABSOLUTE_DEVIATION));
                }

                trimmedHisto.trimByWidth(histogramWidth);

                metrics.MEAN_INSERT_SIZE = trimmedHisto.getMean();
                metrics.STANDARD_DEVIATION = trimmedHisto.getStandardDeviation();

                file.addHistogram(trimmedHisto);
                file.addMetric(metrics);
            }
        }
    }

    /**
     * Combine this PerUnitInsertSizeMetricsCollector with sourceCollector and return a combined
     * PerUnitInsertSizeMetricsCollector.
     * @param sourceCollector PerUnitInsertSizeMetricsCollector to combine in
     * @return PerUnitInsertSizeMetricsCollector representing the combination of the source collector with this collector
     */
    public PerUnitInsertSizeMetricsCollector combine(final PerUnitInsertSizeMetricsCollector sourceCollector) {
        Utils.nonNull(sourceCollector);
        final String validationMessage = "Internal error combining collectors";
        validateEquals(this.sample, sourceCollector.sample, validationMessage);
        validateEquals(this.library, sourceCollector.library, validationMessage);
        validateEquals(this.readGroup, sourceCollector.readGroup, validationMessage);

        PerUnitInsertSizeMetricsCollector combinedCollector = new PerUnitInsertSizeMetricsCollector(
                this.sample,
                this.library,
                this.readGroup,
                this.minimumPct,
                this.deviations,
                this.histogramWidth);
        combinedCollector.totalInserts = this.totalInserts + sourceCollector.totalInserts;

        // combine the histograms from each collector into the new combined collector
        // each collector has an entry for each pair orientation added in the initialization
        // code above; though any given entry may be empty
        this.histograms.forEach(
                (po, targetHist) -> {
                    Histogram<Integer> sourceHist = sourceCollector.histograms.get(po);
                    // as a check, make sure the bin and value labels for the histograms for this PO
                    // are the same in both histograms
                    Utils.validate(targetHist.getBinLabel().equals(sourceHist.getBinLabel()) &&
                            targetHist.getValueLabel().equals(sourceHist.getValueLabel()), "Internal error combining collectors: attempt to combine mismatched histograms");

                    Histogram<Integer> combinedHist = new Histogram<>(targetHist.getBinLabel(), targetHist.getValueLabel());
                    combinedHist.addHistogram(sourceHist);
                    combinedHist.addHistogram(targetHist);

                    combinedCollector.histograms.put(po, combinedHist);
                }
        );
        return combinedCollector;
    }

    // Safely validate that two strings are equal, even if one or both are null.
    private static void validateEquals(final String source, final String target, final String message) {
        Utils.validate( Objects.equals(source, target), () -> String.format("%s (%s : %s)",
                            message,
                            source == null ? "null" : source,
                            target == null ? "null" : target));
    }

    private String createHistogramValuePrefix() {
        String prefix;
        if (this.readGroup != null) {
            prefix = this.readGroup + ".";
        }
        else if (this.library != null) {
            prefix = this.library + ".";
        }
        else if (this.sample != null) {
            prefix = this.sample + ".";
        }
        else {
            prefix = "All_Reads.";
        }
        return prefix;
    }

}
