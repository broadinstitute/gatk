package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.insertsize;

import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.values.KV;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.HistogramDataflow;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.MetricsFileDataflow;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * CombineFN to combine the {@link HistogramDataflow}s of insert sizes into a single {@link MetricsFileDataflow}
 */
final class CombineInsertSizeHistogramsIntoMetricsFile
        extends Combine.CombineFn<KV<InsertSizeAggregationLevel, HistogramDataflow<Integer>>, CombineInsertSizeHistogramsIntoMetricsFile, MetricsFileDataflow<InsertSizeMetrics, Integer>> {
    private final static long serialVersionUID = 1l;
    private static final String BIN_LABEL = "insert_size";

    private final double DEVIATIONS;
    private final Integer HISTOGRAM_WIDTH;
    private final float MINIMUM_PERCENT;

    private final Map<InsertSizeAggregationLevel, HistogramDataflow<Integer>> histograms = new HashMap<>();

    public CombineInsertSizeHistogramsIntoMetricsFile(double deviations, Integer histogramWidth, float minimumPercent) {
        this.DEVIATIONS = deviations;
        this.HISTOGRAM_WIDTH = histogramWidth;
        this.MINIMUM_PERCENT = minimumPercent;
    }

    public void addHistogramToMetricsFile(final InsertSizeAggregationLevel aggregationLevel, final HistogramDataflow<Integer> histogram, final MetricsFileDataflow<InsertSizeMetrics, Integer> metricsFile, final double totalInserts) {
        if (histogramHasEnoughInsertsToInclude(histogram, totalInserts)) {
            final InsertSizeMetrics metrics = initializeInsertSizeMetric(aggregationLevel, histogram);

            // Trim the Histogram down to get rid of outliers that would make the chart useless.
            histogram.trimByWidth(inferHistogramWidth(HISTOGRAM_WIDTH, metrics, DEVIATIONS));

            metrics.MEAN_INSERT_SIZE = histogram.getMean();
            metrics.STANDARD_DEVIATION = histogram.getStandardDeviation();

            histogram.setBinLabel(BIN_LABEL);
            histogram.setValueLabel(aggregationLevel.getValueLabel());
            metricsFile.addHistogram(histogram);
            metricsFile.addMetric(metrics);
        }
    }

    /**
     * Decide if it's worth including this histogram in the output
     * @param totalInserts count of all the inserts from all histograms in this aggregation level
     * @return is the number of inserts in this histogram greater than {@link #MINIMUM_PERCENT} of all the total number of inserts in this aggregation level
     */
    private boolean histogramHasEnoughInsertsToInclude(HistogramDataflow<?> histogram, double totalInserts) {
        return histogram.getCount() > totalInserts * MINIMUM_PERCENT;
    }


    /**
     * Create
     */
    private static InsertSizeMetrics initializeInsertSizeMetric(InsertSizeAggregationLevel aggregationLevel, HistogramDataflow<Integer> histogram) {
        final InsertSizeMetrics metrics = new InsertSizeMetrics();
        metrics.SAMPLE = aggregationLevel.getSample();
        metrics.LIBRARY = aggregationLevel.getLibrary();
        metrics.READ_GROUP = aggregationLevel.getReadGroup();
        metrics.PAIR_ORIENTATION = aggregationLevel.getOrientation();

        metrics.READ_PAIRS = (long) histogram.getCount();
        metrics.MAX_INSERT_SIZE = (int) histogram.getMax();
        metrics.MIN_INSERT_SIZE = (int) histogram.getMin();
        metrics.MEDIAN_INSERT_SIZE = histogram.getMedian();
        metrics.MEDIAN_ABSOLUTE_DEVIATION = histogram.getMedianAbsoluteDeviation();

        setWidthValues(histogram, metrics);

        return metrics;
    }

    private static void setWidthValues(HistogramDataflow<Integer> histogram, InsertSizeMetrics metrics) {
        final double median = histogram.getMedian();
        double covered = 0;
        double low = median;
        double high = median;

        while (low >= histogram.getMin() || high <= histogram.getMax()) {
            @SuppressWarnings("unchecked")
            final Histogram<Integer>.Bin lowBin = histogram.get((int) low);

            if (lowBin != null) {
                covered += lowBin.getValue();
            }

            if (low != high) {
                @SuppressWarnings("unchecked")
                final Histogram<Integer>.Bin highBin = histogram.get((int) high);
                if (highBin != null) {
                    covered += highBin.getValue();
                }
            }

            final double percentCovered = covered / histogram.getCount();
            final int distance = (int) (high - low) + 1;
            updateWidthValues(metrics, percentCovered, distance);

            --low;
            ++high;
        }
    }

    private static void updateWidthValues(InsertSizeMetrics metrics, double percentCovered, int distance) {
        if (percentCovered >= 0.1 && metrics.WIDTH_OF_10_PERCENT == 0) {metrics.WIDTH_OF_10_PERCENT = distance;}
        if (percentCovered >= 0.2 && metrics.WIDTH_OF_20_PERCENT == 0) {metrics.WIDTH_OF_20_PERCENT = distance;}
        if (percentCovered >= 0.3 && metrics.WIDTH_OF_30_PERCENT == 0) {metrics.WIDTH_OF_30_PERCENT = distance;}
        if (percentCovered >= 0.4 && metrics.WIDTH_OF_40_PERCENT == 0) {metrics.WIDTH_OF_40_PERCENT = distance;}
        if (percentCovered >= 0.5 && metrics.WIDTH_OF_50_PERCENT == 0) {metrics.WIDTH_OF_50_PERCENT = distance;}
        if (percentCovered >= 0.6 && metrics.WIDTH_OF_60_PERCENT == 0) {metrics.WIDTH_OF_60_PERCENT = distance;}
        if (percentCovered >= 0.7 && metrics.WIDTH_OF_70_PERCENT == 0) {metrics.WIDTH_OF_70_PERCENT = distance;}
        if (percentCovered >= 0.8 && metrics.WIDTH_OF_80_PERCENT == 0) {metrics.WIDTH_OF_80_PERCENT = distance;}
        if (percentCovered >= 0.9 && metrics.WIDTH_OF_90_PERCENT == 0) {metrics.WIDTH_OF_90_PERCENT = distance;}
        if (percentCovered >= 0.99 && metrics.WIDTH_OF_99_PERCENT == 0) {metrics.WIDTH_OF_99_PERCENT = distance;}
    }

    /**
     * If histogramWidth is null infer a value for it
     */
    private static int inferHistogramWidth(Integer histogramWidth, InsertSizeMetrics metrics, double deviations) {
        if (histogramWidth == null) {
            return (int) (metrics.MEDIAN_INSERT_SIZE + (deviations * metrics.MEDIAN_ABSOLUTE_DEVIATION));
        } else {
            return histogramWidth;
        }
    }

    @Override
    public CombineInsertSizeHistogramsIntoMetricsFile createAccumulator() {
        return new CombineInsertSizeHistogramsIntoMetricsFile(this.DEVIATIONS, this.HISTOGRAM_WIDTH, this.MINIMUM_PERCENT);
    }

    @Override
    public CombineInsertSizeHistogramsIntoMetricsFile addInput(CombineInsertSizeHistogramsIntoMetricsFile accumulator, KV<InsertSizeAggregationLevel, HistogramDataflow<Integer>> input) {
        if(accumulator.histograms.containsKey(input.getKey())){
            throw new GATKException("Cannot merge two histograms with the same key. " +
                    "CombineHistogramsIntoMetrics file assumes that each key is represented only once. " +
                    "Key:" + input.getKey().toString());
        }
        accumulator.histograms.put(input.getKey(), input.getValue());
        return accumulator;
    }

    @Override
    public CombineInsertSizeHistogramsIntoMetricsFile mergeAccumulators(Iterable<CombineInsertSizeHistogramsIntoMetricsFile> accumulators) {
        Map<InsertSizeAggregationLevel, HistogramDataflow<Integer>> histograms = StreamSupport.stream(accumulators.spliterator(), false)
                .flatMap(a -> a.histograms.entrySet().stream())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

        CombineInsertSizeHistogramsIntoMetricsFile accum = createAccumulator();
        accum.histograms.putAll(histograms);
        return accum;
    }

    @Override
    public MetricsFileDataflow<InsertSizeMetrics, Integer> extractOutput(CombineInsertSizeHistogramsIntoMetricsFile accumulator) {
        final MetricsFileDataflow<InsertSizeMetrics, Integer> metricsFile = new MetricsFileDataflow<>();

        double totalInserts = accumulator.histograms.values()
                .stream()
                .mapToDouble(Histogram::getCount)
                .sum();

        Map<InsertSizeAggregationLevel, HistogramDataflow<Integer>> sortedHistograms = new TreeMap<>(Comparator.comparing((InsertSizeAggregationLevel a) -> a.getSample() != null ? a.getSample() : "")
                .thenComparing(a -> a.getLibrary() != null ? a.getLibrary() : "")
                .thenComparing(a -> a.getReadGroup() != null ? a.getReadGroup() : ""));
        sortedHistograms.putAll(accumulator.histograms);
        sortedHistograms.entrySet().forEach(kv -> addHistogramToMetricsFile(kv.getKey(), kv.getValue(), metricsFile, totalInserts));


        return metricsFile;
    }

}
