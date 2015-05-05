package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.insertsize;

import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.common.collect.Sets;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.MetricsFileDataflow;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

/**
 * Combiner to combine multiple MetricsFileDataflow at different aggregation levels into a single {@link MetricsFileDataflow}
 * for output.
 */
class CombineInsertSizeMetricsFiles
        extends Combine.CombineFn<MetricsFileDataflow<InsertSizeMetrics, Integer>, MetricsFileDataflow<InsertSizeMetrics, Integer>, MetricsFileDataflow<InsertSizeMetrics, Integer>> {
    public static final long serialVersionUID = 1l;

    @Override
    public MetricsFileDataflow<InsertSizeMetrics, Integer> createAccumulator() {
        return new MetricsFileDataflow<>();
    }

    @Override
    public MetricsFileDataflow<InsertSizeMetrics, Integer> addInput(MetricsFileDataflow<InsertSizeMetrics, Integer> accumulator, MetricsFileDataflow<InsertSizeMetrics, Integer> input) {
        return combineMetricsFiles(accumulator, input);
    }

    private MetricsFileDataflow<InsertSizeMetrics, Integer> combineMetricsFiles(MetricsFileDataflow<InsertSizeMetrics, Integer> accumulator, MetricsFileDataflow<InsertSizeMetrics, Integer> input) {
        Set<Header> headers = Sets.newLinkedHashSet(accumulator.getHeaders());
        Set<Header> inputHeaders = Sets.newLinkedHashSet(input.getHeaders());
        inputHeaders.removeAll(headers);
        inputHeaders.forEach(accumulator::addHeader);

        accumulator.addAllMetrics(input.getMetrics());
        input.getAllHistograms().forEach(accumulator::addHistogram);
        return accumulator;
    }

    @Override
    public MetricsFileDataflow<InsertSizeMetrics, Integer> mergeAccumulators(Iterable<MetricsFileDataflow<InsertSizeMetrics, Integer>> accumulators) {
        MetricsFileDataflow<InsertSizeMetrics, Integer> base = createAccumulator();
        accumulators.forEach(accum -> combineMetricsFiles(base, accum));
        return base;
    }

    @Override
    public MetricsFileDataflow<InsertSizeMetrics, Integer> extractOutput(MetricsFileDataflow<InsertSizeMetrics, Integer> accumulator) {
        List<InsertSizeMetrics> metrics = new ArrayList<>(accumulator.getMetrics());
        metrics.sort(InsertSizeMetrics.getInsertSizeMetricsComparator()); //sort the metrics
        MetricsFileDataflow<InsertSizeMetrics, Integer> sorted = new MetricsFileDataflow<>();
        sorted.addAllMetrics(metrics);
        accumulator.getAllHistograms().stream().sorted(Comparator.comparing(Histogram::getValueLabel))
                .forEach(sorted::addHistogram); //sort the histograms
        accumulator.getHeaders().stream().forEach(sorted::addHeader);
        return sorted;
    }
}
