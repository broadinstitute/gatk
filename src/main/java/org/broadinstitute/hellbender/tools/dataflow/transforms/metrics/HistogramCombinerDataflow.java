package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics;

import com.google.cloud.dataflow.sdk.transforms.Combine;

import java.io.Serializable;

/**
 * Use with {@link HistogramDataflow} to combine histograms.
 * @param <K>
 */
public class HistogramCombinerDataflow<K extends Comparable<K> & Serializable> extends Combine.AccumulatingCombineFn<K, HistogramDataflow<K>, HistogramDataflow<K>> {
    private static final long serialVersionUID = 1l;

    @Override
    public HistogramDataflow<K> createAccumulator() {
        return new HistogramDataflow<>();
    }
}
