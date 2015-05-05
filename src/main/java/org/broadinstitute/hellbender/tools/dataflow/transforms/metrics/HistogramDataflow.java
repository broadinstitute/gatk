package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics;

import com.google.cloud.dataflow.sdk.coders.DefaultCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import htsjdk.samtools.util.Histogram;

import java.io.Serializable;

/**
 * Mergeable, Serializable {@link Histogram} for use with dataflow.
 * @param <K> the type of the values being histogrammed
 */
@DefaultCoder(SerializableCoder.class)
public class HistogramDataflow<K extends Comparable<K> & Serializable> extends Histogram<K> implements Combine.AccumulatingCombineFn.Accumulator<K, HistogramDataflow<K>, HistogramDataflow<K>>{
    private final static long serialVersionUID = 1l;

    @Override
    public void addInput(K input) {
        this.increment(input);
    }

    @Override
    public void mergeAccumulator(HistogramDataflow<K> other) {
        this.addHistogram(other);
    }

    @Override
    public HistogramDataflow<K> extractOutput() {
        return this;
    }
}

