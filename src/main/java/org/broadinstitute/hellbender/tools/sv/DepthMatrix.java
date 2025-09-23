package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

public final class DepthMatrix {
    private final List<SimpleInterval> bins;
    private final Map<String, double[]> matrix;

    public DepthMatrix() {
        this.bins = new ArrayList<>();
        this.matrix = new HashMap<>();
    }

    public DepthMatrix(List<SimpleInterval> bins, Map<String, double[]> matrix) {
        Utils.nonNull(bins);
        Utils.nonNull(matrix);
        for (final double[] value : matrix.values()) {
            Utils.nonNull(value);
            Utils.validateArg(value.length == bins.size(), "Each vector must have " + value.length + " bins");
        }
        this.bins = bins;
        this.matrix = matrix;
    }

    public double[] slice(final int binIndex, final Collection<String> samples) {
        Utils.validateArg(matrix.keySet().containsAll(samples), "Matrix does not contain all queried samples");
        if (matrix.isEmpty()) {
            return new double[0];
        }
        final double[] slice = new double[samples.size()];
        int j = 0;
        for (final String sample : samples) {
            slice[j++] = getSample(sample)[binIndex];
        }
        return slice;
    }

    public List<SimpleInterval> getBins() {
        return bins;
    }

    public double[] getSample(final String sample) {
        Utils.validateArg(matrix.containsKey(sample), "Matrix does not contain queried sample");
        return matrix.get(sample);
    }

    public boolean isEmpty() {
        return matrix.isEmpty();
    }

    public Set<String> getSampleSet() {
        return matrix.keySet();
    }

    public int getNumBins() {
        return bins.size();
    }
}
