package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

    public List<SimpleInterval> getBins() {
        return bins;
    }

    public Map<String, double[]> getMatrix() {
        return matrix;
    }

    public int getNumBins() {
        return bins.size();
    }
}
