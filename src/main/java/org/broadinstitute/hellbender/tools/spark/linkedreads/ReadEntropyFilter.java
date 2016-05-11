package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Map;
import java.util.stream.Collectors;

public class ReadEntropyFilter implements ReadFilter {
    private static final long serialVersionUID = 1l;

    public double minEntropy = 4.5;

    public ReadEntropyFilter(final double minEntropy) {
        this.minEntropy = minEntropy;
    }

    public static double computeEntropy(String sequence) {
        final Map<SVKmer, Long> frequencies = SVKmerizer.stream(sequence, 3).collect(Collectors.groupingBy(k -> k, Collectors.counting()));
        final int numKmers = sequence.length() - 2;
        final Double entropy = -1 * frequencies.entrySet().stream().collect(Collectors.summingDouble(e -> {
            final double p = (double) e.getValue() / numKmers;
            return p * Math.log(p) / Math.log(2.0d);
        }));
        return entropy;
    }

    @Override
    public boolean test(final GATKRead read) {
        return computeEntropy(read.getBasesString()) > minEntropy;
    }
}
