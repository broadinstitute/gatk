package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import java.nio.file.Path;
import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.Pair;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Helper class used on the final pass of {@link FilterMutectCalls} to record total expected true positives, false positives,
 * and false negatives, as well as false positives and false negatives attributable to each filter
 */
public class FilteringOutputStats {
    private int pass = 0;
    private double TPs = 0;
    private double FPs = 0;
    private double FNs = 0;

    private Map<Mutect2Filter, MutableDouble> filterFPs;
    private Map<Mutect2Filter, MutableDouble> filterFNs;

    private final List<Mutect2Filter> filters;

    public FilteringOutputStats(final List<Mutect2Filter> filters) {
        this.filters = filters;
        filterFPs = makeEmptyFilterCounts();
        filterFNs = makeEmptyFilterCounts();
    }

    public void recordCall(final ErrorProbabilities errorProbabilities, final double threshold) {
        final double errorProbability = errorProbabilities.getCombinedErrorProbabilities();
        final boolean filtered = errorProbability > threshold;

        if (filtered) {
            FNs += 1 - errorProbability;
        } else {
            pass++;
            FPs += errorProbability;
            TPs += 1 - errorProbability;
        }

        for (final Map.Entry<Mutect2Filter, List<Double>> entry : errorProbabilities.getAlleleProbabilitiesByFilter().entrySet()) {
            final List<Double> filterArtifactProbability = entry.getValue();
            if (filterArtifactProbability > Mutect2FilteringEngine.EPSILON && filterArtifactProbability > threshold - Mutect2FilteringEngine.EPSILON) {
                filterFNs.get(entry.getKey()).add(1 - errorProbability);
            } else if (!filtered) {
                filterFPs.get(entry.getKey()).add(filterArtifactProbability);
            }
        }
    }

    public void writeFilteringStats(final Path filteringStats, final double threshold, List<Pair<String, String>> clusteringMetadata) {
        final double totalTrueVariants = TPs + FNs;

        final List<FilterStats> filterStats = filters.stream()
                .map(f -> new FilterStats(f.filterName(), filterFPs.get(f).getValue(), filterFPs.get(f).getValue() / pass,
                        filterFNs.get(f).getValue(), filterFNs.get(f).getValue() / totalTrueVariants))
                .filter(stats -> stats.getFalsePositiveCount() > 0 || stats.getFalseNegativeCount() > 0)
                .collect(Collectors.toList());

        FilterStats.writeM2FilterSummary(filterStats, filteringStats, clusteringMetadata, threshold, pass, TPs, FPs, FNs);
    }

    private Map<Mutect2Filter, MutableDouble> makeEmptyFilterCounts() {
        return filters.stream().collect(Collectors.toMap(f -> f, f -> new MutableDouble(0)));
    }

    public void clear() {
        pass = 0;
        TPs = 0;
        FPs = 0;
        FNs = 0;

        filterFPs = makeEmptyFilterCounts();
        filterFNs = makeEmptyFilterCounts();
    }
}
