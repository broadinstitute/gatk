package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

public class FilterStats {
    public static final String THRESHOLD_METADATA_TAG = "threshold";
    public static final String SENSITIVITY_METADATA_TAG = "sensitivity";
    public static final String FDR_METADATA_TAG = "fdr";

    private final String filterName;
    private final double falsePositiveCount;
    private final double falseDiscoveryRate;
    private final double falseNegativeCount;
    private final double falseNegativeRate;

    public FilterStats(final String filterName, final double falsePositiveCount, final double falseDiscoveryRate,
                       final double falseNegativeCount, final double falseNegativeRate){
        this.filterName = filterName;
        this.falsePositiveCount = falsePositiveCount;
        this.falseDiscoveryRate = falseDiscoveryRate;
        this.falseNegativeCount = falseNegativeCount;
        this.falseNegativeRate = falseNegativeRate;
    }

    public String getFilterName() { return filterName; }

    public double getFalsePositiveCount() { return falsePositiveCount; }

    public double getFalseDiscoveryRate() { return falseDiscoveryRate; }

    public double getFalseNegativeCount() { return falseNegativeCount; }

    public double getFalseNegativeRate() { return falseNegativeRate; }

    private enum M2FilterStatsTableColumn {
        FILTER("filter"),
        FALSE_POSITIVE_COUNT("FP"),
        FALSE_DISCOVERY_RATE("FDR"),
        FALSE_NEGATIVE_COUNT("FN"),
        FALSE_NEGATIVE_RATE("FNR");

        private String columnName;

        M2FilterStatsTableColumn(final String columnName) {
            this.columnName = columnName;
        }

        @Override
        public String toString() { return columnName; }

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static class Mutect2FilterStatsWriter extends TableWriter<FilterStats> {
        private Mutect2FilterStatsWriter(final File output) throws IOException {
            super(output, M2FilterStatsTableColumn.COLUMNS);
        }

        @Override
        protected void composeLine(final FilterStats stats, final DataLine dataLine) {
            dataLine.set(M2FilterStatsTableColumn.FILTER.toString(), stats.getFilterName())
                    .set(M2FilterStatsTableColumn.FALSE_POSITIVE_COUNT.toString(), stats.getFalsePositiveCount(), 2)
                    .set(M2FilterStatsTableColumn.FALSE_DISCOVERY_RATE.toString(), stats.getFalseDiscoveryRate(), 2)
                    .set(M2FilterStatsTableColumn.FALSE_NEGATIVE_COUNT.toString(), stats.getFalseNegativeCount(), 2)
                    .set(M2FilterStatsTableColumn.FALSE_NEGATIVE_RATE.toString(), stats.getFalseNegativeRate(), 2);
        }
    }

    public static void writeM2FilterSummary(final Collection<FilterStats> filterStats, final File outputTable, List<Pair<String, String>> clusteringMetadata,
                                            final double threshold, final double totalCalls, final double expectedTruePositives,
                                            final double expectedFalsePositives, final double expectedFalseNegatives) {
        try (Mutect2FilterStatsWriter writer = new Mutect2FilterStatsWriter(outputTable)) {
            for (final Pair<String, String> pair : clusteringMetadata) {
                writer.writeMetadata(pair.getKey(), pair.getValue());
            }
            writer.writeMetadata(THRESHOLD_METADATA_TAG, Double.toString(round(threshold)));
            writer.writeMetadata(FDR_METADATA_TAG, Double.toString(round(expectedFalsePositives / totalCalls)));
            writer.writeMetadata(SENSITIVITY_METADATA_TAG, Double.toString(round(expectedTruePositives / (expectedTruePositives + expectedFalseNegatives))));
            writer.writeAllRecords(filterStats);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", outputTable), e);
        }
    }

    private static double round(final double x) {
        return MathUtils.roundToNDecimalPlaces(x, 3);
    }

}
