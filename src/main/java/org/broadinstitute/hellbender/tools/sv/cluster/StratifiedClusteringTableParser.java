package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.ImmutableSet;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.Set;
import java.util.function.Function;

public class StratifiedClusteringTableParser {

    // Configuration table column names
    public static final String NAME_COLUMN = "NAME";
    public static final String RECIPROCAL_OVERLAP_COLUMN = "RECIPROCAL_OVERLAP";
    public static final String SIZE_SIMILARITY_COLUMN = "SIZE_SIMILARITY";
    public static final String BREAKEND_WINDOW_COLUMN = "BREAKEND_WINDOW";
    public static final String SAMPLE_OVERLAP_COLUMN = "SAMPLE_OVERLAP";
    protected static final Set<String> COLUMN_NAMES = ImmutableSet.of(NAME_COLUMN, RECIPROCAL_OVERLAP_COLUMN, SIZE_SIMILARITY_COLUMN, BREAKEND_WINDOW_COLUMN, SAMPLE_OVERLAP_COLUMN);

    public static Function<DataLine, StratumParameters> tableParser(TableColumnCollection columns, Function<String, RuntimeException> exceptionFactory) {
        for (final String column : COLUMN_NAMES) {
            if (!columns.contains(column)) {
                throw exceptionFactory.apply("Missing column " + column);
            }
        }
        if (columns.columnCount() != COLUMN_NAMES.size()) {
            throw exceptionFactory.apply("Expected " + columns.columnCount() + " columns but found " + columns.columnCount());
        }
        return StratifiedClusteringTableParser::parseTableLine;
    }

    protected static StratumParameters parseTableLine(final DataLine dataLine) {
        final String name = dataLine.get(NAME_COLUMN);
        final double reciprocalOverlap = dataLine.getDouble(RECIPROCAL_OVERLAP_COLUMN);
        final double sizeSimilarity = dataLine.getDouble(SIZE_SIMILARITY_COLUMN);
        final double sampleOverlap = dataLine.getDouble(SAMPLE_OVERLAP_COLUMN);
        final int breakendWindow = dataLine.getInt(BREAKEND_WINDOW_COLUMN);
        return new StratumParameters(name, reciprocalOverlap, sizeSimilarity, breakendWindow, sampleOverlap);
    }

    public record StratumParameters(String name, double reciprocalOverlap, double sizeSimilarity,
                                    int breakendWindow, double sampleOverlap) {
    }
}
