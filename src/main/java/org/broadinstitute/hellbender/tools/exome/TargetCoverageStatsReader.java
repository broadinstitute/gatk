package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

/**
 * Table reader for {@link TargetCoverageStats} records.
 * <p>
 *     The input file format must follow the one described in {@link TargetCoverageStatsWriter} documentation.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetCoverageStatsReader extends TableReader<TargetCoverageStats> {

    private final Function<DataLine, Target> targetExtractor;
    private final ToDoubleFunction<DataLine> meanExtractor;
    private final ToDoubleFunction<DataLine> varianceExtractor;
    private final ToDoubleFunction<DataLine> interquartileRangeExtractor;

    /**
     * Creates a new reader given the source file.
     *
     * @param outputFile the source file.
     * @throws IOException if such an exception is thrown when opening or reading from the source file.
     * @throws IllegalArgumentException if {@code outputFile} is {@code null}.
     */
    public TargetCoverageStatsReader(final File outputFile) throws IOException {
        super(outputFile);
        targetExtractor = createTargetExtractor(columns());
        meanExtractor = createStatExtractor(columns(), TargetCoverageStats.MEAN_COLUMN_NAME);
        varianceExtractor = createStatExtractor(columns(), TargetCoverageStats.VARIANCE_COLUMN_NAME);
        interquartileRangeExtractor = createStatExtractor(columns(), TargetCoverageStats.INTERQUARTILE_RANGE_COLUMN_NAME);
    }

    private ToDoubleFunction<DataLine> createStatExtractor(final TableColumnCollection columns, final String columnName) {
        final int columnIndex = columns.indexOf(columnName);
        if (columnIndex < 0) {
            throw formatException(String.format("missing mandatory column named '%s'", columnName));
        }
        return dataLine -> dataLine.getDouble(columnIndex);
    }

    private Function<DataLine,Target> createTargetExtractor(final TableColumnCollection columns) {
        if (columns.containsAll(TargetTableColumn.MANDATORY_COLUMNS.names())) {
            final int nameColumnIndex = columns.indexOf(TargetTableColumn.NAME.toString());
            final int contigColumnIndex = columns.indexOf(TargetTableColumn.CONTIG.toString());
            final int startColumnIndex = columns.indexOf(TargetTableColumn.START.toString());
            final int endColumnIndex = columns.indexOf(TargetTableColumn.END.toString());
            return dataLine ->
                    new Target(dataLine.get(nameColumnIndex),
                    new SimpleInterval(dataLine.get(contigColumnIndex),
                            dataLine.getInt(startColumnIndex),
                            dataLine.getInt(endColumnIndex)));

        } else {
            final int nameColumnIndex = columns.indexOf(TargetTableColumn.NAME.toString());
            if (nameColumnIndex == -1) {
                throw formatException(String.format("compulsory column name '%s' is missing",
                        TargetTableColumn.NAME.toString()));
            }
            return dataLine -> new Target(dataLine.get(nameColumnIndex));
        }
    }

    @Override
    protected TargetCoverageStats createRecord(final DataLine dataLine) {
        return new TargetCoverageStats(
                targetExtractor.apply(dataLine),
                meanExtractor.applyAsDouble(dataLine),
                varianceExtractor.applyAsDouble(dataLine),
                interquartileRangeExtractor.applyAsDouble(dataLine));
    }
}
