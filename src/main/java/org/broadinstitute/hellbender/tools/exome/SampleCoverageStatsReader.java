package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Tab-separated table reader for {@link SampleCoverageStats} instances.
 *
 * <p>
 * The expected file format is described in {@link SampleCoverageStatsWriter}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleCoverageStatsReader extends TableReader<SampleCoverageStats> {

    private final int sampleColumnIndex;

    private final int meanColumnIndex;

    private final int varianceColumnIndex;

    /**
     * Creates a new reader.
     *
     * @param file the source file.
     * @throws IOException if such a exception is throw when opening or reading from {@code file}.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     */
    public SampleCoverageStatsReader(final File file) throws IOException {
        super(file);
        final TableColumnCollection columns = columns();
        checkColumnsArePresent(columns);
        sampleColumnIndex = columns.indexOf(SampleCoverageStats.SAMPLE_COLUMN_NAME);
        meanColumnIndex = columns.indexOf(SampleCoverageStats.MEAN_COLUMN_NAME);
        varianceColumnIndex = columns.indexOf(SampleCoverageStats.VARIANCE_COLUMN_NAME);
    }

    private void checkColumnsArePresent(TableColumnCollection columns) {
        if (!columns.containsAll(SampleCoverageStats.SAMPLE_COLUMN_NAME, SampleCoverageStats.MEAN_COLUMN_NAME,
                SampleCoverageStats.VARIANCE_COLUMN_NAME)) {
            throw formatException(String.format("missing mandatory columns: %s",
                    Stream.of(SampleCoverageStats.SAMPLE_COLUMN_NAME, SampleCoverageStats.MEAN_COLUMN_NAME,
                            SampleCoverageStats.VARIANCE_COLUMN_NAME)
                            .filter(n -> !columns.contains(n))
                            .collect(Collectors.joining(", "))
                    ));
        }
    }

    @Override
    protected SampleCoverageStats createRecord(final DataLine dataLine) {
        return new SampleCoverageStats(dataLine.get(sampleColumnIndex), dataLine.getDouble(meanColumnIndex),
                dataLine.getDouble(varianceColumnIndex));
    }
}
