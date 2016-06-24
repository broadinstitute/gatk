package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;

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
        TableUtils.checkMandatoryColumns(columns, SampleCoverageStats.COLUMNS, this::formatException);
        sampleColumnIndex = columns.indexOf(SampleCoverageStats.SAMPLE_COLUMN_NAME);
        meanColumnIndex = columns.indexOf(SampleCoverageStats.MEAN_COLUMN_NAME);
        varianceColumnIndex = columns.indexOf(SampleCoverageStats.VARIANCE_COLUMN_NAME);
    }

    @Override
    protected SampleCoverageStats createRecord(final DataLine dataLine) {
        return new SampleCoverageStats(dataLine.get(sampleColumnIndex), dataLine.getDouble(meanColumnIndex),
                dataLine.getDouble(varianceColumnIndex));
    }
}
