package org.broadinstitute.hellbender.tools.exome.coveragestats;

import org.broadinstitute.hellbender.tools.exome.coveragestats.SampleCoverageStats;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Table writer for {@link SampleCoverageStats} objects.
 * <p>
 *     The file is formatted as a tab-separated three column table with the following columns:
 *     <dl>
 *      <dt>SAMPLE</dt><dd>The sample name.</dd>
 *      <dt>MEAN</dt><dd>The sample mean coverage across targets.</dd>
 *      <dt>VARIATION</dt><dd>The sample coverage population variance across targets.</dd>
 *     </dl>
 *     The order is irrelevant but all those columns must be present.
 * </p>
 *
 * <p>
 *     Example:
 *     <pre>
 *         SAMPLE   MEAN    VARIANCE
 *         sample_1 100.34  34.3
 *         sample_2 45.12   15.1
 *         ...
 *     </pre>
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class SampleCoverageStatsWriter extends TableWriter<SampleCoverageStats> {

    /**
     * Creates a new writer.
     *
     * @param file the destination file.
     * @throws IOException if such an exception is thrown when opening or writing into {@code file}.
     */
    public SampleCoverageStatsWriter(final File file) throws IOException {
        super(file, new TableColumnCollection(
                SampleCoverageStats.SAMPLE_COLUMN_NAME,
                SampleCoverageStats.MEAN_COLUMN_NAME,
                SampleCoverageStats.VARIANCE_COLUMN_NAME));
    }

    @Override
    protected void composeLine(final SampleCoverageStats record, final DataLine dataLine) {
        dataLine.append(record.sample)
                .append(record.mean)
                .append(record.variance);
    }
}
