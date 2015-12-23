package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

/**
 * Table writer for {@link TargetCoverageStats} records.
 * <p>
 *     The output format is a tab separated file with a header line given names to the different
 *     table columns:
 *     <dl>
 *        <dt>NAME</dt><dd>Target name (mandatory)</dd>.
 *        <dt>MEAN</dt><dd>Target coverage mean across samples (mandatory).</dd>
 *        <dt>VARIANCE</dt><dd>Target coverage (sample) variance across samples (mandatory).</dd>
 *        <dt>CONTIG</dt><dd>Target contig name (optional).</dd>
 *        <dt>START</dt><dd>Target start position (optional).</dd>
 *        <dt>END</dt><dd>Target end position (optional).</dd>
 *     </dl>
 * </p>
 * <p>
 *    Targets might be indicated by name only or by name and interval. If the former only the
 *    {@link TargetTableColumns#NAME NAME} column is present. If the latter also the coordinate
 *    columns are present: {@link TargetTableColumns#CONTIG CONTIG}, {@link TargetTableColumns#START START}
 *    and {@link TargetTableColumns#END END}.
 * </p>
 * <p>
 * Example without coordinates:
 * <pre>
 *      NAME     MEAN    VARIANCE
 *      target1  100.3   23.5
 *      target2  45.7    16.1
 *      ...
 * </pre>
 * </p>
 * <p>
 * Example with coordinates:
 * <pre>
 *      CONTIG  START   END     NAME    MEAN    VARIANCE
 *      chr1    1       100     target1 100.3   23.5
 *      chr1    202     301     target2 45.7    16.1
 *      ...
 * </pre>
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetCoverageStatsWriter extends TableWriter<TargetCoverageStats> {

    private final boolean outputIntervals;

    /**
     * Creates a new writer indicating whether the target's interval will be part of the output or not.
     * <p>
     *     At this point the caller must decide whether target intervals are to be included in the output.
     *     If this is the case, {@link #writeRecord} will fail with a {@link IllegalArgumentException} if
     *     the input record does not have an interval.
     * </p>
     * @param file the output file name.
     * @param outputIntervals {@code true} if the target intervals must be output, {@code false} otherwise.
     * @throws IOException if a {@link IOException} was thrown when writing into the output.
     */
    public TargetCoverageStatsWriter(final File file, final boolean outputIntervals) throws IOException {
        super(file, outputIntervals ?
                new TableColumnCollection(
                        TargetTableColumns.CONTIG,
                        TargetTableColumns.START,
                        TargetTableColumns.END,
                        TargetTableColumns.NAME,
                        TargetCoverageStats.MEAN_COLUMN_NAME,
                        TargetCoverageStats.VARIANCE_COLUMN_NAME)
               : new TableColumnCollection(
                        TargetTableColumns.NAME,
                        TargetCoverageStats.MEAN_COLUMN_NAME,
                        TargetCoverageStats.VARIANCE_COLUMN_NAME));
        this.outputIntervals = outputIntervals;
    }

    @Override
    protected void composeLine(final TargetCoverageStats record, final DataLine dataLine) {
        if (outputIntervals) {
            final SimpleInterval interval = record.target.getInterval();
            Utils.nonNull(interval, "the target must have an interval since intervals are output by this writer");
            dataLine.append(record.target.getContig())
                    .append(record.target.getStart())
                    .append(record.target.getEnd());
        }
        dataLine.append(record.target.getName())
                .append(record.mean)
                .append(record.variance);
    }
}
