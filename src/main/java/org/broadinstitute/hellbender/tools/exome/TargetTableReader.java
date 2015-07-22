package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Target table file reader.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetTableReader extends TableReader<Target> {

    /**
     * Enum of the standard column names for a target list file.
     */
    public enum StandardColumn {
        CONTIG, START, END, NAME;

        /**
         * Column name array.
         */
        private static final String[] COLUMN_NAME_ARRAY =
                Stream.of(values()).map(StandardColumn::toString).toArray(String[]::new);

        /**
         * Column name set.
         */
        public static final Set<String> COLUMN_NAME_SET =
                Collections.unmodifiableSet(Stream.of(COLUMN_NAME_ARRAY).collect(Collectors.toSet()));
    }

    public TargetTableReader(final File file) throws IOException {
        super(file);
    }

    @Override
    protected void processColumns(final TableColumnCollection columns) {
        if (!columns.containsAll(StandardColumn.COLUMN_NAME_ARRAY)) {
            throw formatException("Bad header: missing required columns: "
                    + Stream.of(StandardColumn.COLUMN_NAME_ARRAY)
                            .filter(c -> !columns.contains(c)).collect(Collectors.joining(", ")));
        }
    }

    @Override
    protected Target createRecord(final DataLine dataLine) {

        final String contig = dataLine.get(StandardColumn.CONTIG);
        final int start = dataLine.getInt(StandardColumn.START);
        final int end = dataLine.getInt(StandardColumn.END);
        if (start < 0) {
            throw formatException("the start position must not be negative: " + start);
        } else if (start > end) {
            throw formatException(String.format("the start position %d cannot be greater than the end position %d", start, end));
        }
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        final String name = dataLine.get(StandardColumn.NAME);
        return new Target(name, interval);
    }
}
