package org.broadinstitute.hellbender.utils.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Target table file reader.
 *
 * <p>
 * <a href="#format"/>
 * The minimum content of a target table file consists of 4 columns indicate the target's name, that must be unique across
 * targets and its coordinates: enclosing contig, start and end positions. The start and end positions are inclusive
 * and 1-based (the first position index is 1).
 * </p>
 * <p>
 * Example:
 * <pre>
 *         CONTIG   START   END     NAME
 *         1        100     200     target_0
 *         1        300     400     target_1
 *         2        100     500     target_2
 * </pre>
 * </p>
 * <p>
 * There might be additional columns indicating arbitrary target annotations such as per-sample coverage, gc content
 * and so forth.
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetTableReader extends TableReader<Target> {

    private TargetTableAnnotationManager annotationCollection;

    public TargetTableReader(final Reader reader) throws IOException {
        super(reader);
    }

    public TargetTableReader(final File file) throws IOException {
        super(file);
    }

    @Override
    protected void processColumns(final TableColumnCollection columns) {
        if (!columns.containsAll(TargetTableColumns.MANDATORY_COLUMN_NAME_ARRAY)) {
            throw formatException("Bad header: missing required columns: "
                    + Stream.of(TargetTableColumns.MANDATORY_COLUMN_NAME_ARRAY)
                            .filter(c -> !columns.contains(c)).collect(Collectors.joining(", ")));
        }
        annotationCollection = new TargetTableAnnotationManager(getSource(), columns);
    }

    @Override
    protected Target createRecord(final DataLine dataLine) {

        final String contig = dataLine.get(TargetTableColumns.CONTIG);
        final int start = dataLine.getInt(TargetTableColumns.START);
        final int end = dataLine.getInt(TargetTableColumns.END);
        if (start < 0) {
            throw formatException("the start position must not be negative: " + start);
        } else if (start > end) {
            throw formatException(String.format("the start position %d cannot be greater than the end position %d", start, end));
        }
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        final String name = dataLine.get(TargetTableColumns.NAME);
        return new Target(name, interval, annotationCollection.createTargetAnnotationCollection(dataLine));
    }
}