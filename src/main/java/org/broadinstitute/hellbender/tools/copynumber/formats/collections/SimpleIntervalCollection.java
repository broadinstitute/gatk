package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents a collection of {@link SimpleInterval} associated to a sample.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SimpleIntervalCollection extends AbstractLocatableCollection<LocatableMetadata, SimpleInterval> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END
     */
    enum SimpleIntervalTableColumn {
        CONTIG,
        START,
        END;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, SimpleInterval> SIMPLE_INTERVAL_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(SimpleIntervalTableColumn.CONTIG);
        final int start = dataLine.getInt(SimpleIntervalTableColumn.START);
        final int end = dataLine.getInt(SimpleIntervalTableColumn.END);
        return new SimpleInterval(contig, start, end);
    };

    private static final BiConsumer<SimpleInterval, DataLine> SIMPLE_INTERVAL_RECORD_TO_DATA_LINE_ENCODER = (simpleInterval, dataLine) ->
            dataLine.append(simpleInterval.getContig())
                    .append(simpleInterval.getStart())
                    .append(simpleInterval.getEnd());

    public SimpleIntervalCollection(final File inputFile) {
        super(inputFile, SimpleIntervalTableColumn.COLUMNS, SIMPLE_INTERVAL_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }

    public SimpleIntervalCollection(final LocatableMetadata metadata,
                                    final List<SimpleInterval> simpleIntervals) {
        super(metadata, simpleIntervals, SimpleIntervalTableColumn.COLUMNS, SIMPLE_INTERVAL_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }
}
