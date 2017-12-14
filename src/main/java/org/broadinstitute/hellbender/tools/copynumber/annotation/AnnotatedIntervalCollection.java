package org.broadinstitute.hellbender.tools.copynumber.annotation;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AnnotatedIntervalCollection extends LocatableCollection<AnnotatedInterval> {
    enum AnnotatedIntervalTableColumn {
        CONTIG,
        START,
        END,
        GC_CONTENT;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
    
    private static final Function<DataLine, AnnotatedInterval> ANNOTATED_INTERVAL_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(AnnotatedIntervalTableColumn.CONTIG);
        final int start = dataLine.getInt(AnnotatedIntervalTableColumn.START);
        final int end = dataLine.getInt(AnnotatedIntervalTableColumn.END);
        final double gcContent = dataLine.getDouble(AnnotatedIntervalTableColumn.GC_CONTENT);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        final AnnotationSet annotationSet = new AnnotationSet(gcContent);
        return new AnnotatedInterval(interval, annotationSet);
    };

    private static final BiConsumer<AnnotatedInterval, DataLine> ANNOTATED_INTERVAL_RECORD_TO_DATA_LINE_ENCODER = (annotatedInterval, dataLine) ->
            dataLine.append(annotatedInterval.getInterval().getContig())
                    .append(annotatedInterval.getInterval().getStart())
                    .append(annotatedInterval.getInterval().getEnd())
                    .append(annotatedInterval.getAnnotationSet().getGCContent());

    public AnnotatedIntervalCollection(final File inputFile) {
        super(inputFile, AnnotatedIntervalCollection.AnnotatedIntervalTableColumn.COLUMNS, ANNOTATED_INTERVAL_RECORD_FROM_DATA_LINE_DECODER, ANNOTATED_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }

    public AnnotatedIntervalCollection(final List<AnnotatedInterval> annotatedIntervals) {
        super(annotatedIntervals, AnnotatedIntervalCollection.AnnotatedIntervalTableColumn.COLUMNS, ANNOTATED_INTERVAL_RECORD_FROM_DATA_LINE_DECODER, ANNOTATED_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }
}
