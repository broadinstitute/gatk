package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationKey;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.CopyNumberAnnotations;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents a collection of intervals annotated with {@link CopyNumberAnnotations}.
 * Supports {@link AnnotationKey}s of integer, long, double, and string type.
 * Can be constructed from a TSV file that contains the standard interval column headers,
 * any subset of the {@link CopyNumberAnnotations}, and additional columns (which are ignored).
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AnnotatedIntervalCollection extends AbstractLocatableCollection<LocatableMetadata, AnnotatedInterval> {
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END; columns headers for additional annotations can be specified
     */
    enum AnnotatedIntervalTableColumn {
        CONTIG,
        START,
        END;

        static final TableColumnCollection STANDARD_COLUMNS = new TableColumnCollection((Object[]) values());
    }

    enum AnnotationValueType {
        Integer,
        Long,
        Double,
        String
    }
    
    private static final BiConsumer<AnnotatedInterval, DataLine> ANNOTATED_INTERVAL_RECORD_TO_DATA_LINE_ENCODER = (annotatedInterval, dataLine) -> {
        dataLine.append(annotatedInterval.getInterval().getContig())
                .append(annotatedInterval.getInterval().getStart())
                .append(annotatedInterval.getInterval().getEnd());
        final AnnotationMap annotations = annotatedInterval.getAnnotationMap();
        for (final AnnotationKey<?> key : annotations.getKeys()) {
            final AnnotationValueType type = AnnotationValueType.valueOf(key.getType().getSimpleName());
            switch (type) {
                case Integer:
                    dataLine.append((Integer) annotations.getValue(key));
                    break;
                case Long:
                    dataLine.append((Long) annotations.getValue(key));
                    break;
                case Double:
                    dataLine.append(formatDouble((Double) annotations.getValue(key)));
                    break;
                case String:
                    dataLine.append((String) annotations.getValue(key));
                    break;
                default:
                    throw new UserException.BadInput(String.format("Unsupported annotation type: %s", type));
            }
        }
    };

    public AnnotatedIntervalCollection(final File inputFile) {
        this(inputFile, getAnnotationKeys(CopyNumberFormatsUtils.readColumnsFromHeader(inputFile)));
    }

    private AnnotatedIntervalCollection(final File inputFile,
                                        final List<AnnotationKey<?>> annotationKeys) {
        super(
                inputFile,
                getColumns(annotationKeys),
                getAnnotatedIntervalRecordFromDataLineDecoder(annotationKeys),
                ANNOTATED_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }

    public AnnotatedIntervalCollection(final LocatableMetadata metadata,
                                       final List<AnnotatedInterval> annotatedIntervals) {
        super(
                metadata,
                annotatedIntervals,
                getColumns(getAnnotationKeys(annotatedIntervals)),
                getAnnotatedIntervalRecordFromDataLineDecoder(getAnnotationKeys(annotatedIntervals)),
                ANNOTATED_INTERVAL_RECORD_TO_DATA_LINE_ENCODER);
    }

    private static TableColumnCollection getColumns(final List<AnnotationKey<?>> annotationKeys) {
        return new TableColumnCollection(
                ListUtils.union(
                        AnnotatedIntervalTableColumn.STANDARD_COLUMNS.names(),
                        annotationKeys.stream().map(AnnotationKey::getName).collect(Collectors.toList())));
    }

    private static List<AnnotationKey<?>> getAnnotationKeys(final TableColumnCollection columns) {
        Utils.nonNull(columns);
        Utils.validateArg(columns.columnCount() != 0, "TableColumnCollection cannot be empty.");
        Utils.validateArg(columns.containsAll(AnnotatedIntervalTableColumn.STANDARD_COLUMNS.names()),
                String.format("TableColumnCollection must contain standard columns: %s.",
                        AnnotatedIntervalTableColumn.STANDARD_COLUMNS.names()));
        return CopyNumberAnnotations.ANNOTATIONS.stream()
                .filter(a -> columns.contains(a.getName()))
                .collect(Collectors.toList());
    }

    private static List<AnnotationKey<?>> getAnnotationKeys(final List<AnnotatedInterval> annotatedIntervals) {
        return annotatedIntervals.isEmpty() ? new ArrayList<>() : annotatedIntervals.get(0).getAnnotationMap().getKeys();
    }

    private static Function<DataLine, AnnotatedInterval> getAnnotatedIntervalRecordFromDataLineDecoder(
            final List<AnnotationKey<?>> annotationKeys) {
        return dataLine -> {
            final String contig = dataLine.get(AnnotatedIntervalTableColumn.CONTIG);
            final int start = dataLine.getInt(AnnotatedIntervalTableColumn.START);
            final int end = dataLine.getInt(AnnotatedIntervalTableColumn.END);
            final SimpleInterval interval = new SimpleInterval(contig, start, end);
            final List<Pair<AnnotationKey<?>, Object>> annotations = new ArrayList<>(annotationKeys.size());
            for (final AnnotationKey<?> key : annotationKeys) {
                final AnnotationValueType type = AnnotationValueType.valueOf(key.getType().getSimpleName());
                switch (type) {
                    case Integer:
                        annotations.add(Pair.of(key, dataLine.getInt(key.getName())));
                        break;
                    case Long:
                        annotations.add(Pair.of(key, dataLine.getLong(key.getName())));
                        break;
                    case Double:
                        annotations.add(Pair.of(key, dataLine.getDouble(key.getName())));
                        break;
                    case String:
                        annotations.add(Pair.of(key, dataLine.get(key.getName())));
                        break;
                    default:
                        throw new UserException.BadInput(String.format("Unsupported annotation type: %s", type));
                }
            }
            final AnnotationMap annotationMap = new AnnotationMap(annotations);
            return new AnnotatedInterval(interval, annotationMap);
        };
    }

    /**
     * Columns, encoder, and decoder are not used.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final AbstractRecordCollection<?, ?> that = (AbstractRecordCollection<?, ?>) o;
        return getMetadata().equals(that.getMetadata()) &&
                getRecords().equals(that.getRecords());
    }

    /**
     * Columns, encoder, and decoder are not used.
     */
    @Override
    public int hashCode() {
        int result = getMetadata().hashCode();
        result = 31 * result + getRecords().hashCode();
        return result;
    }
}