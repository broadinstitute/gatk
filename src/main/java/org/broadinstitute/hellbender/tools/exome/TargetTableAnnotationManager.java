package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.*;
import java.util.function.Function;

/**
 * Target table annotation collection.
 * <p>
 *     Holds information on what annotations are present in the target table file, what columns contain them and how to
 *     construct {@link TargetAnnotationCollection} instances to query these per target.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
final class TargetTableAnnotationManager {

    private Set<TargetAnnotation> annotationSet;

    private final int[] annotationToColumnIndex;

    private final Function<DataLine, Function<String, RuntimeException>> formatAnnotationExceptionFactory;

    /**
     * Creates a new target-table annotation manager.
     *
     * @param source name of the data-source that contains the target table (e.g. the enclosing file name).
     * @param columnCollection the table column collection for the target table.
     */
    public TargetTableAnnotationManager(final String source, final TableColumnCollection columnCollection) {
        Utils.nonNull(columnCollection);
        if (source == null) {
            formatAnnotationExceptionFactory = (dataLine) -> (message) -> new UserException.BadInput(String.format("in target table line %d: %s", dataLine.getLineNumber(), message));
        } else {
            formatAnnotationExceptionFactory = (dataLine) -> (message) -> new UserException.BadInput(String.format("in target table '%s', line %d: %s", source, dataLine.getLineNumber(), message));
        }

        final EnumSet<TargetAnnotation> annotationsPresent = EnumSet.noneOf(TargetAnnotation.class);
        annotationToColumnIndex = new int[TargetAnnotation.values().length];
        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            if ((annotationToColumnIndex[annotation.ordinal()] = columnCollection.indexOf(annotation.column.toString())) != -1) {
                annotationsPresent.add(annotation);
            }
        }
        annotationSet = Collections.unmodifiableSet(annotationsPresent);
    }

    /**
     * Creates a target annotation collection based on a data-line in the target table source.
     * @param dataLine the target data-line.
     * @return never {@code null}.
     */
    public TargetAnnotationCollection createTargetAnnotationCollection(final DataLine dataLine) {

        if (annotationSet.size() == 0) {
            return TargetAnnotationCollection.NO_ANNOTATION;
        } else {
            return new DataLineTargetAnnotationCollection(dataLine);
        }
    }

    /**
     * Target annotation collection for a data-line.
     */
    private class DataLineTargetAnnotationCollection implements TargetAnnotationCollection {
        private final DataLine dataLine;

        private DataLineTargetAnnotationCollection(DataLine dataLine) {
            this.dataLine = dataLine;
        }

        @Override
        public int size() {
            return annotationSet.size();
        }

        @Override
        public boolean hasAnnotation(final TargetAnnotation annotation) {
            return annotationToColumnIndex[annotation.ordinal()] >= 0;
        }

        @Override
        public double getDouble(final TargetAnnotation annotation) {
            if (hasAnnotation(annotation)) {
                return dataLine.getDouble(annotationToColumnIndex[annotation.ordinal()], formatAnnotationExceptionFactory.apply(dataLine));
            } else {
                throw formatAnnotationExceptionFactory.apply(dataLine).apply(
                    String.format("not found expected target annotation '%s' in input", annotation.column));
            }
        }

        @Override
        public String get(final TargetAnnotation annotation) {
            if (hasAnnotation(annotation)) {
                return dataLine.get(annotationToColumnIndex[annotation.ordinal()]);
            } else {
                throw formatAnnotationExceptionFactory.apply(dataLine).apply(
                        String.format("not found expected target annotation '%s' in input", annotation.column));
            }
        }

        @Override
        public Set<TargetAnnotation> annotationSet() {
            return annotationSet;
        }
    }
}
