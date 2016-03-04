package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.*;

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

    private final EnumMap<TargetAnnotation, Integer> annotationToColumnIndex = new EnumMap<TargetAnnotation, Integer>(TargetAnnotation.class);
    private final String source;

    /**
     * Creates a new target-table annotation manager.
     *
     * @param source name of the data-source that contains the target table (e.g. the enclosing file name).
     * @param columnCollection the table column collection for the target table.
     */
    public TargetTableAnnotationManager(final String source, final TableColumnCollection columnCollection) {
        Utils.nonNull(columnCollection);
        this.source = source;

        for (final TargetAnnotation annotation : TargetAnnotation.values()) {
            final int columnIndex = columnCollection.indexOf(annotation.column.toString());
            if (columnIndex != -1) {
                annotationToColumnIndex.put(annotation, columnIndex);
            }
        }
    }

    /**
     * Creates a target annotation collection based on a data-line in the target table source.
     * @param dataLine the target data-line.
     * @return never {@code null}.
     */
    public TargetAnnotationCollection createTargetAnnotationCollection(final DataLine dataLine) {
        final TargetAnnotationCollection annotationCollection = new TargetAnnotationCollection();

        for (final Map.Entry<TargetAnnotation, Integer> entry : annotationToColumnIndex.entrySet()) {
            final TargetAnnotation annotation = entry.getKey();
            final int columnIndex = entry.getValue();
            try {
                annotationCollection.put(annotation, dataLine.get(columnIndex));
            } catch (final IllegalStateException e) {
                new UserException.BadInput(String.format("Annotation %s not found in line %d oftarget table file %s.", annotation.toString(), dataLine.getLineNumber(), source));
            }
        }
        return annotationCollection;
    }

}
