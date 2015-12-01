package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * Collections of annotations for a given target.
 *
 * <p>
 *     Annotations are properties of a target not including its name and coordinates.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public interface TargetAnnotationCollection {

    /**
     * Number of annotations in the collection.
     * @return 0 or greater.
     */
    int size();

    /**
     * Checks whether an annotation is present.
     * @param annotation the target annotation.
     * @return {@code true} iff the corresponding target has a value for the input annotation.
     */
    boolean hasAnnotation(TargetAnnotation annotation);

    /**
     * Returns the value of the annotation as a double.
     * @param annotation the target annotation.
     * @return any double value.
     * @throws NoSuchElementException if the corresponding target does not have a value for such an annotation.
     * @throws UserException.BadInput if the annotation cannot be transformed into a double (the user has provided an
     *   annotation with the wrong format).
     */
    double getDouble(TargetAnnotation annotation);

    /**
     * Returns the annotation like an string.
     * @param annotation the target annotation.
     * @return never {@code null}.
     * @throws NoSuchElementException if the annotation is not present in the collection.
     */
    String get(final TargetAnnotation annotation);

    Set<TargetAnnotation> annotationSet();

    /**
     * Special annotation collection that represents the absence of all annotations.
     */
    TargetAnnotationCollection NO_ANNOTATION = new TargetAnnotationCollection() {

        @Override
        public int size() {
            return 0;
        }

        @Override
        public boolean hasAnnotation(final TargetAnnotation annotation) {
            return false;
        }

        @Override
        public double getDouble(final TargetAnnotation annotation) {
            Utils.nonNull(annotation, "the input annotation cannot be null");
            throw new NoSuchElementException("there is not such a annotation " + annotation);
        }

        @Override
        public String get(final TargetAnnotation annotation) {
            Utils.nonNull(annotation, "the input annotation cannot be null");
            throw new NoSuchElementException("there is not such a annotation " + annotation);
        }

        @Override
        public Set<TargetAnnotation> annotationSet() {
            return Collections.emptySet();
        }
    };

    static TargetAnnotationCollection fromMap(final Map<TargetAnnotation, String> annotationMap) {
        Utils.nonNull(annotationMap);
        final int size = annotationMap.size();
        if (size == 0) {
            return NO_ANNOTATION;
        } else {
            return new HashTargetAnnotationCollection(annotationMap);
        }
    }
}
