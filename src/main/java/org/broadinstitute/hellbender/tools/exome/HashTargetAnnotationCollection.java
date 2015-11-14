package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Simple target annotation collection based on a hash map from the annotation to its value as an string.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
final class HashTargetAnnotationCollection implements TargetAnnotationCollection {

    private final Map<TargetAnnotation, String> values;

    /**
     * Creates a new annotation collection based on a annotation to string map.
     * <p>
     *     Later changes in the input map won't have any effect on this collection.
     * </p>
     *
     * @param annotationMap map between annotations an values.
     */
    public HashTargetAnnotationCollection(final Map<TargetAnnotation, String> annotationMap) {
        Utils.nonNull(annotationMap);
        if (annotationMap.containsValue(null)) {
            throw new IllegalArgumentException("the input map cannot contain a null value");
        } else if (annotationMap.containsKey(null)) {
            throw new IllegalArgumentException("the input map cannot contain a null key");
        } else if (annotationMap.isEmpty())  {
            values = Collections.emptyMap();
        } else {
            values = new EnumMap<>(annotationMap);
        }
    }

    @Override
    public int size() {
        return values.size();
    }

    @Override
    public boolean hasAnnotation(final TargetAnnotation annotation) {
        return values.containsKey(annotation);
    }

    @Override
    public double getDouble(final TargetAnnotation annotation) {
        final String stringValue = get(annotation);
        try {
            return Double.parseDouble(stringValue);
        } catch (final NumberFormatException ex) {
            throw new IllegalStateException("the annotation value is not a valid double string: " + stringValue);
        }
    }

    @Override
    public String get(final TargetAnnotation annotation) {
        final String result = values.get(annotation);
        if (result == null) {
            throw new NoSuchElementException("there is no value for annotation " + annotation);
        } else {
            return result;
        }
    }

    @Override
    public Set<TargetAnnotation> annotationSet() {
        return Collections.unmodifiableSet(values.keySet());
    }
}
