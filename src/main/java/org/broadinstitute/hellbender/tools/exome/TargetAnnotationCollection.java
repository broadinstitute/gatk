package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Simple target annotation collection based on an EnumMap from the annotation to its value as an string.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetAnnotationCollection {

    private final EnumMap<TargetAnnotation, String> values = new EnumMap<TargetAnnotation, String>(TargetAnnotation.class);

    public TargetAnnotationCollection() { }

    public void put(final TargetAnnotation annotation, final String value) {
        values.put(annotation, value);
    }

    /**
     * Creates a new annotation collection based on a annotation to string map.
     * <p>
     *     Later changes in the input map won't have any effect on this collection.
     * </p>
     *
     * @param annotationMap map between annotations and values.
     */
    public TargetAnnotationCollection(final Map<TargetAnnotation, String> annotationMap) {
        Utils.nonNull(annotationMap);
        Utils.validateArg(!annotationMap.containsValue(null), "the input map cannot contain a null value");
        Utils.validateArg(!annotationMap.containsKey(null), "the input map cannot contain a null key");

        for (final Map.Entry<TargetAnnotation, String> entry : annotationMap.entrySet()) {
            put(entry.getKey(), entry.getValue());
        }
    }

    public int size() {
        return values.size();
    }

    public boolean hasAnnotation(final TargetAnnotation annotation) {
        return values.containsKey(annotation);
    }

    public double getDouble(final TargetAnnotation annotation) {
        final String stringValue = get(annotation);
        try {
            return Double.parseDouble(stringValue);
        } catch (final NumberFormatException ex) {
            throw new IllegalStateException("the annotation value is not a valid double string: " + stringValue);
        }
    }

    public String get(final TargetAnnotation annotation) {
        final String result = values.get(annotation);
        if (result == null) {
            throw new NoSuchElementException("there is no value for annotation " + annotation);
        } else {
            return result;
        }
    }

    public Set<TargetAnnotation> annotationSet() {
        return Collections.unmodifiableSet(values.keySet());
    }
}
