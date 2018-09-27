package org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation;

import com.google.common.collect.ImmutableMap;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Represents an immutable ordered collection of named, typed annotations for an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AnnotationMap {

    private final Map<AnnotationKey<?>, Object> annotationMap;

    public AnnotationMap(final List<Pair<AnnotationKey<?>, Object>> annotations) {
        Utils.nonEmpty(annotations);
        final ImmutableMap.Builder<AnnotationKey<?>, Object> builder = new ImmutableMap.Builder<>();
        annotations.forEach(a -> builder.put(a.getKey(), a.getValue()));
        annotationMap = builder.build();
    }

    public List<AnnotationKey<?>> getKeys() {
        return new ArrayList<>(annotationMap.keySet());
    }

    public <T> T getValue(final AnnotationKey<T> key) {
        Utils.nonNull(key);
        if (!annotationMap.containsKey(key)) {
            throw new IllegalArgumentException(
                    String.format("Annotation %s not contained in AnnotationMap.", key.getName()));
        }
        return key.getType().cast(annotationMap.get(key));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final AnnotationMap that = (AnnotationMap) o;
        return annotationMap.equals(that.annotationMap);
    }

    @Override
    public int hashCode() {
        return annotationMap.hashCode();
    }

    @Override
    public String toString() {
        return "AnnotationMap{" +
                "annotationMap=" + annotationMap +
                '}';
    }
}
