package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval;

import com.google.common.collect.ImmutableSortedMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.SortedMap;

/**
 * Simple class that just has an interval and sorted name-value pairs.
 */
public final class AnnotatedInterval implements Locatable, Feature {

    private final SimpleInterval interval;
    private final SortedMap<String, String> annotations;

    public AnnotatedInterval(final SimpleInterval interval, final SortedMap<String, String> annotations) {
        this.interval = interval;
        this.annotations = annotations;
    }

    /** Returns a copy */
    public SimpleInterval getInterval() {
        return new SimpleInterval(this.interval);
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    /**
     * @param annotationName annotation name to query.
     * @return the value associated with the annotation name.  This can be {@code null}, if the annotation name is not
     * found.
     */
    public String getAnnotationValue(final String annotationName) {
        return annotations.get(annotationName);
    }

    /**
     * @param annotationName annotation name to query.
     * @param defaultValue Value to return if the annotation name is not found in this annotated interval.
     * @return the value associated with the annotation name or the specified default value if the annotation name is
     *  not found.  This can be {@code null}, if the default value is null.
     */
    public String getAnnotationValueOrDefault(final String annotationName, final String defaultValue) {
        return annotations.getOrDefault(annotationName, defaultValue);
    }

    /** Returns a copy of the annotations as a map.
     * Dev note: this does not create a copy, unless necessary.  See {@link ImmutableSortedMap#copyOfSorted(SortedMap)}*/
    public ImmutableSortedMap<String, String> getAnnotations() {
        return ImmutableSortedMap.copyOfSorted(this.annotations);
    }

    /**
     * @param annotationName Name of annotation to query.
     * @return Whether the given annotation is in this annotated interval.
     */
    public boolean hasAnnotation(final String annotationName) {
        return annotations.containsKey(annotationName);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final AnnotatedInterval that = (AnnotatedInterval) o;
        return this.interval.equals(that.getInterval()) && this.getAnnotations().equals(that.getAnnotations());
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + annotations.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "AnnotatedInterval{" +
                "interval=" + interval +
                ", annotations=" + annotations +
                '}';
    }
}
