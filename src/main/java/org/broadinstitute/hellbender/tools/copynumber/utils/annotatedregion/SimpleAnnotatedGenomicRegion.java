package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.List;
import java.util.SortedMap;
import java.util.stream.Collectors;

/**
 * Simple class that just has an interval and name value pairs.
 *
 * When reading a TSV file of simple annotated genomic regions, the genomic region columns for the header when reading
 *  are specified in {@link SimpleAnnotatedGenomicRegion::CONTIG_HEADER},
 * {@link SimpleAnnotatedGenomicRegion::START_HEADER}, and {@link SimpleAnnotatedGenomicRegion::END_HEADER}
 */
public final class SimpleAnnotatedGenomicRegion implements Locatable {
    public final static String CONTIG_HEADER = "CONTIG";
    public final static String START_HEADER = "START";
    public final static String END_HEADER = "END";

    private SimpleInterval interval;
    private final SortedMap<String, String> annotations;

    public SimpleAnnotatedGenomicRegion(final SimpleInterval interval, final SortedMap<String, String> annotations) {
        this.interval = interval;
        this.annotations = annotations;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public SortedMap<String, String> getAnnotations() {
        return annotations;
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

    public String getAnnotationValue(final String annotationName) {
        return annotations.get(annotationName);
    }

    public boolean hasAnnotation(final String annotationName) {
        return annotations.containsKey(annotationName);
    }

    /**
     * Creates the annotation if it does not exist.
     *
     * @param annotationName the name for the annotation
     * @param annotationValue the value
     */
    public void setAnnotation(final String annotationName, final String annotationValue) {
        annotations.put(annotationName, annotationValue);
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final SimpleAnnotatedGenomicRegion that = (SimpleAnnotatedGenomicRegion) o;
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
        return interval.toString() + " :: " + annotations.entrySet().stream()
                .map(e -> e.getKey() + "->" + e.getValue()).collect(Collectors.joining(","));
    }

    /**
     * @return list of all of the annotations (not the values) in this annotated region.
     */
    public List<String> getAnnotationNames() {
        return getAnnotations().keySet().stream().sorted().collect(Collectors.toList());
    }
}
