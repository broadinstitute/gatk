package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.annotation.AnnotationMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents an interval with a set of annotations.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AnnotatedInterval implements Locatable {
    private final SimpleInterval interval;
    private final AnnotationMap annotationMap;

    public AnnotatedInterval(final SimpleInterval interval,
                             final AnnotationMap annotationMap) {
        this.interval = Utils.nonNull(interval);
        this.annotationMap = Utils.nonNull(annotationMap);
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

    public SimpleInterval getInterval() {
        return interval;
    }

    public AnnotationMap getAnnotationMap() {
        return annotationMap;
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
        return interval.equals(that.interval) && annotationMap.equals(that.annotationMap);
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + annotationMap.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "AnnotatedInterval{" +
                "interval=" + interval +
                ", annotationMap=" + annotationMap +
                '}';
    }
}
