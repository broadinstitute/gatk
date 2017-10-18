package org.broadinstitute.hellbender.tools.copynumber.annotation;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents an interval with a set of annotations.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AnnotatedInterval implements Locatable, Feature {
    private final SimpleInterval interval;
    private final AnnotationSet annotationSet;

    public AnnotatedInterval(final SimpleInterval interval,
                             final AnnotationSet annotationSet) {
        this.interval = Utils.nonNull(interval);
        this.annotationSet = Utils.nonNull(annotationSet);
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

    public AnnotationSet getAnnotationSet() {
        return annotationSet;
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
        return interval.equals(that.interval) && annotationSet.equals(that.annotationSet);
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + annotationSet.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "AnnotatedInterval{" +
                "interval=" + interval +
                ", annotationSet=" + annotationSet +
                '}';
    }
}
