package org.broadinstitute.hellbender.tools.copynumber.annotation;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * Represents a set of annotations for an interval.  Currently, only GC content is represented.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AnnotationSet {
    /**
     * If additional annotation fields are added here, then {@link AnnotatedIntervalCollection}
     * should be updated accordingly.
     */
    private final double gcContent;

    public AnnotationSet(final double gcContent) {
        Utils.validateArg((0. <= gcContent && gcContent <= 1.) || Double.isNaN(gcContent),
                "GC content must be in [0, 1] or NaN.");
        this.gcContent = gcContent;
    }

    public double getGCContent() {
        return gcContent;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final AnnotationSet that = (AnnotationSet) o;
        return Double.compare(that.gcContent, gcContent) == 0;
    }

    @Override
    public int hashCode() {
        long temp = Double.doubleToLongBits(gcContent);
        return (int) (temp ^ (temp >>> 32));
    }

    @Override
    public String toString() {
        return "AnnotationSet{" +
                "gcContent=" + gcContent +
                '}';
    }
}
