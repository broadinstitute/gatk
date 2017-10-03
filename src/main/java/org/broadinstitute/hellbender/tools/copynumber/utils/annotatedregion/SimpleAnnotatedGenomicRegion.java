package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;
import java.util.stream.Collectors;

/**
 * Simple class that just has an interval and name value pairs.
 */
public class SimpleAnnotatedGenomicRegion implements Locatable {
    private SimpleInterval interval;
    private Map<String, String> annotations;

    public SimpleAnnotatedGenomicRegion(SimpleInterval interval, Map<String, String> annotations) {
        this.interval = interval;
        this.annotations = annotations;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public Map<String, String> getAnnotations() {
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

    @Override
    public boolean equals(Object o) {
        if (this == o) {return true;}
        if (o == null || getClass() != o.getClass()) return false;
        final SimpleAnnotatedGenomicRegion that = (SimpleAnnotatedGenomicRegion) o;
        return this.interval.equals(that.getInterval()) && this.getAnnotations().equals(that.getAnnotations());

    }

    @Override
    public String toString() {
        return interval.toString() + " :: " + StringUtils.join(Utils.stream(annotations.entrySet().iterator())
                .map(e -> e.getKey() + "->" + e.getValue()).collect(Collectors.toList()), ",");
    }
}
