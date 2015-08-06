package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * A  genomic interval with a copy number call (eg duplication, deletion, copy-neutral. . .)
 *
 * @author David Benjamin;
 */
//TO-DO (maybe) change hellbender so that this class extends SimpleInterval
public final class CalledInterval implements Locatable {

    private final SimpleInterval interval;
    private String call;

    /**
     * Construct a new CalledInterval given all its properties.
     * @param interval the interval.
     * @param call the call
     *
     * @throws IllegalArgumentException if {@code interval} is not a valid interval
     */
    public CalledInterval(final SimpleInterval interval, final String call) {
        Utils.nonNull(interval, "Interval can't be null");
        this.interval = interval;
        this.call = call;
    }

    /**
     * Returns the interval.
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    @Override
    public String getContig() {return interval.getContig(); }

    @Override
    public int getStart() {return interval.getStart(); }

    @Override
    public int getEnd() {return interval.getEnd(); }


    /**
     * Returns the call.  Returns null for uncalled segments.
     *
     * @return maybe {@code null}
     */
    public String getCall() {
        return call;
    }

    /**
     * Sets the call.
     */
    public void setCall(final String call) {
        this.call = call;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof CalledInterval)) {
            return false;
        }

        CalledInterval calledInterval = (CalledInterval) o;
        return interval.equals(calledInterval.interval) && call.equals(calledInterval.call);
    }

    @Override
    public int hashCode() {
        if (call == null) {
            return interval.hashCode();
        } else {
            return 31 * interval.hashCode() + (call != null ? call.hashCode() : 0);
        }
    }
}
