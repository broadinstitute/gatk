package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;

/**
 * Exome analysis segment.
 *
 * A  genomic interval optionally attached to it</p>
 *
 * @author David Benjamin;
 */
//TO-DO (maybe) change hellbender so that this class extends SimpleInterval
public final class Segment implements Locatable {

    private final String sample;

    private final SimpleInterval interval;

    private String call = null;

    /**
     * Construct a new Segment with an empty call
     * @param sample the sample name
     * @param interval the interval.
     *
     * @throws IllegalArgumentException if either {@code interval}
     * does not represent a valid interval
     */
    public Segment(final String sample, final SimpleInterval interval) {
        Utils.nonNull(sample, "Segment needs a sample name.");
        Utils.nonNull(interval, "Can't construct segment with null interval.");
        this.sample = sample;
        this.interval = interval;

    }

    /**
     * Construct a new Segment given all its properties.
     * @param sample the sample name
     * @param interval the interval.
     * @param call the call
     *
     * @throws IllegalArgumentException if either {@code interval} or {@code start}
     *    and {@code end} do not represent a valid interval
     */
    public Segment(final String sample, final SimpleInterval interval, final String call) {
        this(sample, interval);
        this.call = call;
    }

    /**
     * Returns the sample.
     *
     * @return never {@code null}
     */
    public String getSample() {
        return sample;
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
     * Returns the overlapping targets.
     *
     * Delegates to ExonCollection's binary search.
     */
    public List<TargetCoverage> overlappingTargets(final ExonCollection<TargetCoverage> targets) {
        return targets.exons(interval);
    }

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

    /**
     * the mean of all overlapping targets' coverages
     *
     * @throws IllegalStateException if overlapping targets have not been assigned or if no overlapping targets were found.
     */
    public double mean(final ExonCollection<TargetCoverage> targets) {
        final List<TargetCoverage> myTargets = overlappingTargets(targets);

        if (myTargets.size() == 0) {
            throw new IllegalStateException("Empty segment -- no overlapping targets.");
        }
        return myTargets.stream().mapToDouble(TargetCoverage::getCoverage).average().getAsDouble();
    }

    public int numTargets(final ExonCollection<TargetCoverage> targets) {
        return overlappingTargets(targets).size();
    }

    public boolean overlaps(final Segment other) {
        Utils.nonNull(other, "Testing for overlap with null segment.");
        return this.interval.overlaps(other.interval);
    }

}
