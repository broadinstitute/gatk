package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;


/**
 * Extension of Interval that also stores a weight for the interval
 */
public class WeightedInterval extends Interval {
    private final float weight;

    public WeightedInterval(String sequence, int start, int end, float weight) {
        this(sequence, start, end, false, null, weight);
    }

    public WeightedInterval(Locatable locatable, float weight) {
        this(locatable.getContig(), locatable.getStart(), locatable.getEnd(), weight);
    }

    public WeightedInterval(String sequence, int start, int end, boolean negative, String name, float weight) {
        super(sequence, start, end, negative, name);

        if (end - start < 0) {
            throw new IllegalArgumentException("length must be > 0");
        }

        if (weight < 0) {
            throw new IllegalArgumentException("Weight must be > 0");
        }

        this.weight = weight;
    }

    public float getWeight() {
        return weight;
    }

    /**
     * @return the total weight divided by the length of the interval
     */
    public float getWeightPerBase() {
        return weight / (float) this.length();
    }

    /**
     * Splits this interval into two WeightedIntervals where the first
     * interval is from [start, start + basesInFirstInterval - 1] and the
     * second interval is from [start + basesInFirstInterval, end].  The weight is
     * proportionally distributed across the two intervals by their lengths
     *
     * @return an array of two WeightedIntervals
     */
    public WeightedInterval[] split(int basesInFirstInterval) {
        if (basesInFirstInterval < 1 || basesInFirstInterval >= length()) {
            throw new IllegalArgumentException("basesInFirstInterval must be >= 1 and < the length of this interval");
        }

        WeightedInterval left = new WeightedInterval(
                this.getContig(),
                this.getStart(),
                this.getStart() + basesInFirstInterval - 1,
                this.isNegativeStrand(),
                this.getName() == null ? null : this.getName() + "-1", // ensure non-null names are unique
                this.getWeightPerBase() * basesInFirstInterval);

        WeightedInterval right = new WeightedInterval(
                this.getContig(),
                this.getStart() + basesInFirstInterval,
                this.getEnd(),
                this.isNegativeStrand(),
                this.getName() == null ? null : this.getName() + "-2", // ensure non-null names are unique
                this.getWeight() - left.getWeight()); // give remainder to right

        return new WeightedInterval[]{left, right};
    }

    @Override
    public WeightedInterval clone() {
        return (WeightedInterval) super.clone();
    }

    @Override
    public String toString() {
        return getContig() + ":" + getStart() + "-" + getEnd() + "\t" + getStrand().encode() + "\t" + ((null == getName()) ? '.' : getName()) + "\t" + getWeight();
    }
}

