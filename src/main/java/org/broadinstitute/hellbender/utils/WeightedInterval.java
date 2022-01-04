package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;

public class WeightedInterval extends Interval {
    float weight;

    public WeightedInterval(String sequence, int start, int end, float weight) {
        super(sequence, start, end);
        this.weight = weight;
    }

    public WeightedInterval(Locatable locatable, float weight) {
        super(locatable);
        this.weight = weight;
    }

    public WeightedInterval(String sequence, int start, int end, boolean negative, String name, float weight) {
        super(sequence, start, end, negative, name);
        this.weight = weight;
    }

    public float getWeight() {
        return weight;
    }

    @Override
    public String toString() {
        return getContig() + ":" + getStart() + "-" + getEnd() + "\t" + getStrand().encode() + "\t" + ((null == getName()) ? '.' : getName()) + "\t" + getWeight();
    }
}
