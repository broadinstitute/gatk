package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

public class WeightedInterval extends Interval {
    double weight;

    public WeightedInterval(String sequence, int start, int end, double weight) {
        super(sequence, start, end);
        ParamUtils.isPositive(weight, "weight must be > 0 but it was " + weight);
        this.weight = weight;
    }

    public WeightedInterval(Locatable locatable, double weight) {
        super(locatable);
        ParamUtils.isPositive(weight, "weight must be > 0 but it was " + weight);
        this.weight = weight;
    }

    public WeightedInterval(String sequence, int start, int end, boolean negative, String name, double weight) {
        super(sequence, start, end, negative, name);
        ParamUtils.isPositive(weight, "weight must be > 0 but it was " + weight);
        this.weight = weight;
    }

    public double getWeight() {
        return weight;
    }

    public double getPerBaseWeight() {
        return weight / length();
    }

    @Override
    public String toString() {
        return getContig() + ":" + getStart() + "-" + getEnd() + "\t" + getStrand().encode() + "\t" + ((null == getName()) ? '.' : getName()) + "\t" + getWeight();
    }
}
