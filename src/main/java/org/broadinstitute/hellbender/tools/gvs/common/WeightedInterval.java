package org.broadinstitute.hellbender.tools.gvs.common;

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

    public float getWeightPerBase() {
        return weight / (float) this.length();
    }

    public WeightedInterval[] split(int basesInFirstInterval) {

        WeightedInterval left = new WeightedInterval(
                this.getContig(),
                this.getStart(),
                this.getStart() + basesInFirstInterval - 1,
                this.isNegativeStrand(),
                this.getName(),
                this.getWeightPerBase() * basesInFirstInterval);

        WeightedInterval right = new WeightedInterval(
                this.getContig(),
                this.getStart() + basesInFirstInterval,
                this.getEnd(),
                this.isNegativeStrand(),
                this.getName(),
                this.getWeight() - left.getWeight()); // give remainder to right

        return new WeightedInterval[]{left, right};
    }

    @Override
    public Interval clone() {
        return super.clone();
    }

    @Override
    public String toString() {
        return getContig() + ":" + getStart() + "-" + getEnd() + "\t" + getStrand().encode() + "\t" + ((null == getName()) ? '.' : getName()) + "\t" + getWeight();
    }
}

