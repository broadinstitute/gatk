package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.util.Interval;

public class GtfInfo {

    public enum Type {
        GENE,
        TRANSCRIPT
    }

    private final Type type;
    private final String geneName;
    private final Interval interval;

    public GtfInfo(Interval interval, Type type, String geneName) {
        this.interval = interval;
        this.type = type;
        this.geneName = geneName;
    }

    public Type getType() {
        return type;
    }

    public String getGeneName() {
        return geneName;
    }

    public Interval getInterval() {
        return interval;
    }

    public Integer getStart() {
        return interval.getStart();
    }

    public Integer getEnd() {
        return interval.getEnd();
    }

    @Override
    public String toString() {
        return "GtfInfo{ " + "type = " + type + " geneName = " + geneName + "interval = " + interval;
    }

}
