package org.broadinstitute.hellbender.tools.walkers.conversion;

import htsjdk.samtools.util.Interval;

public class GtfInfo {
    public enum Type{
        GENE,
        TRANSCRIPT
    }

    private Type type;
    private String geneName;
    private Interval interval;

    public GtfInfo(Interval interval, Type type, String geneName){
        this.interval = interval;
        this.type = type;
        this.geneName = geneName;
    }

    public Type getType(){
        return type;
    }

    public void setType(Type type){
        this.type = type;
    }

    public String getGeneName(){
        return geneName;
    }

    public void setGeneName(String geneName){
        this.geneName = geneName;
    }

    public Interval getInterval(){
        return interval;
    }

    public Integer getStart(){
        return interval.getStart();
    }

    public Integer getEnd(){
        return interval.getEnd();
    }

    public void setInterval(Interval interval){
        this.interval = interval;
    }

    @Override
    public String toString(){
        return "GtfInfo{ " + "type = " + type + " geneName = " + geneName + "interval = " + interval;
    }

}
