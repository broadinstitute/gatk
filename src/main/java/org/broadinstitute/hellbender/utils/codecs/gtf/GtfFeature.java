package org.broadinstitute.hellbender.utils.codecs.gtf;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;

import java.util.Map;

public class GtfFeature implements Feature {

    private final String contig;
    private final String source;
    private final String type;
    private final int start;
    private final int end;
    private final Strand strand;
    private final int phase;
    private final Map<String, String> attributes;

    protected GtfFeature(final String contig, final String source, final String type, final int start, final int end, final Strand strand, final int phase, final Map<String, String> attributes) {
        this.contig = contig;
        this.source = source;
        this.type = type;
        this.start = start;
        this.end = end;
        this.phase = phase;
        this.strand = strand;
        this.attributes = attributes;
    }

    public String getSource() {
        return source;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public Strand getStrand() {
        return strand;
    }

    public int getPhase() {
        return phase;
    }

    public String getType() {return type;}

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    public String getAttribute(String key) {
        return attributes.get(key);
    }

    @Override
    public boolean equals(Object other) {
        if (!other.getClass().equals(GtfFeature.class)) {
            return false;
        }

        final GtfFeature otherGtfFeature = (GtfFeature) other;

        return (otherGtfFeature.getContig().equals(contig) &&
                otherGtfFeature.getSource().equals(source) &&
                otherGtfFeature.getType().equals(type) &&
                otherGtfFeature.getStart() == start &&
                otherGtfFeature.getEnd()== end &&
                otherGtfFeature.getStrand() == strand &&
                otherGtfFeature.getPhase() == phase &&
                attributes.equals(attributes)
                );
    }

    @Override
    public int hashCode() {
        return contig.hashCode() + source.hashCode() + type.hashCode() + new Integer(start).hashCode() + new Integer(end).hashCode() + strand.hashCode() + new Integer(phase).hashCode() + attributes.hashCode();
    }
}
