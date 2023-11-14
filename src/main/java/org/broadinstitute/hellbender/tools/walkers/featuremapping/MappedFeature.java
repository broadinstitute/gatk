package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.Map;

public class MappedFeature implements Comparable<MappedFeature> {

    GATKRead read;
    FlowFeatureMapperArgumentCollection.MappingFeatureEnum type;
    byte[] readBases;
    byte[] refBases;
    int readBasesOffset; // offset of read bases array
    int start;      // location (on reference)
    int offsetDelta;
    double score;
    int readCount;
    int filteredCount;
    int nonIdentMBasesOnRead;
    int featuresOnRead;
    int refEditDistance;
    int index;
    int smqLeft;
    int smqRight;
    int smqLeftMean;
    int smqRightMean;

    double scoreForBase[];
    boolean adjacentRefDiff;

    Map<String, Object> vcfAttrCache;

    FlowBasedRead flowRead;

    public MappedFeature(GATKRead read, FlowFeatureMapperArgumentCollection.MappingFeatureEnum type, byte[] readBases,
                         byte[] refBases, int readBasesOffset, int start, int offsetDelta) {
        this.read = read;
        this.type = type;
        this.readBases = readBases;
        this.refBases = refBases;
        this.readBasesOffset = readBasesOffset;
        this.start = start;
        this.offsetDelta = offsetDelta;
    }

    static MappedFeature makeSNV(GATKRead read, int offset, byte refBase, int start, int offsetDelta) {
        byte[] readBases = {read.getBasesNoCopy()[offset]};
        byte[] refBases = {refBase};
        return new MappedFeature(
                read,
                FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV,
                readBases,
                refBases,
                offset,
                start,
                offsetDelta);
    }

    @Override
    public String toString() {
        return "Feature{" +
                "read=" + read +
                ", type=" + type +
                ", readBases=" + Arrays.toString(readBases) +
                ", refBases=" + Arrays.toString(refBases) +
                ", readBasesOffset=" + readBasesOffset +
                ", start=" + start +
                '}';
    }

    @Override
    public int compareTo(MappedFeature o) {

        int delta = this.read.getContig().compareTo(o.read.getContig());
        if (delta != 0) {
            return delta;
        } else {
            return this.start - o.start;
        }
    }
}
