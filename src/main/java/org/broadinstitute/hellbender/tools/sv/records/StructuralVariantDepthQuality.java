package org.broadinstitute.hellbender.tools.sv.records;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.sv.EventCopyNumberPosteriors;

public class StructuralVariantDepthQuality implements Locatable {
    private final EventCopyNumberPosteriors.EventRecord eventRecord;
    private final double quality;

    public StructuralVariantDepthQuality(final EventCopyNumberPosteriors.EventRecord eventRecord, final double quality) {
        this.eventRecord = eventRecord;
        this.quality = quality;
    }

    public EventCopyNumberPosteriors.EventRecord getEventRecord() {
        return eventRecord;
    }

    public double getQuality() {
        return quality;
    }

    @Override
    public String getContig() {
        return eventRecord.getInterval().getContig();
    }

    @Override
    public int getStart() {
        return eventRecord.getInterval().getStart();
    }

    @Override
    public int getEnd() {
        return eventRecord.getInterval().getEnd();
    }
}
