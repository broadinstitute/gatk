package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;


import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public final class EventRecord implements Locatable {
    private final String id;
    private final String type;
    private final SimpleInterval interval;
    private final Set<String> samples;

    public EventRecord(final String id, final String type, final SimpleInterval interval, final String[] samples) {
        this.id = id;
        this.type = type;
        this.interval = interval;
        this.samples = new HashSet<>(Arrays.asList(samples));
    }

    public String getId() {
        return id;
    }

    public String getType() {
        return type;
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public Collection<String> getSamples() {
        return samples;
    }

    public boolean isInSample(final String id) {
        return samples.contains(id);
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }
}
