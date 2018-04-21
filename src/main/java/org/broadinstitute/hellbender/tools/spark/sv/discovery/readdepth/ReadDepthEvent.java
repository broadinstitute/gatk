package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;

final class ReadDepthEvent {
    private final int id;
    private final LargeSimpleSV event;
    private double state;

    public ReadDepthEvent(final int id, final LargeSimpleSV event) {
        this.event = event;
        this.id = id;
        this.state = 0;
    }

    public double getState() {
        return state;
    }

    public void setState(final double x) {
        state = x;
    }

    public int getId() {
        return id;
    }

    public LargeSimpleSV getEvent() {
        return event;
    }

    public static String getBedHeader() {
        return LargeSimpleSV.getBedHeader() + "\tID\tSTATE";
    }

    public String toBedString(final SAMSequenceDictionary dictionary, final double counterEvidencePseudocount) {
        return event.toBedString(dictionary, counterEvidencePseudocount) + "\t" + id + "\t" + state;
    }
}
