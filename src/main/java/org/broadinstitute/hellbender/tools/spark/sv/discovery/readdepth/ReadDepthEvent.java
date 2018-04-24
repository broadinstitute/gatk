package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.LargeSimpleSV;

import java.io.Serializable;

final class ReadDepthEvent implements Serializable {

    public static final long serialVersionUID = 1L;
    private final int id;
    private final LargeSimpleSV event;
    private double state;
    public int leftDistance;
    public int rightDistance;

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
        return LargeSimpleSV.getBedHeader() + "\tID\tSTATE\tDIST1\tDIST2";
    }

    public String toBedString(final SAMSequenceDictionary dictionary, final double counterEvidencePseudocount) {
        return event.toBedString(dictionary, counterEvidencePseudocount) + "\t" + id + "\t" + state + "\t" + leftDistance + "\t" + rightDistance;
    }
}
