package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import scala.Tuple2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

final class ReadDepthEvent implements Serializable {

    public static final long serialVersionUID = 1L;
    private final int id;
    private final LargeSimpleSV event;
    private double state;
    public double leftDistance;
    public double rightDistance;
    public double observedCopyNumber;
    public double linkDensity;
    public double mappabilityIndex;
    public final List<Integer> overlappingEventIds;
    public double copyNumberCallOverlap;

    public double copyNumberLikelihood;
    public double distanceLikelihood;
    public double copyNumberCallOverlapLikelihood;
    public double readPairEvidenceLikelihood;
    public double splitReadEvidenceLikelihood;
    public double snpRate;
    public double probability;
    public boolean isTrue;
    public List<ReadDepthModel.OverlapInfo> overlapInfoList;
    public List<Tuple2<Integer,Double>> optimizedOverlapInfoList;

    public ReadDepthEvent(final int id, final LargeSimpleSV event) {
        this.event = event;
        this.id = id;
        this.state = 0;
        this.overlappingEventIds = new ArrayList<>();
        this.isTrue = false;
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
        return LargeSimpleSV.getBedHeader() + "\tID\tSTATE\tMEAN_CN\tDIST1\tDIST2\tLINK_DENSITY\tMAPPABILITY\tCALL_OVERLAP\tSNP_RATE\tCN_LIK\tDIST_LIK\tCALL_OVERLAP_LIK\tRP_LIK\tSR_LIK\tP\tTRUTH";
    }

    public String toBedString(final SAMSequenceDictionary dictionary) {
        return event.toBedString(dictionary) + "\t" + id + "\t" + state + "\t" + observedCopyNumber
                + "\t" + leftDistance + "\t" + rightDistance + "\t" + linkDensity + "\t" + mappabilityIndex + "\t" + copyNumberCallOverlap + "\t" + snpRate
                + "\t" + copyNumberLikelihood + "\t" + distanceLikelihood + "\t" + copyNumberCallOverlapLikelihood + "\t" + readPairEvidenceLikelihood + "\t" + splitReadEvidenceLikelihood + "\t" + probability + "\t" + (isTrue ? 1 : 0);
    }
}
