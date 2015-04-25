package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * Store the information from Illumina files for a single cluster with one or more reads.
 *
 * @author jburke@broadinstitute.org
 */
public class ClusterData {

    private int lane = -1;
    private int tile = -1;
    private int x = -1;
    private int y = -1;
    private final ReadData[] reads;
    private Boolean pf;
    private String matchedBarcode;

    /**
     * Used for testing, reads is set directly with no copying to the input array
     */
    public ClusterData(final ReadData... reads) {
        this.reads = reads;
    }

    /**
     * Creates a ClusterData with one read for each type provided
     */
    public ClusterData(final ReadType[] readTypes) {
        reads = new ReadData[readTypes.length];
        for (int i = 0; i < readTypes.length; i++) {
            reads[i] = new ReadData(readTypes[i]);
        }
    }

    public String toString() {
        return "ClusterData(lane: " + lane + "; tile: " + tile + "; x: " + x + "; y: " + y + "; pf: " + pf +
                "; matchedBarcode: " + matchedBarcode + ")";
    }

    public int getTile() {
        return tile;
    }

    public void setTile(final int tile) {
        this.tile = tile;
    }

    public ReadData getRead(final int index) {
        return reads[index];
    }

    public int getNumReads() {
        return reads.length;
    }

    public int getLane() {
        return lane;
    }

    public void setLane(final int lane) {
        this.lane = lane;
    }

    public int getX() {
        return x;
    }

    public void setX(final int x) {
        this.x = x;
    }

    public int getY() {
        return y;
    }

    public void setY(final int y) {
        this.y = y;
    }

    public Boolean isPf() {
        return pf;
    }

    public void setPf(final boolean pf) {
        this.pf = pf;
    }

    /**
     * @return The barcode matched (not the actual sequence from the read, which may not perfectly match
     * the barcode).
     */
    public String getMatchedBarcode() {
        return matchedBarcode;
    }

    public void setMatchedBarcode(final String matchedBarcode) {
        this.matchedBarcode = matchedBarcode;
    }
}
