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

    public boolean tileIsSet() {
        return tile != -1;
    }

    public ReadData getRead(final int index) {
        return reads[index];
    }

    public int getNumReads() {
        return reads.length;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckTile(final int tile) {
        if (tileIsSet()) {
            if (this.tile != tile) {
                throw new IlluminaParserException("Tile number mismatch for " + this + " : " + this.tile + " != " + tile);
            }
        } else {
            this.tile = tile;
        }
    }

    public int getLane() {
        return lane;
    }

    public void setLane(final int lane) {
        this.lane = lane;
    }

    public boolean laneIsSet() {
        return lane != -1;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckLane(final int lane) {
        if (laneIsSet()) {
            if (this.lane != lane) {
                throw new IlluminaParserException("Lane number mismatch for " + this + " : " + this.lane + " != " + lane);
            }
        } else {
            this.lane = lane;
        }
    }

    public int getX() {
        return x;
    }

    public void setX(final int x) {
        this.x = x;
    }

    public boolean xIsSet() {
        return x != -1;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckX(final int x) {
        if (xIsSet()) {
            if (this.x != x) {
                throw new IlluminaParserException("X value mismatch for " + this + " : " + this.x + " != " + x);
            }
        } else {
            this.x = x;
        }
    }

    public int getY() {
        return y;
    }

    public void setY(final int y) {
        this.y = y;
    }

    public boolean yIsSet() {
        return y != -1;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckY(final int y) {
        if (yIsSet()) {
            if (this.y != y) {
                throw new IlluminaParserException("Y value mismatch for " + this + " : " + this.y + " != " + y);
            }
        } else {
            this.y = y;
        }
    }

    public Boolean isPf() {
        return pf;
    }

    public void setPf(final boolean pf) {
        this.pf = pf;
    }

    /**
     * Either set this value if not already set, or if already set, throw an exception if new value != current value.
     */
    public void setOrCheckPf(final boolean pf) {
        if (this.pf == null) {
            this.pf = pf;
        } else if (this.pf != pf) {
            throw new IlluminaParserException("PF value mismatch for " + this + " : ");
        }
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
