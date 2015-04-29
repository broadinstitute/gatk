package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * Data for a single end of a paired-end read, a barcode read, or for the entire read if not paired end.
 *
 * @author jburke@broadinstitute.org
 */
public class ReadData {
    private byte[] bases;
    private byte[] qualities;
    private FourChannelIntensityData rawIntensities;
    private FourChannelIntensityData noise;
    private ReadType readType;

    public ReadData() {
    }

    public ReadData(ReadType readType) {
        this.readType = readType;
    }

    /**
     * @return ASCII byte representation of bases.
     */
    public byte[] getBases() {
        return bases;
    }

    public void setBases(final byte[] bases) {
        this.bases = bases;
    }

    /**
     * @return Phred-binary scaled qualities.  E.g. Q20 is the byte with value==20.
     */
    public byte[] getQualities() {
        return qualities;
    }

    public void setQualities(final byte[] qualities) {
        this.qualities = qualities;
    }

    public ReadType getReadType() {
        return readType;
    }

    public void setReadType(final ReadType readType) {
        this.readType = readType;
    }
}
