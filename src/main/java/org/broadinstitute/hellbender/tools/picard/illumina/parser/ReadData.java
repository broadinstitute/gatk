package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * Data for a single end of a paired-end read, a barcode read, or for the entire read if not paired end.
 *
 * @author jburke@broadinstitute.org
 */
public final class ReadData {
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
     * @return Noise values as produced by Illumina software, converted to shorts.
     */
    public FourChannelIntensityData getNoise() {
        return noise;
    }

    public void setNoise(final FourChannelIntensityData noise) {
        this.noise = noise;
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

    /**
     * @return Raw intensity values as produced by Illumina software, converted to shorts.
     */
    public FourChannelIntensityData getRawIntensities() {
        return rawIntensities;
    }

    public void setRawIntensities(final FourChannelIntensityData rawIntensities) {
        this.rawIntensities = rawIntensities;
    }

    public ReadType getReadType() {
        return readType;
    }

    public void setReadType(final ReadType readType) {
        this.readType = readType;
    }
}
