package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/** A class that implements the IlluminaData interfaces provided by this parser
 * One BclData object is returned to IlluminaDataProvider per cluster and each
 * first level array in bases and qualities represents a single read in that
 * cluster */
public class BclData implements BaseData, QualityData {
    public final byte [][] bases;
    public final byte [][] qualities;

    public BclData(final int[] outputLengths) {
        bases     = new byte[outputLengths.length][];
        qualities = new byte[outputLengths.length][];

        for(int i = 0; i < outputLengths.length; i++) {
            bases[i]     = new byte[outputLengths[i]];
            qualities[i] = new byte[outputLengths[i]];
        }
    }

    @Override
    public byte[][] getBases() {
        return bases;
    }

    @Override
    public byte[][] getQualities() {
        return qualities;
    }
}
