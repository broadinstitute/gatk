package org.broadinstitute.hellbender.tools.picard.illumina.parser;

/**
 * There is one IlluminaData sub-interface for each IlluminaDataType enum value.
 * IlluminaParsers must return objects implementing at least one of the interfaces below.
 * IlluminaDataProvider will take IlluminaData objects created by IlluminaParsers and cast them to the types they
 * implement and these objects will then be used to populate the ClusterData object.
 *
 * @author jburke@broadinstitute.org
 */
interface IlluminaData {
}

// Note: PositionalData was spun out this round but since every parser has means of retrieving lane/tile from the
// file name, we are going to move lane/tile to be queryable from parsers in future revisions and therefore if you
// want lane/tile info you will NOT have to parse one of the Positional Data formats (pos, locs, clocs, qseqs)
interface PositionalData extends IlluminaData {
    public int getXCoordinate();

    public int getYCoordinate();
}

interface BaseData extends IlluminaData {
    public byte[][] getBases();
}

interface QualityData extends IlluminaData {
    public byte[][] getQualities();
}

interface NoiseData extends IlluminaData {
    public FourChannelIntensityData[] getNoise();
}

interface RawIntensityData extends IlluminaData {
    public FourChannelIntensityData[] getRawIntensities();
}

interface PfData extends IlluminaData {
    public boolean isPf();
}

interface BarcodeData extends IlluminaData {
    public String getBarcode();
}


