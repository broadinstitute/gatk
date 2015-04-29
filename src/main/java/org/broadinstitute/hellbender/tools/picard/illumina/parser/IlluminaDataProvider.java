package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import java.io.File;
import java.util.*;

/**
 * Parse various formats and versions of Illumina Basecall files, and use them the to populate
 * ClusterData objects.  Clients of this code should use IlluminaDataProviderFactory to create an IlluminaDataProvider.
 * IlluminaDataProvider is immutable after construction.
 *
 * @author jburke@broadinstitute.org
 */
public class IlluminaDataProvider implements Iterator<ClusterData>, Iterable<ClusterData> {

    /**
     * contains QSeqs, bcls, or other Illumina file types that will be parsed by this class
     */
    private final File basecallDirectory; //These two are for error reporting only
    private final int lane;

    /**
     * A list of parsers (already initialized) that should output data in a format consistent with readStructure
     */
    private final IlluminaParser[] parsers;

    /**
     * for each parser in this.parsers there is an array of IlluminaDataTypes that specifies what datatypes that parser is providing in
     * this particular run.  A parser may be able to provide data types which may not be listed here because client code may not
     * have specified these data types
     */
    private final IlluminaDataType[][] dataTypes;

    /**
     * Calculated once, outputReadTypes describes the type of read data for each ReadData that will be found in output ClusterData objects
     */
    private final ReadType[] outputReadTypes;

    /**
     * Number of reads in each ClusterData
     */
    private final int numReads;

    /**
     * Create an IlluminaDataProvider given a map of parsersToDataTypes for particular file formats.  Compute once the miscellaneous data for the
     * run that will be passed to each ClusterData.
     *
     * @param basecallDirectory For error reporting only.
     * @param lane              For error reporting only.
     */
    IlluminaDataProvider(final OutputMapping outputMapping,
                         final Map<IlluminaParser, Set<IlluminaDataType>> parsersToDataTypes,
                         final File basecallDirectory, final int lane) {
        this.basecallDirectory = basecallDirectory;
        this.lane = lane;
        numReads = outputMapping.numOutputReads();

        final int numParsers = parsersToDataTypes.size();
        if (numParsers == 0) {
            throw new IlluminaParserException("There were 0 parsers passed to IlluminaDataProvider!");
        }

        int i = 0;
        parsers = new IlluminaParser[numParsers];
        dataTypes = new IlluminaDataType[numParsers][];
        for (final Map.Entry<IlluminaParser, Set<IlluminaDataType>> pToD : parsersToDataTypes.entrySet()) {
            parsers[i] = pToD.getKey();
            final Set<IlluminaDataType> dts = pToD.getValue();
            dataTypes[i] = new IlluminaDataType[dts.size()];
            dts.toArray(dataTypes[i++]);
        }

        this.outputReadTypes = new ReadType[numReads];
        i = 0;
        for (final ReadDescriptor rd : outputMapping.getOutputDescriptors()) {
            outputReadTypes[i++] = rd.type;
        }
    }

    /**
     * @return True if we have more clusters to read
     */
    public boolean hasNext() {
        final boolean more = parsers[0].hasNext();
        if (!more) {
            for (int i = 1; i < parsers.length; i++) {
                if (parsers[i].hasNext()) {
                    throw new IlluminaParserException("Unequal length Illumina files in " + basecallDirectory + ", lane " + lane
                            + ". Failing parser: " + parsers[i].getClass().getName());
                }
            }
        }

        return more;
    }

    /**
     * @return Current cluster data populated with only the data that matches one of the data types in dataTypes.
     */
    public ClusterData next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        final ClusterData cluster = new ClusterData(outputReadTypes);
        cluster.setLane(lane);

        //IMPORTANT NOTE: This assignment to tile MUST happen BEFORE the loop below because getTileOfNextCluster
        //returns the tile for the next cluster and if we call this after the loop then whenever we pass a tile
        //boundary the last cluster in the previous tile will have the wrong tile number
        cluster.setTile(parsers[0].getTileOfNextCluster());

        for (int i = 0; i < parsers.length; i++) {
            final IlluminaData ilData = parsers[i].next();
            for (final IlluminaDataType ilDataType : dataTypes[i]) {
                switch (ilDataType) {
                    case Position:
                        addData(cluster, (PositionalData) ilData);
                        break;

                    case PF:
                        addData(cluster, (PfData) ilData);
                        break;

                    case Barcodes:
                        addData(cluster, (BarcodeData) ilData);
                        break;

                    case BaseCalls:
                        addReadData(cluster, numReads, (BaseData) ilData);
                        break;

                    case QualityScores:
                        addReadData(cluster, numReads, (QualityData) ilData);
                        break;

                    default:
                        throw new IlluminaParserException("Unknown data type " + ilDataType + " requested by IlluminaDataProviderFactory");
                }
            }
        }

        return cluster;
    }

    /*
     * Methods for that transfer data from the IlluminaData objects to the current cluster
     */
    private void addData(final ClusterData clusterData, final PositionalData posData) {
        clusterData.setX(posData.getXCoordinate());
        clusterData.setY(posData.getYCoordinate());
    }

    private void addData(final ClusterData clusterData, final PfData pfData) {
        clusterData.setPf(pfData.isPf());
    }

    private void addData(final ClusterData clusterData, final BarcodeData barcodeData) {
        clusterData.setMatchedBarcode(barcodeData.getBarcode());
    }

    private void addReadData(final ClusterData clusterData, final int numReads, final BaseData baseData) {
        final byte[][] bases = baseData.getBases();
        for (int i = 0; i < numReads; i++) {
            clusterData.getRead(i).setBases(bases[i]);
        }
    }

    private void addReadData(final ClusterData clusterData, final int numReads, final QualityData qualityData) {
        final byte[][] qualities = qualityData.getQualities();
        for (int i = 0; i < numReads; i++) {
            clusterData.getRead(i).setQualities(qualities[i]);
        }
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Jump so that the next record returned will be from the specified tile.
     */
    public void seekToTile(final int oneBasedTileNumber) {
        for (final IlluminaParser parser : parsers) {
            parser.seekToTile(oneBasedTileNumber);
        }
    }

    @Override
    public Iterator<ClusterData> iterator() {
        return this;
    }

    public void close() {
        for (final IlluminaParser parser : parsers) {
            parser.close();
        }
    }
}
