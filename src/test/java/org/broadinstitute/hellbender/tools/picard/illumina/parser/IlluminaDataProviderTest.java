package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import org.broadinstitute.hellbender.tools.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static htsjdk.samtools.util.CollectionUtil.makeList;

/**
 * @author jburke@broadinstitute.org
 */

public final class IlluminaDataProviderTest {

    public static final BclQualityEvaluationStrategy bclQualityEvaluationStrategy = new BclQualityEvaluationStrategy(BclQualityEvaluationStrategy.ILLUMINA_ALLEGED_MINIMUM_QUALITY);
    public static final File BINARY_TD_LOCATION = new File("src/test/resources/org/broadinstitute/hellbender/tools/picard/illumina/25T8B25T/Data/Intensities/BaseCalls");
    private static final IlluminaDataType[] DEFAULT_DATA_TYPES = new IlluminaDataType[]{
            IlluminaDataType.Position, IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores, IlluminaDataType.PF
    };

    private void runTest(
            final String testName, final int size,
            final Map<Integer, ClusterData> readNoToClusterData,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final IlluminaDataProvider dataProvider)
            throws Exception {

        int count = 0;
        int readNum = 0;
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();
            if (readNoToClusterData.containsKey(readNum)) {
                compareReadData(cluster, readNoToClusterData.get(readNum), testName + " cluster num " + readNum);
            }

            if (seekAfterFirstRead != 0 && count == 0) {
                dataProvider.seekToTile(seekAfterFirstRead);
                readNum += seekTestDataReadOffset;
            }

            readNum++;
            count++;
        }
        Assert.assertEquals(count, size, testName);
        dataProvider.close();
    }

    private IlluminaDataType[] getDataTypes(final IlluminaDataType[] extraDataTypes) {
        final IlluminaDataType[] dts;

        if (extraDataTypes == null) {
            dts = DEFAULT_DATA_TYPES;
        } else {
            dts = Arrays.copyOf(DEFAULT_DATA_TYPES, DEFAULT_DATA_TYPES.length + extraDataTypes.length);
            System.arraycopy(extraDataTypes, 0, dts, DEFAULT_DATA_TYPES.length, extraDataTypes.length);
        }
        return dts;
    }


    private void compareBasesAndQuals(final ReadData rd1, final ReadData rd2, final String testName) {
        Assert.assertEquals(rd1.getBases(), rd2.getBases(), testName);
        Assert.assertEquals(rd1.getQualities(), rd2.getQualities(), testName);
        Assert.assertEquals(rd1.getReadType(), rd2.getReadType());
    }

    private void comparePositionalData(final ClusterData cd1, final ClusterData cd2, final String testName) {
        Assert.assertEquals(cd1.getLane(), cd2.getLane(), testName);
        Assert.assertEquals(cd1.getTile(), cd2.getTile(), testName);
        Assert.assertEquals(cd1.getX(), cd2.getX(), testName);
        Assert.assertEquals(cd1.getY(), cd2.getY(), testName);
    }

    //Doesn't compare intensities right now -- Do we want too?
    private void compareReadData(final ClusterData cd1, final ClusterData cd2, final String testName) {
        comparePositionalData(cd1, cd2, testName);
        Assert.assertEquals(cd1.getNumReads(), cd2.getNumReads());
        for (int i = 0; i < cd1.getNumReads(); i++) {
            compareBasesAndQuals(cd1.getRead(i), cd2.getRead(i), testName);
        }

        Assert.assertEquals(cd1.getMatchedBarcode(), cd2.getMatchedBarcode(), testName);
        Assert.assertEquals(cd1.isPf().booleanValue(), cd2.isPf().booleanValue(), testName);
    }

    public void runBarcodeParsingTest(final IlluminaDataProviderFactory factory) {
        int total = 0;
        final IlluminaDataProvider dataProvider = factory.makeDataProvider();
        while (dataProvider.hasNext()) {
            final ClusterData cluster = dataProvider.next();
            final String matchedBarcode = cluster.getMatchedBarcode();
            if (matchedBarcode != null) {
                Assert.assertEquals(matchedBarcode, new String(cluster.getRead(1).getBases()));
            }
            if(total > 10){
                break;
            }
            total++;
        }
        dataProvider.close();
    }

    @Test(enabled = false, description = "bug https://github.com/broadinstitute/hellbender/issues/364")
    public void barcodeParsingTest() {
        runBarcodeParsingTest(new IlluminaDataProviderFactory(BINARY_TD_LOCATION, 1, new ReadStructure("25T8B25T"), bclQualityEvaluationStrategy, IlluminaDataType.BaseCalls,
                IlluminaDataType.Barcodes));
    }

    @DataProvider(name = "binaryData")
    public Object[][] binaryData() {
        return new Object[][]{
                {
                        "Bustard Parsing Test(25T8B25T) w/Clocs", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8B25T",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25T8S25T) w/Clocs", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8S25T",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25T8S25T) w/Clocs with ending skip", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8B1S",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25S8S25T) w/Clocs", 1, 180,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25S8S25T",
                        0, 0,
                        BINARY_TD_LOCATION
                },
                {
                        "Bustard Parsing Test(25T8B25T) w/Clocs And Seeking", 1, 61,
                        makeList(1101, 1201, 2101),
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8B25T",
                        2101, 4631,
                        BINARY_TD_LOCATION
                }
        };
    }

    @Test(dataProvider = "binaryData", enabled = false, description = "bug https://github.com/broadinstitute/hellbender/issues/364")
    public void testIlluminaDataProviderBclMethod(
            final String testName, final int lane, final int size,
            final List<Integer> tiles,
            final IlluminaDataType[] extraDataTypes,
            final String illuminaConfigStr,
            final int seekAfterFirstRead, final int seekTestDataReadOffset,
            final File basecallsDirectory)
            throws Exception {

        final IlluminaDataType[] dts = getDataTypes(extraDataTypes);

        final Map<Integer, ClusterData> readNoToClusterData = BinTdUtil.clusterData(lane, tiles, illuminaConfigStr, dts);
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), bclQualityEvaluationStrategy, dts);
        final IlluminaDataProvider dataProvider = factory.makeDataProvider();

        runTest(testName, size, readNoToClusterData, seekAfterFirstRead, seekTestDataReadOffset, dataProvider);
    }

    //Unlike above, the data types here do not have DEFAULT_DATA_TYPES added before creating the dataProvider
    @DataProvider(name = "badData")
    public Object[][] badData() {
        return new Object[][]{
                {
                        5,
                        new IlluminaDataType[]{IlluminaDataType.Barcodes},
                        "25T8B25T",
                        BINARY_TD_LOCATION
                },
                {
                        4,
                        DEFAULT_DATA_TYPES,
                        "25TB25T",
                        BINARY_TD_LOCATION
                },
                {
                        4,
                        DEFAULT_DATA_TYPES,
                        "25T0B25T",
                        BINARY_TD_LOCATION
                },
                {
                        4,
                        DEFAULT_DATA_TYPES,
                        "-25T0B25T",
                        BINARY_TD_LOCATION
                },
                {
                        9,
                        new IlluminaDataType[]{IlluminaDataType.Position, IlluminaDataType.Barcodes},
                        "25T8B25T",
                        BINARY_TD_LOCATION
                },
                {
                        9,
                        new IlluminaDataType[]{IlluminaDataType.BaseCalls},
                        "25T8B25T",
                        BINARY_TD_LOCATION
                },
                {
                        9,
                        new IlluminaDataType[]{IlluminaDataType.PF, IlluminaDataType.BaseCalls, IlluminaDataType.QualityScores},
                        "25T8B24T",
                        BINARY_TD_LOCATION
                }
        };
    }

    @Test(dataProvider = "badData", expectedExceptions = {IlluminaParserException.class, IllegalArgumentException.class}, enabled = false, description = "bug https://github.com/broadinstitute/hellbender/issues/364")
    public void testIlluminaDataProviderMissingDatas(final int lane,
                                                     final IlluminaDataType[] actualDts,
                                                     final String illuminaConfigStr,
                                                     final File basecallsDirectory)
            throws Exception {
        final IlluminaDataProviderFactory factory = new IlluminaDataProviderFactory(basecallsDirectory, lane, new ReadStructure(illuminaConfigStr), bclQualityEvaluationStrategy, actualDts);
        factory.makeDataProvider();
    }
}
