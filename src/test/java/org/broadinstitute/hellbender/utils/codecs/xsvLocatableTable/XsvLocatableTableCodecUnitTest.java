package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

/**
 * A Class to hold unit tests for {@link XsvLocatableTableCodec}.
 * Created by jonn on 12/18/17.
 */
public class XsvLocatableTableCodecUnitTest extends GATKBaseTest {

    //==================================================================================================================
    // Private Static Members:

    private static final String TEST_RESOURCE_DIR = publicTestDir + "org/broadinstitute/hellbender/utils/codecs/xsvLocatableTable" + File.separator;
    private static final String TEST_FILE1 = TEST_RESOURCE_DIR + "xsv_locatable_test.csv";
    private static final String TEST_FILE2 = TEST_RESOURCE_DIR + "xsv_locatable_test2.tsv";
    private static final String TEST_FILE_NO_CONFIG = TEST_RESOURCE_DIR + "xsv_locatable_test_no_config.csv";

    private static final List<String> file1Headers = Arrays.asList("XSV_LOCATABLE_TEST_NAME_Villain", "XSV_LOCATABLE_TEST_NAME_chr", "XSV_LOCATABLE_TEST_NAME_test_val", "XSV_LOCATABLE_TEST_NAME_start", "XSV_LOCATABLE_TEST_NAME_end", "XSV_LOCATABLE_TEST_NAME_Bond");
    private static final List<String> file1Line1 = Arrays.asList("Blofeld", "chr19", "test_val_chr19", "8959519", "9092018", "Connery");
    private static final List<String> file1Line2 = Arrays.asList("Largo", "chr3", "test_val_chr3", "178866310", "178957882", "Dalton");

    private static final List<String> file2Headers = Arrays.asList("SECOND_XSV_NAME_Car_Maker", "SECOND_XSV_NAME_chr", "SECOND_XSV_NAME_start", "SECOND_XSV_NAME_Tire_Maker", "SECOND_XSV_NAME_end", "SECOND_XSV_NAME_Parent_Company");
    private static final List<String> file2Line1 = Arrays.asList("Lamborghini", "chr3", "178866310", "Michelin", "178957882", "Audi");
    private static final List<String> file2Line2 = Arrays.asList("Ferrari", "chr19", "8959519", "Pirelli", "9092018", "Fiat");

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Helper Methods:

    //==================================================================================================================
    // Data Providers:

    @DataProvider
    private Object[][] provideForTestCanDecode() {
        return new Object[][] {
                { TEST_FILE1, true },
                { TEST_FILE2, true },
                { TEST_FILE_NO_CONFIG, false },
        };
    }

    @DataProvider
    private Object[][] provideForTestDecode() {

        return new Object[][] {
                { TEST_FILE1,
                    Arrays.asList(
                        new XsvTableFeature(1, 3, 4, file1Headers, file1Line1, "XSV_LOCATABLE_TEST_NAME"),
                        new XsvTableFeature(1, 3, 4, file1Headers, file1Line2, "XSV_LOCATABLE_TEST_NAME")
                    )
                },
                { TEST_FILE2,
                    Arrays.asList(
                        new XsvTableFeature(1, 2, 4, file2Headers, file2Line1, "SECOND_XSV_NAME"),
                        new XsvTableFeature(1, 2, 4, file2Headers, file2Line2, "SECOND_XSV_NAME")
                    )
                },
        };
    }

    @DataProvider
    private Object[][] provideForTestReadActualHeader() {
        return new Object[][] {
                { TEST_FILE1, file1Headers },
                { TEST_FILE2, file2Headers },
        };
    }

    @DataProvider
    private Object[][] provideForTestGetConfigFilePath() {
        return new Object[][] {
                { TEST_FILE1, IOUtils.getPath(TEST_RESOURCE_DIR + "xsv_locatable_test.config") },
                { TEST_FILE2, IOUtils.getPath(TEST_RESOURCE_DIR + "xsv_locatable_test2.config") },
        };
    }

    @DataProvider
    private Object[][] provideForTestGetAndValidateConfigFileContents() {

        final Properties configFile1Properties = new Properties();
        configFile1Properties.put("contig_column", "1");
        configFile1Properties.put("start_column", "3");
        configFile1Properties.put("end_column", "4");
        configFile1Properties.put("xsv_delimiter", ",");
        configFile1Properties.put("name", "XSV_LOCATABLE_TEST_NAME");

        final Properties configFile2Properties = new Properties();
        configFile2Properties.put("contig_column", "1");
        configFile2Properties.put("start_column", "2");
        configFile2Properties.put("end_column", "4");
        configFile2Properties.put("xsv_delimiter", "\t");
        configFile2Properties.put("name", "SECOND_XSV_NAME");
        
        return new Object[][] {
                {
                        IOUtils.getPath(TEST_RESOURCE_DIR + "xsv_locatable_test.config"),
                        configFile1Properties
                },
                {
                        IOUtils.getPath(TEST_RESOURCE_DIR + "xsv_locatable_test2.config"),
                        configFile2Properties
                },
        };
    }

    //==================================================================================================================
    // Tests:

    // canDecode
    @Test(dataProvider = "provideForTestCanDecode")
    public void testCanDecode(final String filePath, final boolean expected) {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();
        Assert.assertEquals(xsvLocatableTableCodec.canDecode(filePath), expected);
    }

    // decode
    @Test(dataProvider = "provideForTestDecode")
    public void testDecode(final String filePath, final List<XsvTableFeature> expected) {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();
        if (xsvLocatableTableCodec.canDecode(filePath)) {
            try ( final FileInputStream fileInputStream = new FileInputStream(filePath)) {

                // Lots of scaffolding to do reading here:
                final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));
                final ArrayList<XsvTableFeature> output = new ArrayList<>(expected.size());

                // Read off the header:
                xsvLocatableTableCodec.readActualHeader(lineReaderIterator);

                // Read and decode the lines:
                while ( lineReaderIterator.hasNext() ) {
                    output.add( xsvLocatableTableCodec.decode(lineReaderIterator.next()) );
                }

                // Make sure we have what we expected:
                Assert.assertEquals(output, expected);
            }
            catch ( final FileNotFoundException ex ) {
                throw new GATKException("Error - could not find test file: " + filePath, ex);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Error - IO problem with file " + filePath, ex);
            }
        }
        else {
            throw new GATKException("Error - bad test case.");
        }
    }

    // readActualHeader
    @Test(dataProvider = "provideForTestReadActualHeader")
    public void testReadActualHeader(final String filePath, final List<String> expected) {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();

        if (xsvLocatableTableCodec.canDecode(filePath)) {
            try ( final FileInputStream fileInputStream = new FileInputStream(filePath)) {

                // Lots of scaffolding to do reading here:
                final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));

                // Make sure we have what we expected:
                Assert.assertEquals(xsvLocatableTableCodec.readActualHeader(lineReaderIterator), expected);
            }
            catch ( final FileNotFoundException ex ) {
                throw new GATKException("Error - could not find test file: " + filePath, ex);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Error - IO problem with file " + filePath, ex);
            }
        }
        else {
            throw new GATKException("Error - bad test case.");
        }
    }

    // getConfigFilePath
    @Test(dataProvider = "provideForTestGetConfigFilePath")
    public void testGetConfigFilePath(final String filePath, final Path expected) {
        Assert.assertEquals(XsvLocatableTableCodec.getConfigFilePath(IOUtils.getPath(filePath)), expected);
    }

    // getAndValidateConfigFileContents
    @Test(dataProvider = "provideForTestGetAndValidateConfigFileContents")
    public void testGetAndValidateConfigFileContents(final Path configFilePath, final Properties expected) {
        final Properties properties = XsvLocatableTableCodec.getAndValidateConfigFileContents(configFilePath);
        Assert.assertEquals(properties, expected);
    }

}
