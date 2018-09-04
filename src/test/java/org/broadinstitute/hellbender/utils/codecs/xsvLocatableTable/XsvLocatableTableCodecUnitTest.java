package org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
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

    private static final String TEST_FILE_MIXED_ENCODING = TEST_RESOURCE_DIR + "xsv_locatable_test_mixed_encodings.csv";

    /** Uses column names, instead of index */
    private static final String TEST_FILE3 = TEST_RESOURCE_DIR + "xsv_locatable_test3.csv";
    private static final String TEST_FILE4 = TEST_RESOURCE_DIR + "xsv_locatable_test4.csv";
    private static final String TEST_FILE_NO_CONFIG = TEST_RESOURCE_DIR + "xsv_locatable_test_no_config.csv";

    // Preambles of SAMFileHeaders or just plain ol' comments
    private static final String TEST_FILE_SAMFILEHEADER = TEST_RESOURCE_DIR + "xsv_locatable_test_samfileheader.tsv";
    private static final String TEST_FILE_SAMFILEHEADER_CONFIG = TEST_RESOURCE_DIR + "xsv_locatable_test_samfileheader.config";
    private static final String TEST_FILE_SAMFILEHEADER_CONFIG_MULTIPLE_COLUMNS = TEST_RESOURCE_DIR + "xsv_locatable_test_samfileheader_multiple_columns.config";
    private static final String TEST_FILE_SAMFILEHEADER_CONFIG_MULTIPLE_COLUMNS_NO_NAME = TEST_RESOURCE_DIR + "xsv_locatable_test_samfileheader_multiple_columns_no_name.config";
    private static final String TEST_FILE_SAMFILEHEADER_CONFIG_ERROR_NOTHING_FOUND = TEST_RESOURCE_DIR + "xsv_locatable_test_samfileheader_error_nothing_found.config";

    private static final String TEST_FILE_MIXED_PREAMBLE = TEST_RESOURCE_DIR + "xsv_locatable_test_mixed_preamble.tsv";
    private static final String TEST_FILE_MIXED_PREAMBLE_CONFIG = TEST_RESOURCE_DIR + "xsv_locatable_test_mixed_preamble.config";

    private static final String TEST_FILE_SAME_STARTEND = TEST_RESOURCE_DIR + "xsv_locatable_test_same_startend.csv";
    private static final String TEST_FILE_SAME_STARTEND_CONFIG = TEST_RESOURCE_DIR + "xsv_locatable_test_same_startend.config";
    private static final String TEST_FILE_SAME_STARTEND_CONFIG_NO_NAME = TEST_RESOURCE_DIR + "xsv_locatable_test_same_startend_no_name.config";

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
                { TEST_FILE3, true },
        };
    }

    @DataProvider
    private Object[][] provideForTestDecodeCharsetFailure() {

        return new Object[][]{
                { TEST_FILE_MIXED_ENCODING },
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
                { TEST_FILE3,
                    Arrays.asList(
                            new XsvTableFeature(1, 3, 4, file1Headers, file1Line1, "XSV_LOCATABLE_TEST_NAME"),
                            new XsvTableFeature(1, 3, 4, file1Headers, file1Line2, "XSV_LOCATABLE_TEST_NAME")
                    )
                },
                { TEST_FILE4,
                    Arrays.asList(
                            new XsvTableFeature(1, 3, 4, file1Headers, file1Line1, "XSV_LOCATABLE_TEST_NAME"),
                            new XsvTableFeature(1, 3, 4, file1Headers, file1Line2, "XSV_LOCATABLE_TEST_NAME")
                    )
                }
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

    private void testDecodeHelper(final String filePath, final List<XsvTableFeature> expected) {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();
        if (xsvLocatableTableCodec.canDecode(filePath)) {
            try ( final FileInputStream fileInputStream = new FileInputStream(filePath)) {

                // Lots of scaffolding to do reading here:
                final AsciiLineReaderIterator    lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));
                final ArrayList<XsvTableFeature> output             = new ArrayList<>(expected.size());

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

    // Attempt to decode a malformed file:
    @Test(dataProvider = "provideForTestDecodeCharsetFailure",
            expectedExceptions = {UserException.MalformedFile.class})
    public void testDecodeCharsetFailure(final String filePath ) {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();
        if (xsvLocatableTableCodec.canDecode(filePath)) {
            try ( final FileInputStream fileInputStream = new FileInputStream(filePath)) {

                // Lots of scaffolding to do reading here:
                final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));

                // Read off the header:
                xsvLocatableTableCodec.readActualHeader(lineReaderIterator);

                // Read and decode the lines:
                while ( lineReaderIterator.hasNext() ) {
                    xsvLocatableTableCodec.decode(lineReaderIterator.next());
                }
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

    // decode
    @Test(dataProvider = "provideForTestDecode")
    public void testDecode(final String filePath, final List<XsvTableFeature> expected) {
        testDecodeHelper(filePath, expected);
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

    @Test
    public void testRenderSamFileHeaderFromNoPreamble() {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();
        final String filePath = TEST_FILE3;
        readHeaderOnly(xsvLocatableTableCodec, filePath);

        final SAMFileHeader emptyHeader = xsvLocatableTableCodec.renderSamFileHeader();

        Assert.assertEquals(emptyHeader.getComments().size(), 0);
        Assert.assertEquals(emptyHeader.getReadGroups().size(), 0);
        Assert.assertEquals(emptyHeader.getSequenceDictionary().size(), 0);
    }

    private List<String> readHeaderOnly(final XsvLocatableTableCodec xsvLocatableTableCodec, final String filePath) {
        List<String> header = null;

        if (xsvLocatableTableCodec.canDecode(filePath)) {

            try ( final FileInputStream fileInputStream = new FileInputStream(filePath)) {
                final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));
                header = xsvLocatableTableCodec.readActualHeader(lineReaderIterator);
            }
            catch ( final FileNotFoundException ex ) {
                throw new GATKException("Error - could not find test file: " + filePath, ex);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Error - IO problem with file " + filePath, ex);
            }
        }

        return header;
    }

    @Test
    public void testRenderSamFileHeaderFromSamFileHeaderPreamble() {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();
        final String filePath = TEST_FILE_SAMFILEHEADER;
        Assert.assertNotNull(readHeaderOnly(xsvLocatableTableCodec, filePath), "Header could not be decoded, but it should have been okay.");

        final SAMFileHeader populatedHeader = xsvLocatableTableCodec.renderSamFileHeader();

        Assert.assertEquals(populatedHeader.getSequenceDictionary().size(), 1);
        Assert.assertEquals(populatedHeader.getSequenceDictionary().getSequence(0).getSequenceLength(), 1000000);
        Assert.assertNotNull(populatedHeader.getReadGroup("GATKCopyNumber"));
        Assert.assertEquals(populatedHeader.getReadGroup("GATKCopyNumber").getSample(), "sample_cars");
        Assert.assertEquals(populatedHeader.getComments(), Arrays.asList("@CO\tCars, cars, and cars",
                "@CO\tNo \"family cars\" in this list."));
    }

    @Test
    public void testFalseCanDecodeFromMixedPreambles() {
        // Test that if we read a preamble that is mixed ("#" and "@"), we return a false for canDecode.
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec();
        final String filePath = TEST_FILE_MIXED_PREAMBLE;
        Assert.assertFalse(xsvLocatableTableCodec.canDecode(filePath));
    }

    @Test(expectedExceptions = UserException.BadInput.class, expectedExceptionsMessageRegExp = ".*Input did not contain any headers from the list.*")
    public void testTrueCanDecodeFromMissingColumn() {
        // Can decode should be true, but the parsing of the header should fail.
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec(Paths.get(TEST_FILE_SAMFILEHEADER_CONFIG_ERROR_NOTHING_FOUND));
        Assert.assertTrue(xsvLocatableTableCodec.canDecode(TEST_FILE_SAMFILEHEADER));
        final List<String> header = readHeader(xsvLocatableTableCodec, TEST_FILE_SAMFILEHEADER);
    }

    @DataProvider(name="simpleTestsWithMultipleColumnsInConfig")
    public Object[][] provideSimpleTestsWithMultipleColumnsInConfig() {
        return new Object[][]{
                {TEST_FILE_SAMFILEHEADER_CONFIG_MULTIPLE_COLUMNS, TEST_FILE_SAMFILEHEADER,
                        Arrays.asList("SECOND_XSV_NAME_chr", "SECOND_XSV_NAME_start", "SECOND_XSV_NAME_end")},
                {TEST_FILE_SAMFILEHEADER_CONFIG_MULTIPLE_COLUMNS_NO_NAME, TEST_FILE_SAMFILEHEADER,
                        Arrays.asList("chr", "start", "end")},
                {TEST_FILE_SAME_STARTEND_CONFIG, TEST_FILE_SAME_STARTEND,
                        Arrays.asList("XSV_LOCATABLE_TEST_NAME_chr", "XSV_LOCATABLE_TEST_NAME_Position", "XSV_LOCATABLE_TEST_NAME_Position")},
                {TEST_FILE_SAME_STARTEND_CONFIG_NO_NAME, TEST_FILE_SAME_STARTEND,
                        Arrays.asList("chr", "Position", "Position")}

        };
    }
    @Test(dataProvider = "simpleTestsWithMultipleColumnsInConfig")
    public void testDecodeMultipleChoiceHeaders(final String configFile, final String xsvFile, final List<String> locatableCols) {
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec(Paths.get(configFile));
        Assert.assertTrue(xsvLocatableTableCodec.canDecode(xsvFile));
        final List<String> header = readHeader(xsvLocatableTableCodec, xsvFile);
        Assert.assertNotNull(header);
        Assert.assertEquals(xsvLocatableTableCodec.getContigColumn(), locatableCols.get(0));
        Assert.assertEquals(xsvLocatableTableCodec.getStartColumn(), locatableCols.get(1));
        Assert.assertEquals(xsvLocatableTableCodec.getEndColumn(), locatableCols.get(2));
    }

    private List<String> readHeader(final XsvLocatableTableCodec xsvLocatableTableCodec, final String filePath) {
        List<String> header;
        try ( final InputStream inputStream = new FileInputStream(filePath)) {
            final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(inputStream));
            header = xsvLocatableTableCodec.readActualHeader(lineReaderIterator);
        }
        catch ( final FileNotFoundException ex ) {
            throw new GATKException("Error - could not find test file: " + filePath, ex);
        }
        catch ( final IOException ex ) {
            throw new GATKException("Error - IO problem with test file " + filePath, ex);
        }
        return header;
    }

    @DataProvider(name="contigNameErrors")
    public Object[][] provideContigNameErrors() {
        return new Object[][]{
                {TEST_RESOURCE_DIR + "xsv_locatable_test_samfileheader_error_contig_equals_start.config", TEST_FILE_SAMFILEHEADER},
                {TEST_RESOURCE_DIR + "xsv_locatable_test_samfileheader_error_contig_equals_end.config", TEST_FILE_SAMFILEHEADER}
        };
    }

    @Test(expectedExceptions = UserException.BadInput.class, expectedExceptionsMessageRegExp = ".*is the same as start or end column.*", dataProvider = "contigNameErrors")
    public void testBadContigColumnNames(final String configFile, final String xsvFile) {
        // Failure should happen when trying to get the header.
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec(Paths.get(configFile));
        Assert.assertTrue(xsvLocatableTableCodec.canDecode(xsvFile));
        final List<String> header = readHeader(xsvLocatableTableCodec, xsvFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNotFoundInList() {
        // Failure should happen when trying to get the header.
        final XsvLocatableTableCodec xsvLocatableTableCodec = new XsvLocatableTableCodec(Paths.get(TEST_FILE_SAMFILEHEADER_CONFIG_ERROR_NOTHING_FOUND));
        Assert.assertTrue(xsvLocatableTableCodec.canDecode(TEST_FILE_SAMFILEHEADER));
        final List<String> header = readHeader(xsvLocatableTableCodec, TEST_FILE_SAMFILEHEADER);
    }
}
