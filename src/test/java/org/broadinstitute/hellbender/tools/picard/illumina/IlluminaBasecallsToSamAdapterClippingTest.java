package org.broadinstitute.hellbender.tools.picard.illumina;

import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintWriter;

import static htsjdk.samtools.ReservedTagConstants.XT;
import static htsjdk.samtools.SamReaderFactory.makeDefault;
import static java.lang.System.out;
import static org.testng.Assert.assertEquals;

/**
 * Run IlluminaBasecallsToSam on a sample tests, then sanity-check the generated SAM file
 */
public class IlluminaBasecallsToSamAdapterClippingTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/illumina/125T125T");
    private static final File BASECALLS_DIR = new File(TEST_DATA_DIR, "Data/Intensities/BaseCalls");
    private static final String RUN_BARCODE = "305PJAAXX080716";

    public String getTestedClassName() {
        return IlluminaBasecallsToSam.class.getSimpleName();
    }

    /**
     * Run IlluminaBasecallsToSam on a few test cases, and verify that results agree with hand-checked expectation.
     */
    @Test(dataProvider = "data")
    public void testBasic(final String LANE, final String readStructure) throws Exception {
        final File samFile = File.createTempFile(LANE + ".illuminaBasecallsToSam", ".sam");
        samFile.deleteOnExit();
        final File libraryParams = File.createTempFile(LANE + ".illuminaBasecallsToSam", ".params");
        libraryParams.deleteOnExit();
        final PrintWriter writer = new PrintWriter(libraryParams);
        writer.println("SAMPLE_ALIAS\tLIBRARY_NAME\tOUTPUT");
        writer.println("myalias\tfoo\t" + samFile.getAbsolutePath());
        writer.close();

        final String[] illuminaArgv = {
                "--BASECALLS_DIR", BASECALLS_DIR.getAbsolutePath(),
                "--LANE", LANE,
                "--RUN_BARCODE", RUN_BARCODE,
                "--READ_STRUCTURE", readStructure,
                "--LIBRARY_PARAMS", libraryParams.getAbsolutePath()
        };
        runCommandLine(illuminaArgv);

        // Read the file and confirm it contains what is expected
        final SamReader samReader = makeDefault().open(samFile);

        // look for clipped adaptor attribute in lane 3 PE (2) and in lane 6 (1) non-PE
        int count = 0;
        for (final SAMRecord record : samReader) {
            if (record.getIntegerAttribute(XT) != null) {
                count++;
                if ((count == 1 || count == 2) && LANE.equals("2")) {
                    assertEquals(114, (int) record.getIntegerAttribute(XT));
                } else if (count == 1 || count == 2 && LANE.equals("1")) {
                    assertEquals(68, (int) record.getIntegerAttribute(XT));
                }
            }
        }
        samReader.close();
    }

    @DataProvider(name = "data")
    private Object[][] getIlluminaBasecallsToSamTestData() {
        return new Object[][]{
                {"1", "125T125T"},
                {"2", "125T125T"},
        };
    }
}