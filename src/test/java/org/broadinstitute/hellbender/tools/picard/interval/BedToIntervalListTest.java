package org.broadinstitute.hellbender.tools.picard.interval;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * @author nhomer
 */
public final class BedToIntervalListTest extends CommandLineProgramTest {

    private final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/interval/BedToIntervalListTest");

    private void doTest(final String inputBed, final String header) throws IOException {
        final File outputFile = BaseTest.createTempFile("bed_to_interval_list_test.", ".interval_list");
        final BedToIntervalList program = new BedToIntervalList();
        final File inputBedFile = new File(TEST_DATA_DIR, inputBed);
        program.INPUT = inputBedFile;
        program.SEQUENCE_DICTIONARY = new File(TEST_DATA_DIR, header);
        program.OUTPUT = outputFile;
        program.doWork();

        // Assert they are equal
        SAMException unexpectedException = null;
        try {
            IOUtil.assertFilesEqual(new File(inputBedFile.getAbsolutePath() + ".interval_list"), outputFile);
        } catch (final SAMException e) {
            unexpectedException = e;
        }
        Assert.assertEquals(unexpectedException, null);
    }

    @Test(dataProvider = "testBedToIntervalListDataProvider")
    public void testBedToIntervalList(final String inputBed) throws IOException {
        doTest(inputBed, "header.sam");
    }

    @Test(dataProvider = "testBedToIntervalListOutOfBoundsDataProvider", expectedExceptions = GATKException.class)
    public void testBedToIntervalListOutOfBounds(final String inputBed) throws IOException {
        doTest(inputBed, "header.sam");
    }

    @DataProvider
    public Object[][] testBedToIntervalListDataProvider() {
        return new Object[][]{
                {"simple.bed"},
                {"overlapping.bed"},
                {"extended.bed"},
                {"one_base_interval.bed"},
                {"zero_base_interval.bed"}
        };
    }

    @DataProvider
    public Object[][] testBedToIntervalListOutOfBoundsDataProvider() {
        return new Object[][]{
                {"end_after_chr.bed"},
                {"end_before_chr.bed"},
                {"missing_chr.bed"},
                {"start_after_chr.bed"},
                {"start_before_chr.bed"},
                {"off_by_one_interval.bed"}
        };
    }
}
