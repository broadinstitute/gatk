package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.read.testers.CleanSamTester;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class CleanSamIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/CleanSam");

    public String getTestedClassName() {
        return CleanSam.class.getSimpleName();
    }

    @Test(dataProvider = "testCleanSamDataProvider")
    public void testCleanSam(final String samFile, final String expectedCigar) throws IOException {
        final File inputFile = new File(TEST_DATA_DIR, samFile);
        final File cleanedFile = File.createTempFile(samFile + ".", ".sam");
        cleanedFile.deleteOnExit();
        final String[] args = new String[]{
                "--INPUT", inputFile.getAbsolutePath(),
                "--OUTPUT", cleanedFile.getAbsolutePath()
        };

        runCommandLine(args);

        final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);
        validator.setIgnoreWarnings(true);
        validator.setVerbose(true, 1000);
        validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.MISSING_READ_GROUP));
        SamReader samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(cleanedFile);
        final SAMRecord rec = samReader.iterator().next();
        samReader.close();
        Assert.assertEquals(rec.getCigarString(), expectedCigar);
        samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(cleanedFile);
        final boolean validated = validator.validateSamFileVerbose(samReader, null);
        samReader.close();
        Assert.assertTrue(validated, "ValidateSamFile failed");
    }

    @DataProvider(name = "testCleanSamDataProvider")
    public Object[][] testCleanSamDataProvider() {
        return new Object[][]{
                {"simple_fits.sam", "100M"},
                {"simple_overhang.sam", "99M1S"},
                {"fits_with_deletion.sam", "91M2D9M"},
                {"overhang_with_deletion.sam", "91M2D8M1S"},
                {"trailing_insertion.sam", "99M1I"},
                {"long_trailing_insertion.sam", "90M10I"},
        };
    }

    //identical test case using the SamFileTester to generate that SAM file on the fly
    @Test(dataProvider = "testCleanSamTesterDataProvider")
    public void testCleanSamTester(final String originalCigar, final String expectedCigar, final int defaultChromosomeLength, final int alignStart) throws IOException {
        final CleanSamTester cleanSamTester = new CleanSamTester(expectedCigar, 100, defaultChromosomeLength);
        // NB: this will add in the mate cigar, when enabled in SamPairUtil, for additional validation
        cleanSamTester.addMappedPair(0, alignStart, alignStart, false, false, originalCigar, originalCigar, false, 50);
        cleanSamTester.runTest();
    }

    @DataProvider(name = "testCleanSamTesterDataProvider")
    public Object[][] testCleanSamTesterDataProvider() {
        return new Object[][]{
                {"100M", "100M", 101, 2}, // simple_fits.sam
                {"100M", "99M1S", 101, 3}, // simple_overhang.sam
                {"91M2D9M", "91M2D9M", 102, 1}, // fits_with_deletion.sam
                {"91M2D9M", "91M2D8M1S", 101, 1}, // overhang_with_deletion.sam
                {"99M1I", "99M1I", 101, 3}, // trailing_insertion.sam
                {"90M10I", "90M10I", 101, 3} // long_trailing_insertion.sam
        };
    }
}
