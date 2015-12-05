package org.broadinstitute.hellbender.tools.picard.sam;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;

/**
 * ValidateSamFile is a thin wrapper around {@link htsjdk.samtools.SamFileValidator} which is thoroughly tested in HTSJDK.
 * Thus we don't have much to test here.
 */
public final class ValidateSamFileIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/ValidateSamFile");

    @Override
    public String getTestedClassName() {
        return ValidateSamFile.class.getSimpleName();
    }

    @DataProvider(name="tooltestingData")
    public Object[][] testingData() {
        return new Object[][]{
          {"valid.sam", "SUMMARY", null, true},
          {"valid.sam", "VERBOSE", null, true},
          {"valid.bam", "VERBOSE", null, true},
          {"valid.bam", "VERBOSE", "valid.fasta", false}, // the NM tags in this sam file don't match the reference deltas
          {"valid.cram", "VERBOSE", "valid.fasta", true}, // the NM tags in the CRAM file do match the reference (ty samtools)
          {"invalid_coord_sort_order.sam", "SUMMARY", null, false},
          {"invalid_coord_sort_order.sam", "SUMMARY", "invalid_coord_sort_order.fasta", false},
          {"invalid_coord_sort_order.cram", "SUMMARY", "invalid_coord_sort_order.fasta", false},
        };
    }

    @Test(dataProvider="tooltestingData")
    public void testNoOutputFile(final String input, final String mode,  final String referenceName, final boolean expectedValidity) throws Exception {
        final File testFile = new File(TEST_DATA_DIR, input);
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(testFile.getPath());
        args.add("--MODE");
        args.add(mode);
        if (null != referenceName) {
            final File refFile = new File(TEST_DATA_DIR, referenceName);
            args.add("--R ");
            args.add(refFile.getAbsoluteFile());
        };
        Assert.assertEquals(runCommandLine(args.getArgsArray()), expectedValidity);
    }

}
