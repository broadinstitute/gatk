package org.broadinstitute.hellbender.tools.picard.sam;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public final class CreateSequenceDictionaryIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/CreateSequenceDictionary");
    public static File BASIC_FASTA = new File(TEST_DATA_DIR, "basic.fasta");
    public static File DUPLICATE_FASTA = new File(TEST_DATA_DIR, "duplicate_sequence_names.fasta");

    @Override
    public String getTestedClassName() {
        return CreateSequenceDictionary.class.getSimpleName();
    }

    @Test
    public void testBasic() throws Exception {
        final File outputDict = BaseTest.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.delete();
        final String[] argv = {
                "--reference", BASIC_FASTA.getAbsolutePath(),
                "--output", outputDict.getAbsolutePath(),
                "--URI", BASIC_FASTA.getName()
        };
        runCommandLine(argv);
        IntegrationTestSpec.assertEqualTextFiles(outputDict, new File(BASIC_FASTA.getAbsolutePath() + ".dict"));
    }

    /**
     * Should throw an exception because sequence names are not unique.
     */
    @Test(expectedExceptions = {UserException.MalformedFile.class})
    public void testNonUniqueSequenceName() throws Exception {
        final File outputDict = BaseTest.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.delete();
        final String[] argv = {
                "--reference", DUPLICATE_FASTA.getAbsolutePath(),
                "--output", outputDict.getAbsolutePath(),
        };
        runCommandLine(argv);
        Assert.fail("Exception should have been thrown.");
    }

    // Should throw an exception because no reference file was specified
    @Test(expectedExceptions = {CommandLineException.class})
    public void testNoReferenceSpecified() throws Exception {
        final File output = BaseTest.createTempFile("TestOutput", ".out");
        final String[] argv = {
                "CreateSequenceDictionary",
                "--output", output.getAbsolutePath()
        };
        runCommandLine(argv);
    }
}
