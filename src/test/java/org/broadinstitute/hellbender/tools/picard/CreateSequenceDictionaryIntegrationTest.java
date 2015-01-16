package org.broadinstitute.hellbender.tools.picard;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public class CreateSequenceDictionaryIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/CreateSequenceDictionary");
    public static File BASIC_FASTA = new File(TEST_DATA_DIR, "basic.fasta");
    public static File DUPLICATE_FASTA = new File(TEST_DATA_DIR, "duplicate_sequence_names.fasta");

    public String getCommandLineProgramName() {
        return CreateSequenceDictionary.class.getSimpleName();
    }

    @Test
    public void testBasic() throws Exception {
        final File outputDict = File.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.delete();
        outputDict.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + BASIC_FASTA,
                "OUTPUT=" + outputDict,
        };
        Assert.assertEquals(runCommandLine(new ArgumentsBuilder(argv).getArgsArray()), null);
    }

    /**
     * Should throw an exception because sequence names are not unique.
     */
    @Test(expectedExceptions = {UserException.MalformedFile.class})
    public void testNonUniqueSequenceName() throws Exception {
        final File outputDict = File.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.delete();
        outputDict.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + DUPLICATE_FASTA,
                "OUTPUT=" + outputDict,
        };
        Assert.assertEquals(runCommandLine(new ArgumentsBuilder(argv).getArgsArray()), null);
        Assert.fail("Exception should have been thrown.");
    }
}
