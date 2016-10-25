package org.broadinstitute.hellbender.tools.picard.vcf;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class UpdateVcfSequenceDictionaryIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testUpdateVcfSequenceDictionary() throws Exception {
        final File expectedFile = new File(getTestDataDir(), "count_variants_withSequenceDict.vcf");
        final File tempFile = createTempFile("UpdateVcfSequenceDictionary", ".vcf");
        final ArgumentsBuilder ab = new ArgumentsBuilder()
                .addArgument("SEQUENCE_DICTIONARY", expectedFile.getAbsolutePath())
                .addArgument("INPUT", expectedFile.getAbsolutePath())
                .addArgument("OUTPUT", tempFile.getAbsolutePath());
        runCommandLine(ab);
        // the only part that we consider a "comment" is the VCF format line
        // the rest of the file should be exactly the same
        IntegrationTestSpec.assertEqualTextFiles(tempFile, expectedFile, "##fileformat=VCFv");
    }
}