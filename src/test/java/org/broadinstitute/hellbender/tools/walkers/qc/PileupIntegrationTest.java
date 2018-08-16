package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class PileupIntegrationTest extends CommandLineProgramTest {

    private final String TEST_OUTPUT_DIRECTORY = getToolTestDataDir().toLowerCase() + "/";

    @Test
    public void testSimplePileup() throws IOException {
        // GATK 3.5 code have a the last line with a REDUCE RESULT that was removed in this implementation
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -L 20:9999900-10000000" +
                " -R " + b37_reference_20_21 +
                " -I " + NA12878_20_21_WGS_bam +
                " -O %s",
            Arrays.asList(TEST_OUTPUT_DIRECTORY + "expectedSimplePileup.txt")
        );
        testSpec.executeTest("testSimplePileup", this);
    }

    @Test
    public void testVerbosePileup() throws IOException {
        // GATK 3.5 code have a the last line with a REDUCE RESULT that was removed in this implementation
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -L 20:9999990-10000000" +
                " -verbose " +
                " -R " + b37_reference_20_21 +
                " -I " + NA12878_20_21_WGS_bam +
                " -O %s",
            Arrays.asList(TEST_OUTPUT_DIRECTORY + "expectedVerbosePileup.txt")
        );
        testSpec.executeTest("testVerbosePileup", this);
    }

    @Test
    public void testFeaturesPileup() throws IOException {
        // GATK 3.5 code have a the last line with a REDUCE RESULT that was removed in this implementation
        // GATK 3.5 code have ROD instead of Feature(s) and the source of a VariantContext is set to metadata, thus the output was modified with the following command
        // awk '{gsub("metadata", "Unknown", $0); gsub("\\[ROD: ", "\[Feature(s): ", $0); print $0}'
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -L 20:10000092-10000112" +
                " -R " + b37_reference_20_21 +
                " -I " + NA12878_20_21_WGS_bam +
                " -metadata " + dbsnp_138_b37_20_21_vcf +
                " -O %s",
            Arrays.asList(TEST_OUTPUT_DIRECTORY + "expectedFeaturesPileup.txt")
        );
        testSpec.executeTest("testFeaturesPileup", this);
    }

    @Test
    public void testInsertLengthPileup() throws  Exception {
        // GATK 3.5 code have a the last line with a REDUCE RESULT that was removed in this implementation
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " -L 20:10000092-10000112" +
                        " -R " + b37_reference_20_21 +
                        " -I " + NA12878_20_21_WGS_bam +
                        " -output-insert-length " +
                        " -O %s",
                Arrays.asList(TEST_OUTPUT_DIRECTORY + "expectedInsertLengthPileup.txt")
        );
        testSpec.executeTest("testInsertLengthPileup", this);
    }

    @Test(expectedExceptions = UserException.CouldNotCreateOutputFile.class)
    public void testInvalidOutputFile() throws IOException {
        // GATK 3.5 code have a the last line with a REDUCE RESULT that was removed in this implementation
        final String[] args = new String[] {
                "-L" , "20:9999900-10000000",
                "-R" , b37_reference_20_21,
                "-I" , NA12878_20_21_WGS_bam,
                "-O" , "/foobar/foobar"
        };
        runCommandLine(args);
    }


}
