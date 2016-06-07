package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Run validating pileup across a set of core data as proof of the integrity of the GATK core.
 *
 * Tests only single-sample mpileups (no consensus)
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class CheckPileupIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String TEST_OUTPUT_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/tools/walkers/qc/pileup/";

    /**
     * This test runs on a basic pileup obtained with samtools (version 1.1) and options -B --min-BQ 0
     */
    @Test
    public void testBasicPileup() throws IOException {
        final File emptyTemp = createTempFile("empty", "txt");
        emptyTemp.createNewFile();
        // pileup was generated with "samtools -f hg19MiniReference -B --min-BQ 0 reads_data_source_test1.bam"
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -R " + hg19MiniReference +
            " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
            " -pileup " +  TEST_OUTPUT_DIRECTORY + "reads_data_source_test1.samtools.samp" +
            " -O %s", Arrays.asList(emptyTemp.toString()));

        testSpec.executeTest("testBasicPileup", this);
    }
}
