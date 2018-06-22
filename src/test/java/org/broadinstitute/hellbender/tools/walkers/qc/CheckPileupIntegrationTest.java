package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
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
    private static final String TEST_OUTPUT_DIRECTORY = toolsTestDir + "walkers/qc/pileup/";

    /**
     * This test runs on a basic pileup obtained with samtools (version 1.3.1) and options -B --min-BQ 0
     */
    @Test
    public void testBasicPileup() throws IOException {
        final File emptyTemp = createTempFile("empty", "txt");
        emptyTemp.createNewFile();
        // pileup was generated with "samtools -f hg19MiniReference -B --min-BQ 0 reads_data_source_test1.bam"
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
            " -R " + hg19MiniReference +
            " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
            " -pileup " +  TEST_OUTPUT_DIRECTORY + "reads_data_source_test1.samtools.pileup" +
            " -O %s", Arrays.asList(emptyTemp.toString()));

        testSpec.executeTest("testBasicPileup", this);
    }

    /**
     * This test runs on a basic pileup obtained with samtools (version 1.3.1) and options --min-BQ 0
     * BAQ quality recalibration is activated and should be different
     */
    @Test
    public void testBAQPileup() throws IOException {
        // pileup was generated with "samtools -f hg19MiniReference --min-BQ 0 reads_data_source_test1.bam"
        IntegrationTestSpec testSpec = new IntegrationTestSpec(
                " --continue-after-error " +
                " -R " + hg19MiniReference +
                " -I " + TEST_DATA_DIRECTORY + "reads_data_source_test1.bam" +
                " -pileup " +  TEST_OUTPUT_DIRECTORY + "reads_data_source_test1.samtools.baq.pileup" +
                " -O %s", Arrays.asList(TEST_OUTPUT_DIRECTORY + "reads_data_source_test1.samtools.baq.pileup.diff"));

        testSpec.executeTest("testBasicPileup", this);
    }

}
