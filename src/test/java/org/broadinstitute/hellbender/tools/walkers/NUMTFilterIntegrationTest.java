package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

public final class NUMTFilterIntegrationTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIRECTORY = "/Users/ahaessly/code/gatk3/mito/sample_01C07450/";
    private static final String TEST_OUTPUT_DIRECTORY = "/Users/ahaessly/code/gatk3/mito/sample_01C07450/";

//    File MITOREF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");

    @Test
    public void testNUMTFilter() {
        final String[] args = {
                "-R", toolsTestDir + "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta",
                "-I", TEST_DATA_DIRECTORY + "01C07450.realigned.bam",
                "-V", TEST_DATA_DIRECTORY + "output-nofilter.vcf",
                "-O", TEST_OUTPUT_DIRECTORY + "NUMTFilter.bam"
        };
        runCommandLine(args);
    }

}
