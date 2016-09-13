package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Created by davidben on 9/1/16.
 */
public class Mutect2IntegrationTest extends CommandLineProgramTest {
    private static final String TEST_FILES_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/";

    private static final File CCLE_MICRO_TUMOR_BAM = new File(TEST_FILES_DIR, "HCC1143.cghub.ccle.micro.bam");
    private static final File CCLE_MICRO_NORMAL_BAM = new File(TEST_FILES_DIR, "HCC1143_BL.cghub.ccle.micro.bam");
    private static final File CCLE_MICRO_INTERVALS_FILE = new File(TEST_FILES_DIR, "HCC1143.cghub.ccle.micro.intervals");

    //TODO: this requires Broad server; repalce with the mini reference in large test files directory
    private static final File hg19_REFERENCE = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");

    //TODO: note this test works if you run it locally and mount /seq/references
    //TODO: it's commented out because Travis will fail, but it's useful to have for development until
    //TODO: real integration tests come
    //TODO: currently, results look reasonable on what is admittedly a pretty tiny test
    @Test(enabled = false)
    public void testRunMutect2() throws Exception {
        Utils.resetRandomGenerator();
        final File output = createTempFile("HCC1143_ccle", ".vcf");

        final String[] args = {
                "-I", CCLE_MICRO_TUMOR_BAM.getAbsolutePath(),
                "-tumor", "HCC1143",
                "-I", CCLE_MICRO_NORMAL_BAM.getAbsolutePath(),
                "-normal", "HCC1143 BL",
                "-R", hg19_REFERENCE.getAbsolutePath(),
                "-L", CCLE_MICRO_INTERVALS_FILE.getAbsolutePath(),
                "-O", output.getAbsolutePath()
        };

        runCommandLine(args);
    }

}