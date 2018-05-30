package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Note that realignment filtering requires knowledge of the entire genome, since an alignment artifact at one locus
 * could come from anywhere.  Proper testing requires the BWA mem index image of the entire reference, which is
 * 7 GB and too large for a git repo.  Therefore, this test is absolutely no substitute for testing the tool elsewhere.
 */
public class FilterAlignmentArtifactsIntegrationTest extends CommandLineProgramTest {

    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final String DREAM_VCFS_DIR = toolsTestDir + "mutect/dream/vcfs/";
    private static final String DREAM_MASKS_DIR = toolsTestDir + "mutect/dream/masks/";

    // this is an image made just from bases 20:2,000,000-3,000,000 of the b37 reference
    private static final String BWA_MEM_INDEX_IMAGE = toolsTestDir + "mutect/human_g1k_v37.20.2000000.3000000.img";

    @Test
    public void testDreamTruthData() {
        final File tumorBam = new File(DREAM_BAMS_DIR, "tumor_4.bam");
        final File truthVcf = new File(DREAM_VCFS_DIR, "sample_4.vcf");
        final File mask = new File(DREAM_MASKS_DIR, "mask4.list");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final String[] args = {
                "-I", tumorBam.getAbsolutePath(),
                "-V", truthVcf.getAbsolutePath(),
                "-L", "20",
                "--bwa-mem-index-image", BWA_MEM_INDEX_IMAGE,
                "-XL", mask.getAbsolutePath(),
                "-O", filteredVcf.getAbsolutePath()
        };

        runCommandLine(args);
    }
}