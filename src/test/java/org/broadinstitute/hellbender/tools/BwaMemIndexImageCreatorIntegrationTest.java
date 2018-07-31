package org.broadinstitute.hellbender.tools;

import org.broadinstitute.hellbender.BwaMemTestUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Test if image can be successfully created.
 */
public class BwaMemIndexImageCreatorIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("src/test/resources/large/");
    private static final File testReferenceFasta = new File(TEST_DATA_DIR, "human_g1k_v37.20.21.fasta");

    @Test
    public void testImageFileGeneration() throws Exception {

        final File tempImage = GATKBaseTest.createTempFile("tempBwaMemIndexImage", ".img");
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--" + StandardArgumentDefinitions.INPUT_LONG_NAME, testReferenceFasta.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME, tempImage.getAbsolutePath()));
        runCommandLine(args);

        // piggy-backing on the existing integration test
        try( final BwaMemIndex index = new BwaMemIndex(tempImage.getAbsolutePath()) ){
            BwaMemTestUtils.assertCorrectSingleReadAlignment(index);
            BwaMemTestUtils.assertCorrectChimericContigAlignment(index);
        } finally {
            try { tempImage.delete(); } catch (final Throwable t) {};
        }
    }
}
