package org.broadinstitute.hellbender.utils.bwa;

import org.broadinstitute.hellbender.BwaMemIntegrationTest;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Test if image can be successfully created.
 */
public class BwaMemIndexImageCreatorTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("src/test/resources/large/");
    private static final File testReferenceFasta = new File(TEST_DATA_DIR, "human_g1k_v37.20.21.fasta");

    @Override
    public String getTestedClassName() {
        return BwaMemIndexImageCreator.class.getSimpleName();
    }

    @Test
    public void testImageFileGeneration() throws Exception {

        final File tempImage = BaseTest.createTempFile("tempBwaMemIndexImage", ".img");
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--input", testReferenceFasta.getAbsolutePath(),
                "--output", tempImage.getAbsolutePath()));
        runCommandLine(args);

        // pigg-backing on the existing integration test
        final BwaMemIndex index = new BwaMemIndex(tempImage.getAbsolutePath());
        BwaMemIntegrationTest.singleReadAlingment(index);
        BwaMemIntegrationTest.chimericContigAlignment(index);

        index.close();
    }

}
