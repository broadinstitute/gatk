package org.broadinstitute.hellbender.tools.walkers.longreads;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;

public class BamCustomIndexerPrototypeIntegrationTest extends CommandLineProgramTest {

    public static final String SMALL_TEST_BAM = "src/test/resources/chrM_and_chr20_subset.corrected.bam";


    @Test
    public void testSmallBam() {
        final File outputBam = createTempFile("BamCustomIndexerPrototypeIntegrationTest_testSmallBam", ".bam");
        final File outputIndex = createTempFile("BamCustomIndexerPrototypeIntegrationTest_testSmallBam", ".index");

        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addInput(new File(SMALL_TEST_BAM));
        args.addOutput(outputBam);
        args.addArgument("output-index", outputIndex.getAbsolutePath());

        runCommandLine(args);
    }
}
