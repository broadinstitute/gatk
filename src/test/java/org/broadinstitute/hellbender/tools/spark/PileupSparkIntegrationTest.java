package org.broadinstitute.hellbender.tools.spark;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.List;

public final class PileupSparkIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "walkers/qc/pileup");

    @Test
    public void testSimplePileup() throws Exception {
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        out.deleteOnExit();
        run(out, false);
        File expected = new File(TEST_DATA_DIR, "expectedSimplePileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    @Test
    public void testSimplePileupWithShuffle() throws Exception {
        final File out = File.createTempFile("out", ".txt");
        out.delete();
        out.deleteOnExit();
        run(out, true);
        File expected = new File(TEST_DATA_DIR, "expectedSimplePileup.txt");
        IntegrationTestSpec.assertEqualTextFiles(new File(out, "part-00000"), expected);
    }

    private void run(File out, boolean shuffle) throws Exception {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--input");
        args.add(NA12878_20_21_WGS_bam);
        args.add("--output");
        args.add(out.getAbsolutePath());
        args.add("--reference");
        args.add(b37_reference_20_21);
        args.add("-L 20:9999900-10000000");
        if (shuffle) {
            args.add("--shuffle");
        }
        this.runCommandLine(args.getArgsArray());
    }

}