package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

public final class SystemCommandUtilsTest extends GATKBaseTest {

    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/extract");
    private static final File EXPECTED_TEST_FILES_DIR = new File(TEST_FILES_DIR, "expected");

    // this method is duplicated in the other integration-test classes in this package
    static void runSystemCommand(final String command) {
        logger.debug(String.format("Testing command: %s", command));
        try {
            final ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", command).redirectErrorStream(true);
            final Process process = processBuilder.start();

            final BufferedReader stdInReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String stdInLine;
            while ((stdInLine = stdInReader.readLine()) != null) {
                Assert.fail(String.format("The command \"%s\" resulted in: %s", command, stdInLine));
            }
            stdInReader.close();

        } catch (final IOException e) {
            throw new GATKException.ShouldNeverReachHereException(e.getMessage());
        }
    }

    @Test(groups = {"python"}) // python environment is required to use h5diff
    public void testRunSystemCommand() {
        runSystemCommand(String.format("h5diff %s/extract.AS.indel.pos.annot.hdf5 %s/extract.AS.indel.pos.annot.hdf5",
                EXPECTED_TEST_FILES_DIR, EXPECTED_TEST_FILES_DIR));
        runSystemCommand(String.format("diff %s/extract.AS.indel.pos.vcf %s/extract.AS.indel.pos.vcf",
                EXPECTED_TEST_FILES_DIR, EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class, groups = {"python"}) // python environment is required to use h5diff
    public void testRunSystemCommandH5diffException() {
        runSystemCommand(String.format("h5diff %s/extract.AS.indel.pos.annot.hdf5 %s/extract.AS.snp.pos.annot.hdf5",
                EXPECTED_TEST_FILES_DIR, EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class)
    public void testRunSystemCommandDiffException() {
        runSystemCommand(String.format("diff %s/extract.AS.indel.pos.vcf %s/extract.AS.snp.pos.vcf",
                EXPECTED_TEST_FILES_DIR, EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class)
    public void testRunSystemCommandDiffNoSuchFileException() {
        runSystemCommand(String.format("diff %s/blahblah %s/extract.AS.snp.pos.vcf",
                EXPECTED_TEST_FILES_DIR, EXPECTED_TEST_FILES_DIR));
    }
}
