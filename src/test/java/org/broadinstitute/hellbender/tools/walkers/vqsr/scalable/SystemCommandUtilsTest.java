package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
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

    static void runH5Diff(final String expected, final String actual) {
        final ProcessController controller = ProcessController.getThreadLocal();

        // -r: Report mode. Print the differences.
        // --use-system-epsilon: Return a difference if and only if the difference between two data values exceeds
        //                       the system value for epsilon.
        final String[] command = new String[] { "h5diff", "-r", expected, actual };

        runProcessAndCaptureOutputInExceptionMessage(controller, command);
    }

    static void runDiff(final String expected, final String actual) {
        final ProcessController controller = ProcessController.getThreadLocal();
        final String[] command = new String[] { "diff", expected, actual };
        runProcessAndCaptureOutputInExceptionMessage(controller, command);
    }

    @Test(groups = {"python"}) // python environment is required to use h5diff
    public void testRunH5DiffIdenticalFile() {
        runH5Diff(String.format("%s/extract.AS.indel.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR),
                  String.format("%s/extract.AS.indel.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR));
    }

    @Test
    public void testRunDiffIdenticalFile() {
        runDiff(String.format("%s/extract.AS.indel.pos.vcf", EXPECTED_TEST_FILES_DIR),
                String.format("%s/extract.AS.indel.pos.vcf", EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class, groups = {"python"}) // python environment is required to use h5diff
    public void testRunH5DiffException() {
        runH5Diff(String.format("%s/extract.AS.indel.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR),
                  String.format("%s/extract.AS.snp.pos.annot.hdf5", EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class)
    public void testRunDiffException() {
        runDiff(String.format("%s/extract.AS.indel.pos.vcf", EXPECTED_TEST_FILES_DIR),
                String.format("%s/extract.AS.snp.pos.vcf", EXPECTED_TEST_FILES_DIR));
    }

    @Test(expectedExceptions = AssertionError.class)
    public void testRunDiffNoSuchFile() {
        runDiff(String.format("%s/blahblah", EXPECTED_TEST_FILES_DIR),
                String.format("%s/extract.AS.snp.pos.vcf", EXPECTED_TEST_FILES_DIR));
    }
}
