package org.broadinstitute.hellbender.tools.exome.convertbed;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public class ConvertBedToTargetFileIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File TEST_BED_INPUT = new File(TEST_DIR, "testbedconversion.bed");
    private static final File TEST_BED_INPUT_MORE_ANNOTATIONS = new File(TEST_DIR, "testbedconversion-more-annotations.bed");
    private static final File TARGET_FILE_GT = new File(TEST_DIR, "test_target_file.txt");
    private static final File TEST_TARGET_OUTPUT = createTempFile("test_target_file",".txt");

    @Test
    public void testConvertBedToTargetFile() {
        final String[] arguments = {
                "-" + ConvertBedToTargetFile.BED_INPUT_SHORT_NAME, TEST_BED_INPUT.getAbsolutePath(),
                "-" + ConvertBedToTargetFile.TARGET_OUTPUT_SHORT_NAME, TEST_TARGET_OUTPUT.getAbsolutePath(),
        };
        runCommandLine(arguments);
        Assert.assertTrue(TEST_TARGET_OUTPUT.exists());
        final List<Target> outputTestTargets = TargetTableReader.readTargetFile(TEST_TARGET_OUTPUT);
        Assert.assertNotNull(outputTestTargets);
        final List<Target> gtTargets = TargetTableReader.readTargetFile(TARGET_FILE_GT);
        Assert.assertEquals(outputTestTargets.size(), gtTargets.size());
        for(int i=0; i<gtTargets.size(); i++){
            Assert.assertEquals(outputTestTargets.get(i), gtTargets.get(i));
        }
    }

    @Test
    public void testConvertBedToTargetFileMoreAnnotations() {
        final String[] arguments = {
                "-" + ConvertBedToTargetFile.BED_INPUT_SHORT_NAME, TEST_BED_INPUT_MORE_ANNOTATIONS.getAbsolutePath(),
                "-" + ConvertBedToTargetFile.TARGET_OUTPUT_SHORT_NAME, TEST_TARGET_OUTPUT.getAbsolutePath(),
        };
        runCommandLine(arguments);
        Assert.assertTrue(TEST_TARGET_OUTPUT.exists());
        final List<Target> outputTestTargets = TargetTableReader.readTargetFile(TEST_TARGET_OUTPUT);
        Assert.assertNotNull(outputTestTargets);
        final List<Target> gtTargets = TargetTableReader.readTargetFile(TARGET_FILE_GT);
        Assert.assertEquals(outputTestTargets.size(), gtTargets.size());
        for (int i = 0; i < gtTargets.size(); i++) {
            Assert.assertEquals(outputTestTargets.get(i), gtTargets.get(i));
        }
    }
}