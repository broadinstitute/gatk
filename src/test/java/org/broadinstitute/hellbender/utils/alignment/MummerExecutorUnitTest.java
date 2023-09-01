package org.broadinstitute.hellbender.utils.alignment;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class MummerExecutorUnitTest extends GATKBaseTest {

    private final String COMPARE_REFERENCES_TEST_FILES_DIRECTORY = toolsTestDir + "/reference/CompareReferences/";

    @Test
    public void testExecuteMummer() throws IOException {
        File fasta1 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta");
        File fasta2 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_chr2multiplesnps.fasta");
        File expectedOutputDir = new File(getToolTestDataDir());
        File actualOutputDir = IOUtils.createTempDir("testMummer");
        MummerExecutor exec = new MummerExecutor();

        exec.executeMummer(fasta1, fasta2, actualOutputDir);
        File expectedOutput = new File(expectedOutputDir, "expected.snps_output.snps");
        File actualOutput = new File(actualOutputDir, "snps_output.snps");
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @DataProvider
    public Object[][] dataForTestPrepareMUMmerExecutionDirectory() {
        return new Object[][] {
                { MummerExecutor.MUMMER_X86_64_MAC_BINARIES_ZIPFILE },
                { MummerExecutor.MUMMER_X86_64_LINUX_BINARIES_ZIPFILE }
        };
    }

    @Test(dataProvider = "dataForTestPrepareMUMmerExecutionDirectory")
    public void testPrepareMUMmerExecutionDirectory(final String mummerZipToExtract) {
        final Resource mummerDistribution = new Resource(mummerZipToExtract, getClass());
        final File extractedMUMmerDir = MummerExecutor.prepareMUMmerExecutionDirectory(mummerDistribution);
        checkExtractedMummerDistribution(extractedMUMmerDir);
    }

    @Test
    public void testPrepareMUMmerExecutionDirectoryWithDefaultConstructor(){
        MummerExecutor exec = new MummerExecutor();
        File executableDirectory = exec.getMummerExecutableDirectory();
        checkExtractedMummerDistribution(executableDirectory);
    }

    private void checkExtractedMummerDistribution(final File extractedMUMmerDir) {
        Assert.assertTrue(new File(extractedMUMmerDir, "mummer").exists());
        Assert.assertTrue(new File(extractedMUMmerDir, "nucmer").exists());
        Assert.assertTrue(new File(extractedMUMmerDir, "delta-filter").exists());
        Assert.assertTrue(new File(extractedMUMmerDir, "show-snps").exists());

        for(File file : extractedMUMmerDir.listFiles()){
            Assert.assertTrue(file.getTotalSpace() > 0);
            Assert.assertTrue(file.canExecute());
        }
    }
}
