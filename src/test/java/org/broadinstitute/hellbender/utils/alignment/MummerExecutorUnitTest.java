package org.broadinstitute.hellbender.utils.alignment;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class MummerExecutorUnitTest extends GATKBaseTest {

    private final String COMPARE_REFERENCES_TEST_FILES_DIRECTORY = toolsTestDir + "/reference/CompareReferences/";

    // need intel build for MUMmer
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

    @Test
    public void testPrepareMUMmerExecutionDirectory(){
        MummerExecutor exec = new MummerExecutor();
        File executableDirectory = exec.getMummerExecutableDirectory();
        File[] executables = executableDirectory.listFiles();
        Assert.assertEquals(executables.length, 4);

        for(File file : executables){
            Assert.assertTrue(file.getTotalSpace() > 0);
            Assert.assertTrue(file.canExecute());
            Assert.assertTrue(file.getName().equals("mummer") || file.getName().equals("nucmer") ||
                    file.getName().equals("delta-filter") || file.getName().equals("show-snps"));
        }
    }

   // TODO: test that spits out a useful exception when no MUMmer build compatible for User's machine
    @Test
    public void testNoMUMmerPackageForMachine(){

    }


}
