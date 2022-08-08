package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.alignment.MummerExecutor;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CompareReferencesUnitTest extends CommandLineProgramTest {

    @Test
    public void testGenerateFastaForSequence() throws IOException {
        File ref = new File(getToolTestDataDir() + "hg19mini.fasta");
        File expectedOutput = new File("/Users/ocohen/workingcode/gatk/tempreferences/1.fasta");
        String sequenceName = "1";
        File output = createTempFile("example_chr1", ".fasta");

        ReferenceDataSource source = ReferenceDataSource.of(ref.toPath());
        CompareReferences.generateFastaForSequence(source, sequenceName, new GATKPath(output.toString()));

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testRunShellCommand(){
        String[] command = {"echo", "hello"};

        MummerExecutor.runShellCommand(command, null, new File("/Users/ocohen/workingcode/hello.output"), true);
    }

  /*  @Test
    public void testExecuteMummer() {
        File fasta1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        File fasta2 = new File(getToolTestDataDir() + "hg19mini_chr2multiplesnps.fasta");
        File outputDirectory = new File("/Users/ocohen/workingcode/gatk/tempreferences/");
        MummerExecutor exec = new MummerExecutor();

        exec.executeMummer(fasta1, fasta2, outputDirectory);
        //File expectedOutput = new File(outputDirectory, "snps_output.snps");
    }*/

    @Test
    public void testPrepareMUMmerExecutionDirectory(){
        MummerExecutor exec = new MummerExecutor();
        File executableDirectory = exec.getMummerExecutableDirectory();
        Assert.assertEquals(executableDirectory.listFiles().length, 4);

        for(File file : executableDirectory.listFiles()){
            Assert.assertTrue(file.getTotalSpace() > 0);
            Assert.assertTrue(file.canExecute());
        }
    }
}