package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.runtime.ProcessController;
import org.broadinstitute.hellbender.utils.runtime.ProcessOutput;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

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

    @Test
    public void testPythonCommand() throws InterruptedException{
        //MummerExecutor exec = new MummerExecutor(new File("/Users/ocohen/all2vcf/src/mummer"));
        String script = "/Users/ocohen/all2vcf/src/mummer";
        List<String> command = Arrays.asList("-h");
        MummerExecutor.runPythonCommand(script, command, null, null, true);
    }

    @Test
    public void testExecuteMummer() {
        /*File fasta1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        File fasta2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        MummerExecutor exec = new MummerExecutor(new File("/Users/ocohen/workingcode/MUMmer3.23/"));

        File deltaFile = new File(getToolTestDataDir() + "nucmerOutput");
        String[] nucmerArgs = { exec.getMummerExecutableDirectory() +"/nucmer", "--mum", "-p", deltaFile.getAbsolutePath(), fasta1.getAbsolutePath(), fasta2.getAbsolutePath()};
        ProcessOutput nucmer = MummerExecutor.runShellCommand(nucmerArgs, null, null,false);

        File deltaFilterFile = new File(getToolTestDataDir() + "deltaFilterOutput.delta"); // file for delta filter output --> input to show-snps
        String[] deltaFilterArgs = {exec.getMummerExecutableDirectory() + "/delta-filter", "-1", deltaFile.getAbsolutePath() + ".delta"*//*, ">", deltaFilterFile.getAbsolutePath()*//*};
        ProcessOutput deltaFilter = MummerExecutor.runShellCommand(deltaFilterArgs, null, deltaFilterFile, false);

        // show snps
        File showSNPSOutput = new File(getToolTestDataDir() + "showSNPsOutput.snps");
        String[] showSNPsArgs = {exec.getMummerExecutableDirectory() + "/show-snps", "-rlTH", deltaFilterFile.getAbsolutePath()};
        ProcessOutput showSNPs = MummerExecutor.runShellCommand(showSNPsArgs, null, showSNPSOutput, false);

        File vcf = new File(getToolTestDataDir() + "all2vcfOutput.vcf");
        String script = exec.getMummerExecutableDirectory() + "/all2vcf";
        List<String> command = Arrays.asList("--snps", showSNPSOutput.getAbsolutePath(), "--reference", fasta1.getAbsolutePath(), "--output-header");
        ProcessOutput all2vcf = MummerExecutor.runPythonCommand(script, command, null, vcf, false);
        */

        File fasta1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        File fasta2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        File outputDirectory = new File("/Users/ocohen/workingcode/gatk/tempreferences/testingMummerExecution/");
        MummerExecutor exec = new MummerExecutor(new File("/Users/ocohen/workingcode/MUMmer3.23/"));

        exec.executeMummer(fasta1, fasta2, outputDirectory);


    }

}