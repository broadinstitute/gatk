package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;
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
    public void testPythonCommand(){
        String script = "/Users/ocohen/all2vcf/src/mummer";
        List<String> command = Arrays.asList("-h");
        MummerExecutor.runPythonCommand(script, command, null, null, true);
    }

    @Test
    public void testFindSNPS(){

    }


}