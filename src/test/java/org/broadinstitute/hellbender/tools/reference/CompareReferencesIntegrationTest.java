package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CompareReferencesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testCompareReferencesRenamedSequence() throws IOException{
        final File ref1 = new File("/Users/ocohen/workingcode/gatk/tempreferences/hg19mini.fasta");
        final File ref2 = new File("/Users/ocohen/workingcode/gatk/tempreferences/hg19mini_1renamed.fasta");
        //final File output = createTempFile("testCompareReferencesRenamedSequence", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesRenamedSequence.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", expectedOutput.getAbsolutePath()};
        runCommandLine(args);

        //IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesMissingValue() throws IOException{
        final File ref1 = new File("/Users/ocohen/workingcode/gatk/tempreferences/hg19mini.fasta");
        final File ref2 = new File("/Users/ocohen/workingcode/gatk/tempreferences/hg19mini_chr2snp.fasta");
        final File expectedOutput = new File("expected.testCompareReferencesMissingValue", ".tsv");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", expectedOutput.getAbsolutePath()};
        runCommandLine(args);
    }

}
