package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CompareReferencesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testCompareReferences() throws IOException{
        final File ref1 = new File("/Users/ocohen/workingcode/gatk/tempreferences/hg19mini.fasta");
        final File ref2 = new File("/Users/ocohen/workingcode/gatk/tempreferences/hg19mini_1renamed.fasta");
        //final File outputFile = createTempFile("", ".table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
               /* "-O", outputFile.getAbsolutePath()*/};
        runCommandLine(args);
    }
}
