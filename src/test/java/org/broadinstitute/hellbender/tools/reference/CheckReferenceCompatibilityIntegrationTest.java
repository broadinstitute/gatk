package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

import static org.testng.Assert.*;

public class CheckReferenceCompatibilityIntegrationTest extends CommandLineProgramTest {

    private final String COMPARE_REFERENCES_FILES = "src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/";
    @Test
    public void testReferenceCompatibilityWithMD5s() throws IOException {
        final File ref1 = new File(COMPARE_REFERENCES_FILES + "hg19mini.fasta");
        final File dict = new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s.bam");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-I", dict.getAbsolutePath()};
        runCommandLine(args);
    }

    @Test
    public void testReferenceCompatibilityWithoutMD5s() throws IOException {
        final File ref1 = new File(COMPARE_REFERENCES_FILES + "hg19mini.fasta");
        final File dict = new File(getToolTestDataDir() + "reads_data_source_test1_withoutmd5s.bam");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-I", dict.getAbsolutePath()};
        runCommandLine(args);
    }

}