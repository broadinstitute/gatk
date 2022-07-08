package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CompareReferencesIntegrationTest extends CommandLineProgramTest {

    @Test
    public void testCompareReferencesRenamedSequence() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");
        final File output = createTempFile("testCompareReferencesRenamedSequence", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesRenamedSequence.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesMissingValue() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesMissingValue", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesMissingValue.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesMissingValueUseDictMD5() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesUseDictMD5", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesMissingValue.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath(), "--md5-calculation-mode", "USE_DICT"};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCompareReferencesUseDictMD5MissingValue() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesUseDictMD5", ".table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath(), "--md5-calculation-mode", "USE_DICT"};
        runCommandLine(args);
    }

    @Test
    public void testCompareReferencesMissingValueRecalculateIfMissingMD5() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesRecalculateIfMissingMD5", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesRecalculateIfMissing.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath(), "--md5-calculation-mode", "RECALCULATE_IF_MISSING"};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesMissingValueAlwaysRecalculateMD5() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesAlwaysRecalculateMD5", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesAlwaysRecalculateMD5.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath(), "--md5-calculation-mode", "ALWAYS_RECALCULATE"};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesToStdOutput() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath()};
        runCommandLine(args);
    }

}
