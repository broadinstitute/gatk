package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
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
    public void testCompareReferencesMultipleReferences() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");
        final File ref3 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesMultipleReferences", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesMultipleReferences.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-refcomp", ref3.getAbsolutePath(),
                "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test
    public void testCompareReferencesMissingValue() throws IOException{
        // SNP on chr2 causes the sequence to have a different MD5 than the hg19mini sequence of the same name
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File output = createTempFile("testCompareReferencesMissingValue", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testCompareReferencesMissingValue.table");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }



    @DataProvider(name = "md5ArgumentData")
    public Object[][] md5ArgumentData() {
        return new Object[][]{
                // ref1, ref2, expected output, md5 calculation mode
                new Object[]{ new File(getToolTestDataDir() + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesMissingValue.table"),
                        "USE_DICT"
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesRecalculateIfMissing.table"),
                        "RECALCULATE_IF_MISSING"
                },
                new Object[]{ new File(getToolTestDataDir() + "hg19mini_missingmd5.fasta"),
                        new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta"),
                        new File(getToolTestDataDir(), "expected.testCompareReferencesAlwaysRecalculateMD5.table"),
                        "ALWAYS_RECALCULATE"
                },
        };
    }

    @Test(dataProvider = "md5ArgumentData")
    public void testCompareReferencesMD5CalculationModes(File ref1, File ref2, File expected, String md5Mode) throws IOException {
        final File output = createTempFile("testCompareReferencesMD5CalculationModes", ".table");
        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(),
                "-O", output.getAbsolutePath(), "--md5-calculation-mode", md5Mode};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
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

    // The following three tests run the tool on different combinations of reference files
    // and produce output to stdout for the sake of manually inspecting outputs.
    // Disabled, as no actual assertions made.
    @Test/*(enabled = false)*/
    public void testCompareReferencesToStdOutput() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath()};
        runCommandLine(args);
    }

    @Test(enabled = false)
    public void testCompareReferencesMultipleReferencesStdOut() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");
        final File ref3 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File ref4 = new File(getToolTestDataDir() + "hg19mini_missingchr1.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-refcomp", ref3.getAbsolutePath(),
                "-refcomp", ref4.getAbsolutePath()};
        runCommandLine(args);
    }

    @Test(enabled = false)
    public void testCompareReferencesMissingSequencesStdOut() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_missingchr3.fasta");
        final File ref3 = new File(getToolTestDataDir() + "hg19mini_missingchr1.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-refcomp", ref3.getAbsolutePath()};
        runCommandLine(args);
    }

}
