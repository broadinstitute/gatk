package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.io.IOUtils;
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

    // no assertions made, testing tool runs successfully without -O argument
    @Test
    public void testCompareReferencesToStdOutput() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath()};
        runCommandLine(args);
    }

    // The following test runs the tool on different combinations of reference files
    // and produce output to stdout for the sake of manually inspecting outputs.
    // Disabled, as no actual assertions made.
    @Test(enabled = false)
    public void testCompareReferencesMultipleReferencesStdOut() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_1renamed.fasta");
        final File ref3 = new File(getToolTestDataDir() + "hg19mini_chr2snp.fasta");
        final File ref4 = new File(getToolTestDataDir() + "hg19mini_missingchr1.fasta");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref3.getAbsolutePath(), "-refcomp", ref3.getAbsolutePath(),
                "-refcomp", ref4.getAbsolutePath(), "-display-sequences-by-name", "-display-only-differing-sequences"};
        runCommandLine(args);
    }

    // FIND_SNPS_ONLY tests:
    @Test
    public void testFindSNPsMultipleSNPs() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2multiplesnps.fasta");
        final File output = IOUtils.createTempDir("tempFindSNPs");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FIND_SNPS_ONLY", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File actualOutput = new File(output, "hg19mini.fasta_hg19mini_chr2multiplesnps.fasta_snps.tsv");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_chr2multiplesnps.fasta_snps.tsv");
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test
    public void testFindSNPsIUPACBases() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2iupacsnps.fasta");
        final File output = IOUtils.createTempDir("tempFindSNPs");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FIND_SNPS_ONLY", "-base-comparison-output", output.toPath().toString()};
        runCommandLine(args);

        final File expectedOutput = new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_chr2iupacsnps.fasta_snps.tsv");
        final File actualOutput = new File(output, "hg19mini.fasta_hg19mini_chr2iupacsnps.fasta_snps.tsv");

        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    // FULL_ALIGNMENT tests:
    @Test
    public void testFullAlignmentModeMultipleSNPs() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr2multiplesnps.fasta");
        final File output = IOUtils.createTempDir("tempFullAlignmentSNPs");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FULL_ALIGNMENT", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File expectedOutput = new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_chr2multiplesnps.fasta.vcf");
        final File actualOutput = new File(output, "hg19mini.fasta_hg19mini_chr2multiplesnps.fasta.vcf");
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test
    public void testFullAlignmentModeDeletion() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_chr1indel.fasta");
        final File output = IOUtils.createTempDir("tempFullAlignmentIndel");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FULL_ALIGNMENT", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File expectedOutput = new File(getToolTestDataDir(), "expected.testDeletion.hg19mini.fasta_hg19mini_chr1indel.fasta.vcf");
        final File actualOutput = new File(output, "hg19mini.fasta_hg19mini_chr1indel.fasta.vcf");
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test
    public void testFullAlignmentModeInsertion() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini_chr1indel.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File output = IOUtils.createTempDir("tempFullAlignmentIndel");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FULL_ALIGNMENT", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File expectedOutput = new File(getToolTestDataDir(), "expected.testInsertion.hg19mini_chr1indel.fasta_hg19mini.fasta.vcf");
        final File actualOutput = new File(output, "hg19mini_chr1indel.fasta_hg19mini.fasta.vcf");
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test
    public void testFullAlignmentSNPsOnMultipleContigs() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_snpsmultiplecontigs.fasta");
        final File output = IOUtils.createTempDir("tempFullAlignmentSNPsMultipleContigs");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FULL_ALIGNMENT", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File expectedOutput = new File(getToolTestDataDir(), "expected.hg19mini.fasta_hg19mini_snpsmultiplecontigs.fasta.vcf");
        final File actualOutput = new File(output, "hg19mini.fasta_hg19mini_snpsmultiplecontigs.fasta.vcf");
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }

    @Test
    public void testFullAlignmentModeSNPAndINDEL() throws IOException{
        final File ref1 = new File(getToolTestDataDir() + "hg19mini.fasta");
        final File ref2 = new File(getToolTestDataDir() + "hg19mini_snpandindel.fasta");
        final File output = IOUtils.createTempDir("tempFullAlignmentIndel");

        final String[] args = new String[] {"-R", ref1.getAbsolutePath() , "-refcomp", ref2.getAbsolutePath(), "-base-comparison", "FULL_ALIGNMENT", "-base-comparison-output", output.getAbsolutePath()};
        runCommandLine(args);

        final File expectedOutput = new File(getToolTestDataDir(), "expected.SNPandINDEL.hg19mini.fasta_hg19mini_snpandindel.fasta.vcf");
        final File actualOutput = new File(output, "hg19mini.fasta_hg19mini_snpandindel.fasta.vcf");
        IntegrationTestSpec.assertEqualTextFiles(actualOutput, expectedOutput);
    }
}
