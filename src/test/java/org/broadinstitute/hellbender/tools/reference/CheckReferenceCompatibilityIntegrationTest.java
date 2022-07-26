package org.broadinstitute.hellbender.tools.reference;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class CheckReferenceCompatibilityIntegrationTest extends CommandLineProgramTest {

    private final String COMPARE_REFERENCES_TEST_FILES_DIRECTORY = "src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences/";

    @DataProvider(name = "testReferenceCompatibilityBAMWithMD5sData")
    public Object[][] testReferenceCompatibilityBAMWithMD5sData() {
        return new Object[][]{
                // reference, dictionary for comparison, expected output
                new Object[]{
                        // exact match
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s.bam"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityBAMWithMD5s_exactmatch.table")
                },
                new Object[]{
                        // subset
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s_missingchr1.bam"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityBAMWithMD5s_subset.table")
                },
                new Object[]{
                        // not compatible
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_missingchr3.fasta"),
                        new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s.bam"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityBAMWithMD5s_notcompatible.table")
                },
        };
    }

    @Test(dataProvider = "testReferenceCompatibilityBAMWithMD5sData")
    public void testReferenceCompatibilityBAMWithMD5s(File ref1, File dict, File expectedOutput) throws IOException {
        final File output = createTempFile("testReferenceCompatibilityBAMsWithMD5s", ".table");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-I", dict.getAbsolutePath(), "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @DataProvider(name = "testReferenceCompatibilityBAMWithoutMD5sData")
    public Object[][] testReferenceCompatibilityBAMWithoutMD5sData() {
        return new Object[][]{
                // reference, dictionary for comparison, expected output
                new Object[]{
                        // exact match
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "reads_data_source_test1_withoutmd5s.bam"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityBAMWithoutMD5s_compatible.table")
                },
                new Object[]{
                        // subset
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "reads_data_source_test1_withoutmd5s_missingchr1.bam"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityBAMWithoutMD5s_subset.table")
                },
                new Object[]{
                        // not compatible
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_missingchr3.fasta"),
                        new File(getToolTestDataDir() + "reads_data_source_test1_withoutmd5s.bam"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityBAMWithoutMD5s_notcompatible.table")
                },
        };
    }

    @Test(dataProvider = "testReferenceCompatibilityBAMWithoutMD5sData")
    public void testReferenceCompatibilityBAMWithoutMD5s(File ref1, File dict, File expectedOutput) throws IOException {
        final File output = createTempFile("testReferenceCompatibilityBAMsWithoutMD5s", ".table");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-I", dict.getAbsolutePath(), "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReferenceCompatibilityMultipleBAMs() throws IOException {
        final File ref1 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta");
        final File dict1 = new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s.bam");
        final File dict2 = new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s_missingchr1.bam");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-I", dict1.getAbsolutePath(), "-I", dict2.getAbsolutePath()};
        runCommandLine(args);
    }

    @Test
    public void testReferenceCompatibilityMultipleReferencesBAMWithMD5s() throws IOException {
        final File ref1 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta");
        final File ref2 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_1renamed.fasta");
        final File ref3 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_chr2snp.fasta");
        final File dict = new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s_missingchr1.bam");
        final File output = createTempFile("testReferenceCompatibilityMultipleReferencesWithMD5s", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testReferenceCompatibilityMultipleReferencesBAMWithMD5s.table");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath(), "-refcomp", ref2.getAbsolutePath(),
                "-refcomp", ref3.getAbsolutePath(), "-I", dict.getAbsolutePath(), "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    // Note that hg19mini_chr2snp.fasta is NOT a valid reference, but tool renders them a match due to compatibility by name and length.
    // Included to demonstrate that this is the best the tool can do without MD5s present
    @Test
    public void testReferenceCompatibilityMultipleReferencesBAMWithoutMD5s() throws IOException {
        final File ref1 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta");
        final File ref2 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_1renamed.fasta");
        final File ref3 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_chr2snp.fasta");
        final File dict = new File(getToolTestDataDir() + "reads_data_source_test1_withoutmd5s.bam");
        final File output = createTempFile("testReferenceCompatibilityMultipleReferencesWithoutMD5s", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testReferenceCompatibilityMultipleReferencesBAMWithoutMD5s.table");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath(), "-refcomp", ref2.getAbsolutePath(),
                "-refcomp", ref3.getAbsolutePath(), "-I", dict.getAbsolutePath(), "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    // References are not compatible since one has SNP on chr 2, but tool renders them seemingly a match
    // Included to demonstrate that this is the best the tool can do without MD5s present
    @Test
    public void testReferenceCompatibilityWithoutMD5sMismatch() throws IOException {
        final File ref1 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_chr2snp.fasta");
        final File dict = new File(getToolTestDataDir() + "reads_data_source_test1_withoutmd5s.bam");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-I", dict.getAbsolutePath()};
        runCommandLine(args);
    }

    // TODO: compatibility based on MD5 faulty since MD5s not in sequence dictionary (see ticket #730 "VCFHeader drops sequence dictionary attributes")
    @Test(enabled = false)
    public void testReferenceCompatibilityVCFWithMD5s() throws IOException {
        final File ref1 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_missingchr1.fasta");
        final File dict = new File(getToolTestDataDir() + "example_variants_withSequenceDict.vcf");
        //final File output = createTempFile("testReferenceCompatibilityVCFWithMD5s", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testReferenceCompatibilityVCFWithMD5s.table");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-V", dict.getAbsolutePath()/*, "-O", expectedOutput.getAbsolutePath()*/};
        runCommandLine(args);

        //IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    @DataProvider(name = "testReferenceCompatibilityVCFWithoutMD5sData")
    public Object[][] testReferenceCompatibilityVCFWithoutMD5sData() {
        return new Object[][]{
                // reference, dictionary for comparison, expected output
                new Object[]{
                        // exact match
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "example_variants_withSequenceDict_withoutmd5s.vcf"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityVCFWithoutMD5s_match.table")
                },
                new Object[]{
                        // subset
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta"),
                        new File(getToolTestDataDir() + "example_variants_withSequenceDict_withoutmd5s_missingchr1.vcf"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityVCFWithoutMD5s_subset.table")
                },
                new Object[]{
                        // not compatible
                        new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_missingchr3.fasta"),
                        new File(getToolTestDataDir() + "example_variants_withSequenceDict_withoutmd5s.vcf"),
                        new File(getToolTestDataDir(), "expected.testReferenceCompatibilityVCFWithoutMD5s_notcompatible.table")
                },
        };
    }

    @Test(dataProvider = "testReferenceCompatibilityVCFWithoutMD5sData")
    public void testReferenceCompatibilityVCFWithoutMD5s(File ref1, File dict, File expectedOutput) throws IOException {
        final File output = createTempFile("testReferenceCompatibilityVCFWithoutMD5s", ".table");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath() , "-V", dict.getAbsolutePath(), "-O", output.getAbsolutePath()};
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expectedOutput);
    }

    // for quick stdout testing
    @Test(enabled = false)
    public void testStdOutput() throws IOException{
        final File ref1 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini.fasta");
        final File ref2 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_1renamed.fasta");
        final File ref3 = new File(COMPARE_REFERENCES_TEST_FILES_DIRECTORY + "hg19mini_chr2snp.fasta");
        final File dict = new File(getToolTestDataDir() + "reads_data_source_test1_withmd5s_missingchr1.bam");
        //final File output = createTempFile("testReferenceCompatibilityMultipleReferencesWithMD5s", ".table");
        final File expectedOutput = new File(getToolTestDataDir(), "expected.testReferenceCompatibilityMultipleReferencesBAMWithMD5s.table");

        final String[] args = new String[] {"-refcomp", ref1.getAbsolutePath(), "-refcomp", ref2.getAbsolutePath(),
                "-refcomp", ref3.getAbsolutePath(), "-I", dict.getAbsolutePath(), "-O", expectedOutput.getAbsolutePath()};
        runCommandLine(args);
    }

}