package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;


public class LeftAlignAndTrimVariantsIntegrationTest extends CommandLineProgramTest {
    final Path testDataDir = Paths.get(getToolTestDataDir());

    // note: this test file has one particularly tricky case of a deletion AAA->A at chr21:13255301	that left aligns to CAA->A at chr21:13255289 by
    // jumping past a SNP A->G at chr21:13255296
    private Object[] getTestSet(String expectedVcf, String Options) {
        return new Object[]{testDataDir.resolve("test_left_align_hg38.vcf"), Paths.get(b38_reference_20_21), testDataDir.resolve(expectedVcf), Options};
    }

    @DataProvider(name = "LeftAlignDataProvider")
    public Object[][] LeftAlignTestData() {
        return new Object[][]{getTestSet("expected_left_align_hg38.vcf", ""),
                getTestSet("expected_left_align_hg38_split_multiallelics.vcf", " --" + LeftAlignAndTrimVariants.SPLIT_MULTIALLELEICS_LONG_NAME),
                getTestSet("expected_left_align_hg38_notrim.vcf", " --" + LeftAlignAndTrimVariants.DONT_TRIM_ALLELES_LONG_NAME),
                getTestSet("expected_left_align_hg38_notrim_split_multiallelics.vcf", " --" + LeftAlignAndTrimVariants.DONT_TRIM_ALLELES_SHORT_NAME + " --" + LeftAlignAndTrimVariants.SPLIT_MULTIALLELEICS_LONG_NAME),
                getTestSet("expected_left_align_hg38_split_multiallelics_keepOrigAC.vcf", " --" + LeftAlignAndTrimVariants.SPLIT_MULTIALLELEICS_LONG_NAME + " --" + LeftAlignAndTrimVariants.KEEP_ORIGINAL_AC_LONG_NAME),
                getTestSet("expected_left_align_hg38_maxIndelSize296.vcf", " --" + LeftAlignAndTrimVariants.MAX_INDEL_LENGTH_LONG_NAME + " 296"),
                getTestSet("expected_left_align_hg38_maxIndelSize342.vcf", " --" + LeftAlignAndTrimVariants.MAX_INDEL_LENGTH_LONG_NAME + " 342")
        };
    }

    @Test(dataProvider = "LeftAlignDataProvider")
    public void testLeftAlignment(Path inputFile, Path ref, Path expectedOutputFile, String options) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + ref.toString()
                        + " -V " + inputFile
                        + " -O %s"
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false"
                        + " --suppress-reference-path "
                        + options,
                Collections.singletonList(expectedOutputFile.toString())
        );
        spec.executeTest("testLeftAlignment--" + expectedOutputFile.toString(), this);
    }

    @DataProvider(name = "LeftAlignRequireReferenceDataProvider")
    public Object[][] LeftAlignRequireReferenceData() {
        return new Object[][]{{testDataDir.resolve("test_left_align_hg38.vcf")}};
    }

    @Test(dataProvider = "LeftAlignRequireReferenceDataProvider")
    public void testLefAlignRequireReference(Path inputFile) throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -V " + inputFile
                        + " -O %s"
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false"
                        + " --suppress-reference-path ",
                1, CommandLineException.MissingArgument.class
        );
        spec.executeTest("testLeftAlignment--requireReference", this);
    }

    @Test
    public void testSplitAllelesWithASFilters() throws IOException {
        Path inputFile = testDataDir.resolve("test_split_with_AS_filters.vcf");
        final File MITO_REF = new File(toolsTestDir, "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta");

        Path expectedOutputFile = testDataDir.resolve("expected_split_with_AS_filters.vcf");

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -R " + MITO_REF.getAbsolutePath()
                        + " -V " + inputFile
                        + " -O %s"
                        + " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false"
                        + " --suppress-reference-path "
                        + " --" + LeftAlignAndTrimVariants.SPLIT_MULTIALLELEICS_LONG_NAME
                        + " --" + LeftAlignAndTrimVariants.KEEP_ORIGINAL_AC_LONG_NAME,
                Collections.singletonList(expectedOutputFile.toString())
        );
        spec.executeTest("testSplitAllelesWithASFilters--" + expectedOutputFile.toString(), this);

    }
}
