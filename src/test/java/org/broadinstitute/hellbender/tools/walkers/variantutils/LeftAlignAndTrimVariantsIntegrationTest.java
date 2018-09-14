package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.apache.log4j.WriterAppender;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.IOException;
import java.io.ByteArrayOutputStream;
import java.nio.file.*;
import java.util.*;

import org.apache.log4j.*;


public class LeftAlignAndTrimVariantsIntegrationTest extends CommandLineProgramTest {
    final Path testDataDir = Paths.get(getToolTestDataDir());

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
}
