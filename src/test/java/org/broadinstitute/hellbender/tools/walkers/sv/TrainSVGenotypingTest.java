package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.GatkToolIntegrationTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class TrainSVGenotypingTest extends GatkToolIntegrationTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TOOL_TEST_DIR = toolsTestDir + "walkers/sv/TrainSVGenotyping/";
    public static final String LARGE_FILE_TEST_DIR = largeFileTestDir + "sv/TrainSVGenotyping/";
    public static final String AGGREGATE_SV_TEST_DIR = toolsTestDir + "walkers/sv/AggregateSVEvidence/";

    public static final double FLOAT_TOLERANCE = 1e-4;

    public static final List<String> FORMAT_ATTRIBUTES_WITH_JITTER = Arrays.asList(
            GATKSVVCFConstants.DEPTH_MEDIAN_COPY_RATIO
    );

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS,
                "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testTrainSVGenotyping() throws IOException {
        final File vcfFile = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.vcf.gz");
        final File peFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.pe.txt.gz");
        final File srFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.sr.txt.gz");
        final File rdFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.rd.txt.gz");
        final File trainingIntervals = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.training_intervals.bed");
        final File depthExclude = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.depth_exclude.bed.gz");
        final File pesrExclude = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.pesr_exclude.bed.gz");

        final File ploidyTable = new File(AGGREGATE_SV_TEST_DIR, "1kg_ref_panel_v1.ploidy_table.tsv");
        final File coverageFile = new File(AGGREGATE_SV_TEST_DIR, "aggregate_sv_test.medianCov.tsv");

        final String outputName = "train_sv_test";
        final Path outputDir;
        final String outputVcfPath;
        if (UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            outputDir = Path.of(TOOL_TEST_DIR);
            outputVcfPath = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.vcf.gz").getAbsolutePath();
        } else {
            outputDir = createTempDir("trainSVGenotyping").toPath();
            outputVcfPath = createTempFile("train_sv_genotyping_output", ".vcf.gz").getAbsolutePath();
        }

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(vcfFile)
                .addOutput(outputVcfPath)
                .addReference(hg38Reference)
                .add(TrainSVGenotyping.DEPTH_EVIDENCE_FILE_PATH_LONG_NAME, rdFile)
                .add(AggregateSVEvidence.DISCORDANT_PAIRS_LONG_NAME, peFile)
                .add(AggregateSVEvidence.SPLIT_READ_LONG_NAME, srFile)
                .add(TrainSVGenotyping.MEDIAN_COUNTS_FILE_PATH_LONG_NAME, coverageFile)
                .add(SVClusterWalker.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(TrainSVGenotyping.TRAINING_INTERVALS_LONG_NAME, trainingIntervals)
                .add(TrainSVGenotyping.DEPTH_EXCLUSION_INTERVALS_LONG_NAME, depthExclude)
                .add(TrainSVGenotyping.PESR_EXCLUSION_INTERVALS_LONG_NAME, pesrExclude)
                .add(TrainSVGenotyping.TABLES_DIR_LONG_NAME, outputDir.toString())
                .add(TrainSVGenotyping.TABLES_NAME_LONG_NAME, outputName)
                .add(TrainSVGenotyping.OUTPUT_TRAINING_VCF_LONG_NAME, true);

        runCommandLine(args, TrainSVGenotyping.class.getSimpleName());

        final File rdDepthParamsFile = outputDir.resolve(outputName + ".rd_depth_geno_params.tsv").toFile();
        final File rdPesrParamsFile = outputDir.resolve(outputName + ".rd_pesr_geno_params.tsv").toFile();
        final File peParamsFile = outputDir.resolve(outputName + ".pe_geno_params.tsv").toFile();
        final File srParamsFile = outputDir.resolve(outputName + ".sr_geno_params.tsv").toFile();
        Assert.assertTrue(rdDepthParamsFile.exists(), "RD depth genotype parameters file should exist");
        Assert.assertTrue(rdPesrParamsFile.exists(), "RD PESR genotype parameters file should exist");
        Assert.assertTrue(peParamsFile.exists(), "PE genotype parameters file should exist");
        Assert.assertTrue(srParamsFile.exists(), "SR genotype parameters file should exist");

        if (UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            Files.copy(rdDepthParamsFile.toPath(), Path.of(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.rd_geno_params.tsv"), StandardCopyOption.REPLACE_EXISTING);
            Files.copy(peParamsFile.toPath(), Path.of(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.pe_geno_params.tsv"), StandardCopyOption.REPLACE_EXISTING);
            Files.copy(srParamsFile.toPath(), Path.of(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.sr_geno_params.tsv"), StandardCopyOption.REPLACE_EXISTING);
            return;
        }

        final File expectedVcfFile = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.vcf.gz");
        final Pair<VCFHeader, List<VariantContext>> expected = VariantContextTestUtils.readEntireVCFIntoMemory(expectedVcfFile.getPath());
        final Pair<VCFHeader, List<VariantContext>> output = VariantContextTestUtils.readEntireVCFIntoMemory(outputVcfPath);

        Assert.assertEquals(output.getValue().size(), expected.getValue().size(),
                "Output and expected VCFs should have the same number of variants");
        final Iterator<VariantContext> expectedIter = expected.getValue().iterator();
        final Iterator<VariantContext> outputIter = output.getValue().iterator();
        while (expectedIter.hasNext()) {
            final VariantContext expectedVariant = expectedIter.next();
            final VariantContext actualVariant = outputIter.next();
            VariantContextTestUtils.assertVariantContextsAreEqual(actualVariant, expectedVariant,
                    Collections.emptyList(), FORMAT_ATTRIBUTES_WITH_JITTER);

            for (final Genotype expectedGenotype : expectedVariant.getGenotypes()) {
                final Genotype actualGenotype = actualVariant.getGenotype(expectedGenotype.getSampleName());
                for (final String attribute : FORMAT_ATTRIBUTES_WITH_JITTER) {
                    final Object expectedObj = expectedGenotype.getExtendedAttribute(attribute);
                    final Object actualObj = actualGenotype.getExtendedAttribute(attribute);
                    if (expectedObj == null && actualObj == null) {
                        continue;
                    }
                    Assert.assertNotNull(expectedObj, "Expected attribute " + attribute + " should not be null");
                    Assert.assertNotNull(actualObj, "Attribute " + attribute + " should not be null when expected is non-null");
                    final double expectedValue = Double.parseDouble(String.valueOf(expectedObj));
                    final double actualValue = Double.parseDouble(String.valueOf(actualObj));
                    if (Double.isNaN(expectedValue)) {
                        Assert.assertTrue(Double.isNaN(actualValue),
                                "Attribute " + attribute + " should be NaN for sample " + expectedGenotype.getSampleName());
                    } else {
                        Assert.assertTrue(Math.abs(expectedValue - actualValue) < FLOAT_TOLERANCE,
                                "Attribute " + attribute + " for sample " + expectedGenotype.getSampleName()
                                        + " expected " + expectedValue + " but got " + actualValue);
                    }
                }
            }
        }

        assertTsvFilesEqual(rdDepthParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.rd_geno_params.tsv"));
        assertTsvFilesEqual(rdPesrParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.rd_geno_params.tsv"));
        assertTsvFilesEqual(peParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.pe_geno_params.tsv"));
        assertTsvFilesEqual(srParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.sr_geno_params.tsv"));
    }

    @Test
    public void testTrainSVGenotypingTablesOnly() throws IOException {
        // Run without --output-training-vcf (defaults to false) to verify that
        // parameter tables are produced correctly via the SR-only histogram path.
        final File vcfFile = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.vcf.gz");
        final File peFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.pe.txt.gz");
        final File srFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.sr.txt.gz");
        final File rdFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.rd.txt.gz");
        final File trainingIntervals = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.training_intervals.bed");
        final File depthExclude = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.depth_exclude.bed.gz");
        final File pesrExclude = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.pesr_exclude.bed.gz");

        final File ploidyTable = new File(AGGREGATE_SV_TEST_DIR, "1kg_ref_panel_v1.ploidy_table.tsv");
        final File coverageFile = new File(AGGREGATE_SV_TEST_DIR, "aggregate_sv_test.medianCov.tsv");

        final String outputName = "train_sv_test_tables_only";
        final Path outputDir = createTempDir("trainSVGenotypingTablesOnly").toPath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(vcfFile)
                .addReference(hg38Reference)
                .add(TrainSVGenotyping.DEPTH_EVIDENCE_FILE_PATH_LONG_NAME, rdFile)
                .add(AggregateSVEvidence.DISCORDANT_PAIRS_LONG_NAME, peFile)
                .add(AggregateSVEvidence.SPLIT_READ_LONG_NAME, srFile)
                .add(TrainSVGenotyping.MEDIAN_COUNTS_FILE_PATH_LONG_NAME, coverageFile)
                .add(SVClusterWalker.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(TrainSVGenotyping.TRAINING_INTERVALS_LONG_NAME, trainingIntervals)
                .add(TrainSVGenotyping.DEPTH_EXCLUSION_INTERVALS_LONG_NAME, depthExclude)
                .add(TrainSVGenotyping.PESR_EXCLUSION_INTERVALS_LONG_NAME, pesrExclude)
                .add(TrainSVGenotyping.TABLES_DIR_LONG_NAME, outputDir.toString())
                .add(TrainSVGenotyping.TABLES_NAME_LONG_NAME, outputName);

        runCommandLine(args, TrainSVGenotyping.class.getSimpleName());

        final File rdDepthParamsFile = outputDir.resolve(outputName + ".rd_depth_geno_params.tsv").toFile();
        final File rdPesrParamsFile = outputDir.resolve(outputName + ".rd_pesr_geno_params.tsv").toFile();
        final File peParamsFile = outputDir.resolve(outputName + ".pe_geno_params.tsv").toFile();
        final File srParamsFile = outputDir.resolve(outputName + ".sr_geno_params.tsv").toFile();
        Assert.assertTrue(rdDepthParamsFile.exists(), "RD depth genotype parameters file should exist");
        Assert.assertTrue(rdPesrParamsFile.exists(), "RD PESR genotype parameters file should exist");
        Assert.assertTrue(peParamsFile.exists(), "PE genotype parameters file should exist");
        Assert.assertTrue(srParamsFile.exists(), "SR genotype parameters file should exist");

        assertTsvFilesEqual(rdDepthParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.rd_geno_params.tsv"));
        assertTsvFilesEqual(rdPesrParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.rd_geno_params.tsv"));
        assertTsvFilesEqual(peParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.pe_geno_params.tsv"));
        assertTsvFilesEqual(srParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.sr_geno_params.tsv"));
    }

    private static void assertTsvFilesEqual(final File actual, final File expected) throws IOException {
        final List<String> actualLines = Files.readAllLines(actual.toPath());
        final List<String> expectedLines = Files.readAllLines(expected.toPath());
        Assert.assertEquals(actualLines.size(), expectedLines.size(),
                "TSV files should have the same number of lines: " + actual.getName() + " vs " + expected.getName());
        for (int i = 0; i < expectedLines.size(); i++) {
            final String[] actualTokens = actualLines.get(i).split("\t");
            final String[] expectedTokens = expectedLines.get(i).split("\t");
            Assert.assertEquals(actualTokens.length, expectedTokens.length,
                    "Line " + (i + 1) + " should have the same number of columns");
            for (int j = 0; j < expectedTokens.length; j++) {
                try {
                    final double expectedVal = Double.parseDouble(expectedTokens[j]);
                    final double actualVal = Double.parseDouble(actualTokens[j]);
                    if (Double.isNaN(expectedVal)) {
                        Assert.assertTrue(Double.isNaN(actualVal),
                                "Expected NaN at line " + (i + 1) + " column " + (j + 1));
                    } else if (expectedVal != 0) {
                        Assert.assertTrue(Math.abs(expectedVal - actualVal) / Math.abs(expectedVal) < FLOAT_TOLERANCE,
                                "Line " + (i + 1) + " column " + (j + 1) + " expected " + expectedVal + " but got " + actualVal);
                    } else {
                        Assert.assertTrue(Math.abs(actualVal) < FLOAT_TOLERANCE,
                                "Line " + (i + 1) + " column " + (j + 1) + " expected 0 but got " + actualVal);
                    }
                } catch (final NumberFormatException e) {
                    Assert.assertEquals(actualTokens[j], expectedTokens[j],
                            "Line " + (i + 1) + " column " + (j + 1) + " string mismatch");
                }
            }
        }
    }

}
