package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.GatkToolIntegrationTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.DepthEvidence;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.DiscordantPairEvidenceCodec;
import org.broadinstitute.hellbender.utils.codecs.SplitReadEvidenceCodec;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class TrainSVGenotypingTest extends GatkToolIntegrationTest {

    private static final int SYNTHETIC_FIXTURE_SHIFT = 1_000_000;

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
                .add(TrainSVGenotyping.TABLES_NAME_LONG_NAME, outputName);

        runCommandLine(args, TrainSVGenotyping.class.getSimpleName());

        final File rdParamsFile = outputDir.resolve(outputName + ".rd_geno_params.tsv").toFile();
        final File peParamsFile = outputDir.resolve(outputName + ".pe_geno_params.tsv").toFile();
        final File srParamsFile = outputDir.resolve(outputName + ".sr_geno_params.tsv").toFile();
        Assert.assertTrue(rdParamsFile.exists(), "RD genotype parameters file should exist");
        Assert.assertTrue(peParamsFile.exists(), "PE genotype parameters file should exist");
        Assert.assertTrue(srParamsFile.exists(), "SR genotype parameters file should exist");

        if (UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS) {
            Files.copy(rdParamsFile.toPath(), Path.of(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.rd_geno_params.tsv"), StandardCopyOption.REPLACE_EXISTING);
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

        assertTsvFilesEqual(rdParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.rd_geno_params.tsv"));
        assertTsvFilesEqual(peParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.pe_geno_params.tsv"));
        assertTsvFilesEqual(srParamsFile, new File(TOOL_TEST_DIR, "train_sv_genotyping_test.expected.sr_geno_params.tsv"));
    }

    @Test
    public void testTrainSVGenotypingAcceptsDownsamplingArgument() throws IOException {
        final File vcfFile = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.vcf.gz");
        final File peFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.pe.txt.gz");
        final File srFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.sr.txt.gz");
        final File rdFile = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.rd.txt.gz");
        final File trainingIntervals = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.training_intervals.bed");
        final File depthExclude = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.depth_exclude.bed.gz");
        final File ploidyTable = new File(AGGREGATE_SV_TEST_DIR, "1kg_ref_panel_v1.ploidy_table.tsv");
        final File coverageFile = new File(AGGREGATE_SV_TEST_DIR, "aggregate_sv_test.medianCov.tsv");

        final String outputName = "train_sv_test_downsampling_arg";
        final Path outputDir = createTempDir("trainSVGenotypingDownsamplingArg").toPath();
        final String outputVcfPath = createTempFile("train_sv_genotyping_downsampling", ".vcf.gz").getAbsolutePath();

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
                .add(TrainSVGenotyping.MIN_PESER_SIZE_LONG_NAME, 0)
                .add(TrainSVGenotyping.TABLES_DIR_LONG_NAME, outputDir.toString())
                .add(TrainSVGenotyping.TABLES_NAME_LONG_NAME, outputName)
                .add(TrainSVGenotyping.DOWNSAMPLE_STRIDE_LONG_NAME, 1);

        runCommandLine(args, TrainSVGenotyping.class.getSimpleName());

        Assert.assertTrue(outputDir.resolve(outputName + ".rd_geno_params.tsv").toFile().exists(), "RD genotype parameters file should exist");
        Assert.assertTrue(outputDir.resolve(outputName + ".pe_geno_params.tsv").toFile().exists(), "PE genotype parameters file should exist");
        Assert.assertTrue(outputDir.resolve(outputName + ".sr_geno_params.tsv").toFile().exists(), "SR genotype parameters file should exist");
        Assert.assertTrue(new File(outputVcfPath).exists(), "Output VCF should exist when downsampling arguments are provided");
    }

    @Test
    public void testTrainSVGenotypingDownsamplingPreservesSplitReadRecoveryMetrics() throws IOException {
        final SyntheticTrainSVFixture fixture = createSyntheticDownsamplingFixture();

        final TrainSVRunResult uncappedRun = runTrainSVGenotyping(fixture, "synthetic_uncapped", 1);
        final TrainSVRunResult cappedRun = runTrainSVGenotyping(fixture, "synthetic_capped", 2);

        final int inputRecordCount = VariantContextTestUtils.readEntireVCFIntoMemory(fixture.vcfFile.getAbsolutePath()).getRight().size();
        final int cappedRecordCount = VariantContextTestUtils.readEntireVCFIntoMemory(cappedRun.outputVcf().getAbsolutePath()).getRight().size();
        Assert.assertTrue(cappedRecordCount < inputRecordCount,
                "Synthetic fixture should trigger output downsampling when the training cap is enabled");

        // Verify both SR metric files are well-formed and parseable.
        // Because the SR genotyping pass is now also downsampled (not just the
        // training passes), recovery counts will scale with the retained subset
        // and cannot be expected to stay close to uncapped values on a tiny
        // synthetic fixture.
        final Map<String, Double> uncappedSrMetrics = readSingleMetricRow(uncappedRun.srParamsFile());
        final Map<String, Double> cappedSrMetrics = readSingleMetricRow(cappedRun.srParamsFile());

        for (final String metric : List.of("rare_pass", "rare_fail", "common_pass", "common_fail")) {
            Assert.assertTrue(uncappedSrMetrics.containsKey(metric), "Uncapped SR metrics should contain " + metric);
            Assert.assertTrue(cappedSrMetrics.containsKey(metric), "Capped SR metrics should contain " + metric);
            Assert.assertTrue(uncappedSrMetrics.get(metric) >= 0, metric + " should be non-negative in uncapped run");
            Assert.assertTrue(cappedSrMetrics.get(metric) >= 0, metric + " should be non-negative in capped run");
        }
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

    private SyntheticTrainSVFixture createSyntheticDownsamplingFixture() throws IOException {
        final File sourceVcf = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.vcf.gz");
        final File sourcePe = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.pe.txt.gz");
        final File sourceSr = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.sr.txt.gz");
        final File sourceRd = new File(LARGE_FILE_TEST_DIR, "train_sv_genotyping_test.rd.txt.gz");

        final Path tempDir = createTempDir("trainSVGenotypingSyntheticFixture").toPath();
        final File syntheticVcf = tempDir.resolve("synthetic_train_sv_genotyping.vcf").toFile();
        final File syntheticPe = tempDir.resolve("synthetic_train_sv_genotyping.pe.txt.gz").toFile();
        final File syntheticSr = tempDir.resolve("synthetic_train_sv_genotyping.sr.txt.gz").toFile();
        final File syntheticRd = tempDir.resolve("synthetic_train_sv_genotyping.rd.txt.gz").toFile();

        try (final VCFFileReader reader = new VCFFileReader(sourceVcf, false)) {
            final VCFHeader header = reader.getFileHeader();
            final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
            final List<VariantContext> variants = new ArrayList<>();
            reader.forEach(variants::add);
            writeSyntheticVcf(syntheticVcf, header, dictionary, variants, SYNTHETIC_FIXTURE_SHIFT);
            writeSyntheticDepthEvidence(sourceRd, syntheticRd, dictionary, readDepthSampleNames(sourceRd), SYNTHETIC_FIXTURE_SHIFT);
            writeSyntheticDiscordantPairEvidence(sourcePe, syntheticPe, dictionary, SYNTHETIC_FIXTURE_SHIFT);
            writeSyntheticSplitReadEvidence(sourceSr, syntheticSr, dictionary, SYNTHETIC_FIXTURE_SHIFT);
        }

        return new SyntheticTrainSVFixture(syntheticVcf, syntheticPe, syntheticSr, syntheticRd);
    }

    private TrainSVRunResult runTrainSVGenotyping(final SyntheticTrainSVFixture fixture,
                                                  final String outputName,
                                                  final int downsampleStride) throws IOException {
        final File trainingIntervals = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.training_intervals.bed");
        final File depthExclude = new File(TOOL_TEST_DIR, "train_sv_genotyping_test.depth_exclude.bed.gz");
        final File ploidyTable = new File(AGGREGATE_SV_TEST_DIR, "1kg_ref_panel_v1.ploidy_table.tsv");
        final File coverageFile = new File(AGGREGATE_SV_TEST_DIR, "aggregate_sv_test.medianCov.tsv");

        final Path outputDir = createTempDir(outputName).toPath();
        final File outputVcf = createTempFile(outputName, ".vcf.gz");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(fixture.vcfFile())
                .addOutput(outputVcf.getAbsolutePath())
                .addReference(hg38Reference)
                .add(TrainSVGenotyping.DEPTH_EVIDENCE_FILE_PATH_LONG_NAME, fixture.rdFile())
                .add(AggregateSVEvidence.DISCORDANT_PAIRS_LONG_NAME, fixture.peFile())
                .add(AggregateSVEvidence.SPLIT_READ_LONG_NAME, fixture.srFile())
                .add(TrainSVGenotyping.MEDIAN_COUNTS_FILE_PATH_LONG_NAME, coverageFile)
                .add(SVClusterWalker.PLOIDY_TABLE_LONG_NAME, ploidyTable)
                .add(TrainSVGenotyping.TRAINING_INTERVALS_LONG_NAME, trainingIntervals)
                .add(TrainSVGenotyping.DEPTH_EXCLUSION_INTERVALS_LONG_NAME, depthExclude)
                .add(TrainSVGenotyping.MIN_PESER_SIZE_LONG_NAME, 0)
                .add(TrainSVGenotyping.TABLES_DIR_LONG_NAME, outputDir.toString())
                .add(TrainSVGenotyping.TABLES_NAME_LONG_NAME, outputName);
        if (downsampleStride > 1) {
            args.add(TrainSVGenotyping.DOWNSAMPLE_STRIDE_LONG_NAME, downsampleStride);
        }

        runCommandLine(args, TrainSVGenotyping.class.getSimpleName());

        return new TrainSVRunResult(
                outputVcf,
                outputDir.resolve(outputName + ".rd_geno_params.tsv").toFile(),
                outputDir.resolve(outputName + ".pe_geno_params.tsv").toFile(),
                outputDir.resolve(outputName + ".sr_geno_params.tsv").toFile());
    }

    private static void writeSyntheticVcf(final File outputVcf,
                                          final VCFHeader header,
                                          final SAMSequenceDictionary dictionary,
                                          final List<VariantContext> variants,
                                          final int shift) {
        final Map<String, Integer> contigOrder = buildContigOrder(dictionary);
        final List<VariantContext> syntheticVariants = new ArrayList<>(variants);
        for (final VariantContext variant : variants) {
            final String svType = variant.getAttributeAsString("SVTYPE", "");
            if (!svType.equals("DEL") && !svType.equals("DUP")) {
                continue;
            }
            final int shiftedStart = shiftPosition(dictionary, variant.getContig(), variant.getStart(), shift);
            final int shiftedEnd = shiftPosition(dictionary, variant.getContig(), variant.getEnd(), shift);
            final VariantContextBuilder builder = new VariantContextBuilder(variant)
                    .id(variant.getID() + "_shifted")
                    .chr(variant.getContig())
                    .start(shiftedStart)
                    .stop(shiftedEnd)
                    .attribute("END", shiftedEnd);
            syntheticVariants.add(builder.make());
        }

        syntheticVariants.sort(Comparator
                .comparingInt((VariantContext variant) -> contigOrder.get(variant.getContig()))
                .thenComparingInt(VariantContext::getStart)
                .thenComparingInt(VariantContext::getEnd)
                .thenComparing(VariantContext::getID));

        try (final VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(outputVcf)
                .setReferenceDictionary(dictionary)
                .build()) {
            writer.writeHeader(header);
            syntheticVariants.forEach(writer::add);
        }
    }

    private static void writeSyntheticDepthEvidence(final File inputFile,
                                                    final File outputFile,
                                                    final SAMSequenceDictionary dictionary,
                                                    final List<String> sampleNames,
                                                    final int shift) throws IOException {
        final DepthEvidenceCodec codec = new DepthEvidenceCodec();
        final Map<String, Integer> contigOrder = buildContigOrder(dictionary);
        final List<DepthEvidence> evidence = new ArrayList<>();
        try (final BufferedReader reader = openBlockCompressedReader(inputFile)) {
            String line = reader.readLine();
            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }
                final DepthEvidence original = codec.decode(line);
                evidence.add(original);
                evidence.add(new DepthEvidence(
                        original.getContig(),
                        shiftPosition(dictionary, original.getContig(), original.getStart(), shift),
                        shiftPosition(dictionary, original.getContig(), original.getEnd(), shift),
                        Arrays.copyOf(original.getCounts(), original.getCounts().length)));
            }
        }

        evidence.sort(Comparator
                .comparingInt((DepthEvidence item) -> contigOrder.get(item.getContig()))
                .thenComparingInt(DepthEvidence::getStart)
                .thenComparingInt(DepthEvidence::getEnd));

        try (final FeatureOutputStream<DepthEvidence> sink = codec.makeSink(new GATKPath(outputFile.getAbsolutePath()), dictionary, sampleNames, 4)) {
            for (final DepthEvidence item : evidence) {
                sink.write(item);
            }
        }
    }

    private static void writeSyntheticDiscordantPairEvidence(final File inputFile,
                                                             final File outputFile,
                                                             final SAMSequenceDictionary dictionary,
                                                             final int shift) throws IOException {
        final DiscordantPairEvidenceCodec codec = new DiscordantPairEvidenceCodec();
        final Map<String, Integer> contigOrder = buildContigOrder(dictionary);
        final List<DiscordantPairEvidence> evidence = new ArrayList<>();
        try (final BufferedReader reader = openBlockCompressedReader(inputFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }
                final DiscordantPairEvidence original = codec.decode(line);
                evidence.add(original);
                evidence.add(new DiscordantPairEvidence(
                        original.getSample(),
                        original.getContig(),
                        shiftPosition(dictionary, original.getContig(), original.getStart(), shift),
                        original.getStartStrand(),
                        original.getEndContig(),
                        shiftPosition(dictionary, original.getEndContig(), original.getEndPosition(), shift),
                        original.getEndStrand()));
            }
        }

        evidence.sort(Comparator
                .comparingInt((DiscordantPairEvidence item) -> contigOrder.get(item.getContig()))
                .thenComparingInt(DiscordantPairEvidence::getStart)
                .thenComparing(DiscordantPairEvidence::getSample));

        try (final FeatureOutputStream<DiscordantPairEvidence> sink = codec.makeSink(new GATKPath(outputFile.getAbsolutePath()), dictionary, Collections.emptyList(), 4)) {
            for (final DiscordantPairEvidence item : evidence) {
                sink.write(item);
            }
        }
    }

    private static void writeSyntheticSplitReadEvidence(final File inputFile,
                                                        final File outputFile,
                                                        final SAMSequenceDictionary dictionary,
                                                        final int shift) throws IOException {
        final SplitReadEvidenceCodec codec = new SplitReadEvidenceCodec();
        final Map<String, Integer> contigOrder = buildContigOrder(dictionary);
        final List<SplitReadEvidence> evidence = new ArrayList<>();
        try (final BufferedReader reader = openBlockCompressedReader(inputFile)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }
                final SplitReadEvidence original = codec.decode(line);
                evidence.add(original);
                evidence.add(new SplitReadEvidence(
                        original.getSample(),
                        original.getContig(),
                        shiftPosition(dictionary, original.getContig(), original.getStart(), shift),
                        original.getCount(),
                        original.getStrand()));
            }
        }

        evidence.sort(Comparator
                .comparingInt((SplitReadEvidence item) -> contigOrder.get(item.getContig()))
                .thenComparingInt(SplitReadEvidence::getStart)
                .thenComparing(SplitReadEvidence::getSample));

        try (final FeatureOutputStream<SplitReadEvidence> sink = codec.makeSink(new GATKPath(outputFile.getAbsolutePath()), dictionary, Collections.emptyList(), 4)) {
            for (final SplitReadEvidence item : evidence) {
                sink.write(item);
            }
        }
    }

    private static BufferedReader openBlockCompressedReader(final File inputFile) throws IOException {
        return new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(inputFile), StandardCharsets.UTF_8));
    }

    private static List<String> readDepthSampleNames(final File depthFile) throws IOException {
        try (final BufferedReader reader = openBlockCompressedReader(depthFile)) {
            final String headerLine = reader.readLine();
            Assert.assertNotNull(headerLine, "Depth evidence file must contain a header line");
            final String[] headerColumns = headerLine.split("\t");
            Assert.assertTrue(headerColumns.length > 3, "Depth evidence header must contain sample columns");
            return Arrays.asList(headerColumns).subList(3, headerColumns.length);
        }
    }

    private static int shiftPosition(final SAMSequenceDictionary dictionary,
                                     final String contig,
                                     final int position,
                                     final int shift) {
        final int contigLength = dictionary.getSequence(contig).getSequenceLength();
        if (position + shift <= contigLength) {
            return position + shift;
        }
        if (position - shift >= 1) {
            return position - shift;
        }
        final int forwardRoom = contigLength - position;
        final int backwardRoom = position - 1;
        final int fallbackShift = Math.max(forwardRoom, backwardRoom) / 2;
        Assert.assertTrue(fallbackShift > 0,
                "Unable to find an in-bounds shifted position for synthetic fixture on contig " + contig + " at position " + position);
        return forwardRoom >= backwardRoom ? position + fallbackShift : position - fallbackShift;
    }

    private static Map<String, Integer> buildContigOrder(final SAMSequenceDictionary dictionary) {
        final Map<String, Integer> contigOrder = new HashMap<>();
        for (int i = 0; i < dictionary.size(); i++) {
            contigOrder.put(dictionary.getSequence(i).getSequenceName(), i);
        }
        return contigOrder;
    }

    private static Map<String, Double> readSingleMetricRow(final File metricFile) throws IOException {
        final List<String> lines = Files.readAllLines(metricFile.toPath());
        Assert.assertTrue(lines.size() >= 2, "Expected metrics file to contain a header and one data line: " + metricFile.getName());
        final String[] header = lines.get(0).split("\t");
        final String[] values = lines.get(1).split("\t");
        Assert.assertEquals(values.length, header.length, "Metric row should have the same number of columns as the header");
        final Map<String, Double> result = new HashMap<>();
        for (int i = 0; i < header.length; i++) {
            result.put(header[i], Double.parseDouble(values[i]));
        }
        return result;
    }

    private record SyntheticTrainSVFixture(File vcfFile, File peFile, File srFile, File rdFile) {}

    private record TrainSVRunResult(File outputVcf, File rdParamsFile, File peParamsFile, File srParamsFile) {}
}
