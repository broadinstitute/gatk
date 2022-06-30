package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.ClusteringParameters;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.concordance.ClosestSVFinder;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceAnnotator;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceLinkage;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.utils.MathUtilsUnitTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class SVConcordanceIntegrationTest extends CommandLineProgramTest {

    private static final double TOLERANCE = 0.001;

    @Test
    public void testRefPanel() {
        final File output = createTempFile("concord", ".vcf");
        final String evalVcfPath = getToolTestDataDir() + "ref_panel_1kg.cleaned.gatk.supported_svtypes.chr22_chrY.vcf.gz";
        /**
         * Test file produced from raw standardized VCFs from manta, melt, wham, and cnv callers with SVCluster parameters:
         * --algorithm SINGLE_LINKAGE
         * --pesr-interval-overlap 0.9
         * --pesr-breakend-window 50
         */
        final String truthVcfPath = getToolTestDataDir() + "ref_panel_1kg.raw_calls.chr22_chrY.vcf.gz";

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addOutput(output)
                .add(StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME, GATKBaseTest.FULL_HG38_DICT)
                .add(SVClusterEngineArgumentsCollection.DEPTH_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.DEPTH_INTERVAL_OVERLAP_FRACTION_NAME, 0.5)
                .add(SVClusterEngineArgumentsCollection.DEPTH_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.MIXED_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.MIXED_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.MIXED_BREAKEND_WINDOW_NAME, 2000)
                .add(SVClusterEngineArgumentsCollection.PESR_SAMPLE_OVERLAP_FRACTION_NAME, 0)
                .add(SVClusterEngineArgumentsCollection.PESR_INTERVAL_OVERLAP_FRACTION_NAME, 0.1)
                .add(SVClusterEngineArgumentsCollection.PESR_BREAKEND_WINDOW_NAME, 500)
                .add(SVConcordance.BIALLELIC_DUPLICATIONS_LONG_NAME, true)
                .add(AbstractConcordanceWalker.TRUTH_VARIANTS_LONG_NAME, truthVcfPath)
                .add(AbstractConcordanceWalker.EVAL_VARIANTS_SHORT_NAME, evalVcfPath);


        runCommandLine(args, SVConcordance.class.getSimpleName());

        final Pair<VCFHeader, List<VariantContext>> outputVcf = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        final SAMSequenceDictionary dictionary = SVTestUtils.hg38Dict;
        final ClusteringParameters depthParameters = ClusteringParameters.createDepthParameters(0.5, 2000, 0);
        final ClusteringParameters mixedParameters = ClusteringParameters.createMixedParameters(0.1, 2000, 0);
        final ClusteringParameters pesrParameters = ClusteringParameters.createPesrParameters(0.1, 500, 0);
        final SVConcordanceLinkage linkage = new SVConcordanceLinkage(dictionary);
        linkage.setDepthOnlyParams(depthParameters);
        linkage.setMixedParams(mixedParameters);
        linkage.setEvidenceParams(pesrParameters);
        final SVConcordanceAnnotator annotator = new SVConcordanceAnnotator(false);

        final List<SVCallRecord> inputEvalVariants = VariantContextTestUtils.readEntireVCFIntoMemory(evalVcfPath).getValue()
                .stream().map(SVCallRecordUtils::create).map(r -> SVCallRecordUtils.convertToBiallelicDupGenotypes(r, null)).collect(Collectors.toList());
        final List<SVCallRecord> inputTruthVariants = VariantContextTestUtils.readEntireVCFIntoMemory(truthVcfPath).getValue()
                .stream().map(SVCallRecordUtils::create).map(r -> SVCallRecordUtils.convertToBiallelicDupGenotypes(r, null)).collect(Collectors.toList());

        final ClosestSVFinder finder = new ClosestSVFinder(linkage, annotator::annotate, dictionary);
        final List<SVCallRecord> expectedRecords = new ArrayList<>(inputEvalVariants.size());
        final Set<Map.Entry<Long, SVCallRecord>> truthRecords = new HashSet<>();
        Long itemId = 0L;
        for (final SVCallRecord record : inputTruthVariants) {
            truthRecords.add(new AbstractMap.SimpleImmutableEntry<>(itemId++, record));
        }
        for (final SVCallRecord evalRecord : inputEvalVariants) {
            final Map.Entry<Long, SVCallRecord> closestRecord = finder.getClosestItem(evalRecord, truthRecords);
            final ClosestSVFinder.ClosestPair cluster =
                    new ClosestSVFinder.ClosestPair(evalRecord, closestRecord == null ? null : closestRecord.getValue());
            expectedRecords.add(annotator.annotate(cluster));
        }

        final Comparator<VariantContext> idComparator = Comparator.comparing(VariantContext::getID);
        final List<VariantContext> expectedVariants = expectedRecords.stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .sorted(idComparator)
                .collect(Collectors.toList());
        final List<VariantContext> outputVariants = outputVcf.getValue().stream().sorted(idComparator).collect(Collectors.toList());
        final List<VariantContext> idSortedInputEvalVariants = inputEvalVariants.stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .sorted(idComparator)
                .collect(Collectors.toList());

        Assert.assertEquals(outputVariants.size(), idSortedInputEvalVariants.size());
        Assert.assertEquals(outputVariants.size(), expectedVariants.size());
        final Set<String> checkedVariantsSet = new HashSet<>();
        for (int i = 0; i < outputVariants.size(); i++) {
            final VariantContext outputVariant = outputVariants.get(i);
            final VariantContext expectedVariant = expectedVariants.get(i);
            Assert.assertEquals(outputVariant.getID(), idSortedInputEvalVariants.get(i).getID());
            Assert.assertEquals(outputVariant.getID(), expectedVariant.getID());
            Assert.assertEquals(outputVariant.getContig(), expectedVariant.getContig());
            Assert.assertEquals(outputVariant.getStart(), expectedVariant.getStart());
            Assert.assertEquals(outputVariant.getEnd(), expectedVariant.getEnd());
            Assert.assertEquals(outputVariant.getAlleles(), expectedVariant.getAlleles());
            Assert.assertEquals(outputVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "test_default"), expectedVariant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, "expected_default"));
            checkTruthVariantId(outputVariant, expectedVariant);
            Assert.assertEquals(outputVariant.getAttributeAsString(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, "test_default"), expectedVariant.getAttributeAsString(Concordance.TRUTH_STATUS_VCF_ATTRIBUTE, "expected_default"));
            if (outputVariant.getID().equals("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_1")) {
                checkVariant(outputVariant, expectedVariant, "ref_panel_1kg_raw_00000062",
                        0.8461538461538461, 0.21794871794871795, checkedVariantsSet);
            } else if (outputVariant.getID().equals("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_140")) {
                checkVariant(outputVariant, expectedVariant, "ref_panel_1kg_raw_0000184c",
                        0.9230769230769231, 0.122, checkedVariantsSet);
            } else if (outputVariant.getID().equals("ref_panel_1kg.chr22.final_cleanup_DUP_chr22_152")) {
                checkVariant(outputVariant, expectedVariant, "ref_panel_1kg_raw_00002321",
                        0.717948717948718, 0.256, checkedVariantsSet);
            } else if (outputVariant.getID().equals("ref_panel_1kg.chr22.final_cleanup_INS_chr22_100")) {
                checkVariant(outputVariant, expectedVariant, "ref_panel_1kg_raw_00001ba0",
                        0.9871794871794872, 0.077, checkedVariantsSet);
            }
        }
        Assert.assertEquals(checkedVariantsSet.size(), 4);
        Assert.assertTrue(checkedVariantsSet.contains("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_1"));
        Assert.assertTrue(checkedVariantsSet.contains("ref_panel_1kg.chr22.final_cleanup_DEL_chr22_140"));
        Assert.assertTrue(checkedVariantsSet.contains("ref_panel_1kg.chr22.final_cleanup_DUP_chr22_152"));
        Assert.assertTrue(checkedVariantsSet.contains("ref_panel_1kg.chr22.final_cleanup_INS_chr22_100"));
    }

    private static void checkTruthVariantId(final VariantContext actual, final VariantContext expected) {
        final String expectedTruthId = expected.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, VCFConstants.EMPTY_INFO_FIELD);
        final String actualTruthId = actual.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, VCFConstants.EMPTY_INFO_FIELD);
        Assert.assertEquals(actualTruthId, expectedTruthId);
    }

    private static void checkVariant(final VariantContext actual, final VariantContext expected,
                                     final String expectedTruthId, final double expectedGenotypeConcordance,
                                     final double expectedAF,
                                     final Set<String> checkedVariantsSet) {
        Assert.assertFalse(checkedVariantsSet.contains(actual.getID()));
        checkedVariantsSet.add(actual.getID());
        Assert.assertEquals(actual.getAttributeAsString(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, ""), expectedTruthId);
        MathUtilsUnitTest.assertEqualsDoubleSmart(actual.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), expectedGenotypeConcordance, TOLERANCE);
        MathUtilsUnitTest.assertEqualsDoubleSmart(actual.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), expected.getAttributeAsDouble(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, -1.), TOLERANCE);
        MathUtilsUnitTest.assertEqualsDoubleSmart(actual.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.), expectedAF, TOLERANCE);
    }
}
