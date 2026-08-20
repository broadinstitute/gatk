package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class RefineComplexVariantsUnitTest {

    private static final SAMSequenceDictionary DICTIONARY = new SAMSequenceDictionary(Arrays.asList(
            new SAMSequenceRecord("chr1", 250_000_000),
            new SAMSequenceRecord("chr2", 250_000_000),
            new SAMSequenceRecord("chr6", 250_000_000),
            new SAMSequenceRecord("chr8", 250_000_000)
    ));

    @Test
    public void testCoverageFractionUsesUnionAcrossOverlaps() {
        final List<RefineComplexVariants.DepthInterval> intervals = Arrays.asList(
                new RefineComplexVariants.DepthInterval(100, 150),
                new RefineComplexVariants.DepthInterval(140, 160),
                new RefineComplexVariants.DepthInterval(170, 200)
        );

        final double coverage = RefineComplexVariants.computeCoverageFraction(intervals, 120, 200);
        Assert.assertEquals(coverage, 0.875, 1e-10);
    }

    @Test
    public void testSinkPositionRelativeToSourceTreatsContigOrderAndTouchingIntervalsAsSeparated() {
        Assert.assertEquals(
                RefineComplexVariants.getSinkPositionRelativeToSource(
                        new SimpleInterval("chr1", 100, 200),
                        new SimpleInterval("chr2", 300, 400),
                        DICTIONARY),
                RefineComplexVariants.SinkPositionRelativeToSource.BEFORE_SOURCE);
        Assert.assertEquals(
                RefineComplexVariants.getSinkPositionRelativeToSource(
                        new SimpleInterval("chr8", 100, 200),
                        new SimpleInterval("chr6", 300, 400),
                        DICTIONARY),
                RefineComplexVariants.SinkPositionRelativeToSource.AFTER_SOURCE);
        Assert.assertEquals(
                RefineComplexVariants.getSinkPositionRelativeToSource(
                        new SimpleInterval("chr1", 100, 200),
                        new SimpleInterval("chr1", 200, 300),
                        DICTIONARY),
                RefineComplexVariants.SinkPositionRelativeToSource.BEFORE_SOURCE);
        Assert.assertEquals(
                RefineComplexVariants.getSinkPositionRelativeToSource(
                        new SimpleInterval("chr1", 300, 400),
                        new SimpleInterval("chr1", 100, 300),
                        DICTIONARY),
                RefineComplexVariants.SinkPositionRelativeToSource.AFTER_SOURCE);
    }

    @Test
    public void testSinkPositionRelativeToSourceDetectsOverlappingSinkAsWithinSource() {
        Assert.assertEquals(
                RefineComplexVariants.getSinkPositionRelativeToSource(
                        new SimpleInterval("chr1", 150, 250),
                        new SimpleInterval("chr1", 100, 300),
                        DICTIONARY),
                RefineComplexVariants.SinkPositionRelativeToSource.WITHIN_SOURCE);
    }

    @Test
    public void testIntrachromosomalPlanForSmallDispersedDuplicationUsesSinkBreakpoints() {
        final VariantContext variant = new VariantContextBuilder("test", "chr1", 1000, 1000,
                Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                .id("cpx-small-ddup")
                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(GATKSVVCFConstants.CPX_TYPE, "dDUP")
                .attribute(GATKSVVCFConstants.CPX_INTERVALS, Collections.singletonList("DUP_chr1:2000-2100"))
                .make();

        final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
        Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CPX);
        Assert.assertEquals(plan.queries.size(), 2);

        final RefineComplexVariants.DiscordantPairQuery firstQuery = plan.queries.get(0);
        final RefineComplexVariants.DiscordantPairQuery secondQuery = plan.queries.get(1);
        Assert.assertTrue(firstQuery.startStrand);
        Assert.assertFalse(firstQuery.endStrand);
        Assert.assertEquals(firstQuery.startMin, 0);
        Assert.assertEquals(firstQuery.startMax, 1100);
        Assert.assertEquals(firstQuery.endLowerExclusive, 1900);
        Assert.assertEquals(firstQuery.endUpperExclusive, 3000);
        Assert.assertFalse(secondQuery.startStrand);
        Assert.assertTrue(secondQuery.endStrand);
        Assert.assertEquals(secondQuery.startMin, 900);
        Assert.assertEquals(secondQuery.startMax, 2000);
        Assert.assertEquals(secondQuery.endLowerExclusive, 1100);
        Assert.assertEquals(secondQuery.endUpperExclusive, 2200);
    }

    @Test
    public void testCarrierDecisionMatchesWorkflowOutcome() {
        final RefineComplexVariants.CarrierRefinement partialWithoutDepth = RefineComplexVariants.evaluateCarrierSupport(
                4, 0, 3, Collections.emptyList());
        Assert.assertTrue(partialWithoutDepth.reviseGenotype);
        Assert.assertTrue(partialWithoutDepth.countTowardsUnresolved);

        final RefineComplexVariants.CarrierRefinement lowPeRescuedByDepth = RefineComplexVariants.evaluateCarrierSupport(
                2, 2, 3, Collections.singletonList(Boolean.TRUE));
        Assert.assertFalse(lowPeRescuedByDepth.reviseGenotype);
        Assert.assertFalse(lowPeRescuedByDepth.countTowardsUnresolved);

        final RefineComplexVariants.CarrierRefinement highPeWithDepthFailure = RefineComplexVariants.evaluateCarrierSupport(
                4, 5, 3, Arrays.asList(Boolean.TRUE, Boolean.FALSE));
        Assert.assertTrue(highPeWithDepthFailure.reviseGenotype);
        Assert.assertFalse(highPeWithDepthFailure.countTowardsUnresolved);
    }

    @Test
    public void testPlanForInvDup() {
        final VariantContext variant = new VariantContextBuilder("test", "chr1", 7443032, 7443928,
                Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                .id("cpx1")
                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(GATKSVVCFConstants.CPX_TYPE, "INVdup")
                .attribute(GATKSVVCFConstants.CPX_INTERVALS,
                        Arrays.asList("INV_chr1:7443032-7443928", "DUP_chr1:7443630-7443928"))
                .make();

        final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
        Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CPX);
        Assert.assertEquals(plan.queries.size(), 2);
        Assert.assertTrue(plan.requiresDepthEvidence);

        final RefineComplexVariants.DiscordantPairQuery firstQuery = plan.queries.get(0);
        final RefineComplexVariants.DiscordantPairQuery secondQuery = plan.queries.get(1);
        Assert.assertTrue(firstQuery.startStrand);
        Assert.assertTrue(firstQuery.endStrand);
        Assert.assertFalse(secondQuery.startStrand);
        Assert.assertFalse(secondQuery.endStrand);
    }

    @Test
    public void testPlanForDelInvDupUsesDupStartBreakpoint() {
            final VariantContext variant = new VariantContextBuilder("test", "chr1", 157917449, 157918126,
                            Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                            .id("cpx-delINVdup")
                            .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                            .attribute(GATKSVVCFConstants.CPX_TYPE, "delINVdup")
                            .attribute(GATKSVVCFConstants.CPX_INTERVALS,
                                            Arrays.asList(
                                                            "DEL_chr1:157917449-157917619",
                                                            "INV_chr1:157917619-157918126",
                                                            "DUP_chr1:157917760-157918126"))
                            .make();

            final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
            Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CPX);
            Assert.assertEquals(plan.queries.size(), 2);
            Assert.assertEquals(plan.queries.get(1).endLowerExclusive, 157917660);
            Assert.assertEquals(plan.queries.get(1).endUpperExclusive, 157918760);
    }

    @Test
    public void testCtxPqQpPlanUsesBreakpointQueryWindowsInContigOrder() {
            final VariantContext variant = new VariantContextBuilder("test", "chr8", 7000, 7200,
                            Arrays.asList(Allele.REF_N, Allele.create("<CTX>", false)))
                            .id("ctx-pq-qp")
                            .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX)
                    .attribute(GATKSVVCFConstants.CPX_TYPE, "CTX_PQ/QP")
                            .attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, "chr2")
                            .attribute(GATKSVVCFConstants.END2_ATTRIBUTE, 4000)
                            .make();

            final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
            Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CTX);
            Assert.assertEquals(plan.queries.size(), 2);

            final RefineComplexVariants.DiscordantPairQuery firstQuery = plan.queries.get(0);
            final RefineComplexVariants.DiscordantPairQuery secondQuery = plan.queries.get(1);
            Assert.assertEquals(firstQuery.startContig, "chr2");
            Assert.assertEquals(firstQuery.endContig, "chr8");
            Assert.assertTrue(firstQuery.startStrand);
            Assert.assertTrue(firstQuery.endStrand);
            Assert.assertEquals(firstQuery.startMin, 3000);
            Assert.assertEquals(firstQuery.startMax, 4100);
            Assert.assertEquals(firstQuery.endLowerExclusive, 6000);
            Assert.assertEquals(firstQuery.endUpperExclusive, 7100);
            Assert.assertFalse(secondQuery.startStrand);
            Assert.assertFalse(secondQuery.endStrand);
            Assert.assertEquals(secondQuery.startMin, 3900);
            Assert.assertEquals(secondQuery.startMax, 5000);
            Assert.assertEquals(secondQuery.endLowerExclusive, 7100);
            Assert.assertEquals(secondQuery.endUpperExclusive, 8200);
    }

    @Test
    public void testCtxPpQqPlanUsesOpposingStrandsAtOrderedBreakpoints() {
            final VariantContext variant = new VariantContextBuilder("test", "chr1", 7000, 7200,
                            Arrays.asList(Allele.REF_N, Allele.create("<CTX>", false)))
                            .id("ctx-pp-qq")
                            .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.StructuralVariantAnnotationType.CTX)
                    .attribute(GATKSVVCFConstants.CPX_TYPE, "CTX_PP/QQ")
                            .attribute(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, "chr6")
                            .attribute(GATKSVVCFConstants.END2_ATTRIBUTE, 4000)
                            .make();

            final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
            Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CTX);
            Assert.assertEquals(plan.queries.size(), 2);

            final RefineComplexVariants.DiscordantPairQuery firstQuery = plan.queries.get(0);
            final RefineComplexVariants.DiscordantPairQuery secondQuery = plan.queries.get(1);
            Assert.assertEquals(firstQuery.startContig, "chr1");
            Assert.assertEquals(firstQuery.endContig, "chr6");
            Assert.assertTrue(firstQuery.startStrand);
            Assert.assertFalse(firstQuery.endStrand);
            Assert.assertEquals(firstQuery.startMin, 6000);
            Assert.assertEquals(firstQuery.startMax, 7100);
            Assert.assertEquals(firstQuery.endLowerExclusive, 3900);
            Assert.assertEquals(firstQuery.endUpperExclusive, 5000);
            Assert.assertFalse(secondQuery.startStrand);
            Assert.assertTrue(secondQuery.endStrand);
            Assert.assertEquals(secondQuery.startMin, 7100);
            Assert.assertEquals(secondQuery.startMax, 8200);
            Assert.assertEquals(secondQuery.endLowerExclusive, 3000);
            Assert.assertEquals(secondQuery.endUpperExclusive, 4100);
    }

    @Test
    public void testInterchromosomalPlanForDispersedDuplication() {
        final VariantContext variant = new VariantContextBuilder("test", "chr1", 1000, 1200,
                Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                .id("cpx2")
                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(GATKSVVCFConstants.CPX_TYPE, "dDUP")
                .attribute(GATKSVVCFConstants.CPX_INTERVALS, Collections.singletonList("DUP_chr2:3000-3200"))
                .make();

        final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
        Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CPX);
        Assert.assertEquals(plan.queries.size(), 2);

        final RefineComplexVariants.DiscordantPairQuery firstQuery = plan.queries.get(0);
        Assert.assertEquals(firstQuery.startContig, "chr1");
        Assert.assertEquals(firstQuery.endContig, "chr2");
        Assert.assertTrue(firstQuery.startStrand);
        Assert.assertFalse(firstQuery.endStrand);
    }

    @Test
    public void testInterchromosomalPlanForDispersedDuplicationWhenSinkComesAfterSource() {
            final VariantContext variant = new VariantContextBuilder("test", "chr8", 7000, 7200,
                            Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                            .id("cpx2-after")
                            .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                            .attribute(GATKSVVCFConstants.CPX_TYPE, "dDUP")
                            .attribute(GATKSVVCFConstants.CPX_INTERVALS, Collections.singletonList("DUP_chr6:3000-3200"))
                            .make();

            final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
            Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CPX);
            Assert.assertEquals(plan.queries.size(), 2);

            final RefineComplexVariants.DiscordantPairQuery firstQuery = plan.queries.get(0);
            Assert.assertEquals(firstQuery.startContig, "chr6");
            Assert.assertEquals(firstQuery.endContig, "chr8");
            Assert.assertTrue(firstQuery.startStrand);
            Assert.assertFalse(firstQuery.endStrand);
            Assert.assertEquals(firstQuery.startMin, 2200);
            Assert.assertEquals(firstQuery.endLowerExclusive, 7100);
    }


    @Test
    public void testDDupIdelPlanUsesDeletionIntervalFromCpxIntervals() {
        final VariantContext variant = new VariantContextBuilder("test", "chr1", 1000, 1000,
                Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                .id("cpx3")
                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(GATKSVVCFConstants.CPX_TYPE, "dDUP_iDEL")
                .attribute(GATKSVVCFConstants.CPX_INTERVALS,
                        Arrays.asList("DUP_chr1:2000-2100", "DEL_chr1:5000-6500"))
                .make();

        final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
        Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CPX);
        Assert.assertEquals(plan.queries.size(), 2);
        Assert.assertEquals(plan.queries.get(0).endLowerExclusive, 6400);
        Assert.assertEquals(plan.queries.get(1).endLowerExclusive, 4000);
    }

    @Test
    public void testInsIdelPlanUsesDeletionIntervalFromCpxIntervals() {
        final VariantContext variant = new VariantContextBuilder("test", "chr1", 1000, 1000,
                Arrays.asList(Allele.REF_N, Allele.create("<CPX>", false)))
                .id("cpx4")
                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.CPX_SV_SYB_ALT_ALLELE_STR)
                .attribute(GATKSVVCFConstants.CPX_TYPE, "INS_iDEL")
                .attribute(GATKSVVCFConstants.CPX_INTERVALS, Collections.singletonList("DEL_chr1:5000-6500"))
                .attribute("SOURCE", "INS_chr1:2000-2100")
                .make();

        final RefineComplexVariants.EvaluationPlan plan = RefineComplexVariants.createEvaluationPlan(variant, DICTIONARY);
        Assert.assertEquals(plan.variantType, RefineComplexVariants.EvaluatedVariantType.CPX);
        Assert.assertEquals(plan.queries.size(), 2);
        Assert.assertEquals(plan.queries.get(0).endLowerExclusive, 6400);
        Assert.assertEquals(plan.queries.get(1).endLowerExclusive, 4000);
    }

    @Test
    public void testApplyFormattingTransformForInsertionWithInvertedSource() {
        final VariantContext variant = new VariantContextBuilder("test", "chr1", 6633858, 6633900,
                Arrays.asList(Allele.REF_N, Allele.create("<INS>", false)))
                .id("ins1")
                .attribute(GATKSVVCFConstants.SVTYPE, GATKSVVCFConstants.SYMB_ALT_STRING_INS)
                .attribute(SOURCE_ATTRIBUTE, "INV_chr6:34880888-34881035")
                .make();
        final Map<String, Object> attributes = new LinkedHashMap<>(variant.getAttributes());

        final Allele transform = RefineComplexVariants.applyFormattingTransform(variant, attributes);
        Assert.assertEquals(transform, Allele.create("<CPX>", false));
        Assert.assertEquals(attributes.get(GATKSVVCFConstants.SVTYPE), GATKSVVCFConstants.StructuralVariantAnnotationType.CPX);
        Assert.assertEquals(attributes.get(GATKSVVCFConstants.CPX_TYPE), "dDUP");
        Assert.assertEquals(attributes.get(GATKSVVCFConstants.CPX_INTERVALS),
                Arrays.asList("DUP_chr6:34880888-34881035", "INV_chr6:34880888-34881035"));
        Assert.assertFalse(attributes.containsKey(SOURCE_ATTRIBUTE));
    }

    private static final String SOURCE_ATTRIBUTE = "SOURCE";
}
