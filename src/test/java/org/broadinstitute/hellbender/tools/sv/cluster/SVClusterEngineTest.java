package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
import static org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE;

public class SVClusterEngineTest {

    private final SVClusterEngine engine = SVTestUtils.defaultSingleLinkageEngine;

    private static final ClusteringParameters depthOnlyParametersSizeSimilarity = ClusteringParameters.createDepthParameters(0.1, 0.5, 0, 0);
    private static final ClusteringParameters mixedParametersSizeSimilarity = ClusteringParameters.createMixedParameters(0.1, 0.5, 5000, 0);
    private static final ClusteringParameters evidenceParametersSizeSimilarity = ClusteringParameters.createPesrParameters(0.1, 0.5, 5000, 0);
    private final CanonicalSVLinkage<SVCallRecord> linkageSizeSimilarity = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, false);
    private final SVClusterEngine engineSizeSimilarity = new SVClusterEngine(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE, SVTestUtils.defaultCollapser::collapse, linkageSizeSimilarity, SVTestUtils.hg38Dict);

    // Assign length according to coordinate scheme
    private static Integer inferLength(final String contigA, final int posA, final String contigB, final int posB) {
        if (contigA.equals(contigB)) {
            if (posA == posB) {
                // no interval
                return null;
            } else {
                // intrachromosomal and intervaled
                return posB - posA;
            }
        } else {
            // interchromosomal
            return null;
        }
    }

    @BeforeTest
    public void initializeClusterEngine() {
        engine.add(SVTestUtils.call1);
        linkageSizeSimilarity.setDepthOnlyParams(depthOnlyParametersSizeSimilarity);
        linkageSizeSimilarity.setMixedParams(mixedParametersSizeSimilarity);
        linkageSizeSimilarity.setEvidenceParams(evidenceParametersSizeSimilarity);
    }

    @Test
    public void testCollapser() {
        //depth only and depthAndStuff have same bounds, less than call2
        final SVClusterEngine.OutputCluster testCluster = new SVClusterEngine.OutputCluster(Arrays.asList(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff, SVTestUtils.call2));
        final SVCallRecord flattened = engine.getCollapser().apply(testCluster);
        Assert.assertEquals(flattened.getPositionA(), SVTestUtils.depthAndStuff.getPositionA());
        Assert.assertEquals(flattened.getPositionB(), SVTestUtils.depthAndStuff.getPositionB());
        //should have all the algs
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthAndStuff.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.depthOnly.getAlgorithms()));
        Assert.assertTrue(flattened.getAlgorithms().containsAll(SVTestUtils.call2.getAlgorithms()));

        SVTestUtils.assertContainsAllIgnoreRefAlleleBase(flattened.getGenotypes(), SVTestUtils.depthAndStuff.getGenotypes(), true);
        SVTestUtils.assertContainsAllIgnoreRefAlleleBase(flattened.getGenotypes(), SVTestUtils.depthOnly.getGenotypes(), true);
        SVTestUtils.assertContainsAllIgnoreRefAlleleBase(flattened.getGenotypes(), SVTestUtils.call2.getGenotypes(), true);

        // Test subtyped insertions
        final List<SVCallRecord> subtypedRecords = Lists.newArrayList(
                SVTestUtils.newCallRecordWithAlleles(
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                        Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_INS),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        2,
                        null),
                SVTestUtils.newCallRecordWithAlleles(
                        Lists.newArrayList(Allele.REF_N, Allele.create("<INS:MEI>")),
                        Lists.newArrayList(Allele.REF_N, Allele.create("<INS:MEI>")),
                        GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                        2,
                        null)
        );
        final SVCallRecord subtypedFlattened = engine.getCollapser().apply(new SVClusterEngine.OutputCluster(subtypedRecords));
        Assert.assertEquals(subtypedFlattened.getAlleles().size(), 2);
        Assert.assertEquals(subtypedFlattened.getAltAlleles().size(), 1);
        final List<Allele> collapsedAlleles = subtypedFlattened.getGenotypes().stream().map(Genotype::getAlleles)
                .flatMap(Collection::stream).distinct().sorted().collect(Collectors.toList());
        Assert.assertEquals(subtypedFlattened.getAlleles(), collapsedAlleles);
    }

    @Test
    public void testClusterTogether() {
        Assert.assertTrue(engine.getLinkage().areClusterable(SVTestUtils.depthOnly, SVTestUtils.depthAndStuff));
        Assert.assertFalse(engine.getLinkage().areClusterable(SVTestUtils.depthOnly, SVTestUtils.inversion));
        Assert.assertFalse(engine.getLinkage().areClusterable(SVTestUtils.call1, SVTestUtils.call2));
        Assert.assertTrue(engine.getLinkage().areClusterable(SVTestUtils.call1, SVTestUtils.overlapsCall1));
    }

    @DataProvider(name = "testClusterTogetherWithSizeSimilarityDataProvider")
    public Object[][] testClusterTogetherWithSizeSimilarityDataProvider() {
        return new Object[][]{
                {1001, 2000, 1001, 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, true},
                {1001, 2000, 1501, 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, true},
                {1001, 2000, 1900, 2900, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, true},
                // Fails size similarity
                {1001, 2000, 1502, 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, false},
                {1001, 2000, 1901, 2399, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, false},
                // Fails reciprocal overlap
                {1001, 2000, 1902, 2900, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, false}
        };
    }

    @Test(dataProvider= "testClusterTogetherWithSizeSimilarityDataProvider")
    public void testClusterTogetherWithSizeSimilarity(final int start1, final int end1,
                                                      final int start2, final int end2,
                                                      final GATKSVVCFConstants.StructuralVariantAnnotationType svtype,
                                                      final boolean expectedResult) {
        // PESR
        final SVCallRecord record1 = SVTestUtils.newPESRCallRecordWithIntervalAndType(start1, end1, svtype);
        final SVCallRecord record2 = SVTestUtils.newPESRCallRecordWithIntervalAndType(start2, end2, svtype);
        Assert.assertEquals(engineSizeSimilarity.getLinkage().areClusterable(record1, record2), expectedResult);
        Assert.assertEquals(engineSizeSimilarity.getLinkage().areClusterable(record2, record1), expectedResult);

        // Depth only
        final SVCallRecord record3 = SVTestUtils.newDepthCallRecordWithIntervalAndType(start1, end1, svtype);
        final SVCallRecord record4 = SVTestUtils.newDepthCallRecordWithIntervalAndType(start2, end2, svtype);
        Assert.assertEquals(engineSizeSimilarity.getLinkage().areClusterable(record3, record4), expectedResult);
        Assert.assertEquals(engineSizeSimilarity.getLinkage().areClusterable(record4, record3), expectedResult);
    }

    @DataProvider(name = "testClusterTogetherWithSizeSimilarityInsertionsDataProvider")
    public Object[][] testClusterTogetherWithSizeSimilarityInsertionsDataProvider() {
        return new Object[][]{
                {1000, 1000, 1000, 1000, true},
                {1000, 1000, 1000, 500, true},
                {1000, 1000, 1000, 2000, true},
                {1000, 1000, 100, 1000, true},
                {1000, 1000, 1900, 1000, true},
                // Fails reciprocal overlap
                {1000, 1000, 99, 1000, false},
                {1000, 1000, 1901, 1000, false},
                // Fails size similarity
                {1000, 1000, 1000, 499, false},
                {1000, 2001, 1000, 1000, false},
        };
    }

    @Test(dataProvider= "testClusterTogetherWithSizeSimilarityInsertionsDataProvider")
    public void testClusterTogetherWithSizeSimilarityInsertions(final int start1, final int length1,
                                                                final int start2, final int length2,
                                                                final boolean expectedResult) {
        final SVCallRecord record1 = SVTestUtils.newInsertionWithPositionAndLength(start1, length1);
        final SVCallRecord record2 = SVTestUtils.newInsertionWithPositionAndLength(start2, length2);
        Assert.assertEquals(engineSizeSimilarity.getLinkage().areClusterable(record1, record2), expectedResult);
        Assert.assertEquals(engineSizeSimilarity.getLinkage().areClusterable(record2, record1), expectedResult);
    }

    @Test(expectedExceptions = { IllegalArgumentException.class })
    public void testClusterTogetherInvalidInterval() {
        // End position beyond contig end after padding
        final SVCallRecord deletion1 = new SVCallRecord("test_del", "chr1", 1000, true, "chr1", 248956423 + SVTestUtils.defaultEvidenceParameters.getWindow(), false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                null, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord deletion2 = new SVCallRecord("test_del", "chr1", 1000, true, "chr1", 248956422, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                null, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        engine.getLinkage().areClusterable(deletion1, deletion2);
        Assert.fail("Expected exception not thrown");
    }

    @Test
    public void testGetClusteringIntervalEdge() {
        //edge case - end of contig
        Assert.assertTrue(engine.getLinkage().getMaxClusterableStartingPosition(SVTestUtils.rightEdgeCall) <= SVTestUtils.chr1Length);
    }

    @DataProvider(name = "maxPositionIntervals")
    public Object[][] recordPairs() {
        return new Object[][]{
                {100, 200},
                {50000, 500000},
                {1, 1},
                {1, 2}
        };
    }

    @Test(dataProvider= "maxPositionIntervals")
    public void testGetMaxClusterableStartingPosition(final int start, final int end) {
        testGetMaxClusterableStartingPositionWithAlgorithm(start, end, GATKSVVCFConstants.DEPTH_ALGORITHM);
        testGetMaxClusterableStartingPositionWithAlgorithm(start, end, SVTestUtils.PESR_ALGORITHM);
    }

    private void testGetMaxClusterableStartingPositionWithAlgorithm(final int start, final int end, final String algorithm) {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", start, true, "chr1", end, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                end - start + 1, Collections.singletonList(algorithm),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final int maxClusterableStart = engine.getLinkage().getMaxClusterableStartingPosition(call1);

        final int call2Start = maxClusterableStart;
        final SVCallRecord call2Depth = new SVCallRecord("call2", "chr1", call2Start, true, "chr1", call2Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2Pesr = new SVCallRecord("call2", "chr1", call2Start, true, "chr1", call2Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertTrue(engine.getLinkage().areClusterable(call1, call2Depth) || engine.getLinkage().areClusterable(call1, call2Pesr));

        final int call3Start = maxClusterableStart + 1;
        final SVCallRecord call3Depth = new SVCallRecord("call2", "chr1", call3Start, true, "chr1", call3Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call3Pesr = new SVCallRecord("call2", "chr1", call3Start, true, "chr1", call3Start + call1.getLength() - 1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, Collections.emptyList(),
                call1.getLength(), SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertFalse(engine.getLinkage().areClusterable(call1, call3Depth) || engine.getLinkage().areClusterable(call1, call3Pesr));
    }

    @Test
    public void testGetClusteringIntervalLists() {
        //test overloaded function with List
        final List<SVCallRecord> pesrClusterList = new ArrayList<>();
        pesrClusterList.add(SVTestUtils.depthAndStuff);
        final int pesrCluster1 = engine.getLinkage().getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster1 >= SVTestUtils.depthAndStuff.getPositionA() + SVTestUtils.defaultEvidenceParameters.getWindow());
        Assert.assertTrue(pesrCluster1 >= SVTestUtils.depthAndStuff.getPositionA() + SVTestUtils.defaultMixedParameters.getWindow());
        //add an upstream variant
        pesrClusterList.add(SVTestUtils.depthAndStuff2);
        final int pesrCluster2 = engine.getLinkage().getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster2 >= pesrCluster1);
        //add a downstream variant
        pesrClusterList.add(SVTestUtils.depthAndStuff3);
        final int pesrCluster3 = engine.getLinkage().getMaxClusterableStartingPosition(pesrClusterList);
        Assert.assertTrue(pesrCluster3 >= SVTestUtils.depthAndStuff3.getPositionA() + SVTestUtils.defaultEvidenceParameters.getWindow());
        Assert.assertTrue(pesrCluster3 >= SVTestUtils.depthAndStuff3.getPositionA() + SVTestUtils.defaultMixedParameters.getWindow());
    }

    @Test
    public void testIsDepthOnlyCall() {
        Assert.assertTrue(SVTestUtils.call1.isDepthOnly());
        Assert.assertFalse(SVTestUtils.depthAndStuff.isDepthOnly());
        Assert.assertFalse(SVTestUtils.inversion.isDepthOnly());
    }

    @DataProvider(name = "clusterTogetherVaryPositionsProvider")
    public Object[][] clusterTogetherVaryPositionsProvider() {
        return new Object[][]{
                {500, 1001, 1001, 1502, false},  // abutting
                {500, 1001, 500, 1001, true},  // exactly equal
                {500, 1001, 600, 1101, true},  // barely meet reciprocal overlap
                {500, 1001, 601, 1102, false}, // call2 shifted slightly up
                {500, 1000, 600, 1101, false}, // call1 slightly larger
                {500, 501, 500, 501, true}, // tiny but equal
                {500, 501, 500, 502, false}, // tiny but call2 twice as big
                {500, 500, 500, 500, true}, // 0-length and equal
                {500, 500, 501, 501, false}, // 0-length and not equal
                {1, SVTestUtils.chr1Length, 1, SVTestUtils.chr1Length, true}, // really big
                {1, 10001, 1, 10001, true}, // left contig edge
                {SVTestUtils.chr1Length - 10000, SVTestUtils.chr1Length, SVTestUtils.chr1Length - 10000, SVTestUtils.chr1Length, true}, // right contig edge
                {100000, 200000, 101001, 201001, false}, // window test fail
                {100000, 200000, 101000, 201000, true} // window test success
        };
    }

    @Test(dataProvider= "clusterTogetherVaryPositionsProvider")
    public void testClusterTogetherVaryPositions(final int start1, final int end1, final int start2, final int end2, final boolean result) {
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", start1, true,
                "chr1", end1, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), end1 - start1 + 1, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                SVTestUtils.threeGenotypes, Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", start2, true,
                "chr1", end2, false,
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(), end2 - start2 + 1, Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                SVTestUtils.threeGenotypes, Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertEquals(engine.getLinkage().areClusterable(call1, call2), result);
    }

    @Test
    public void testClusterTogetherVaryTypes() {
        for (final GATKSVVCFConstants.StructuralVariantAnnotationType type1 : GATKSVVCFConstants.StructuralVariantAnnotationType.values()) {
            // Pass in null strands to let them be determined automatically
            final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, SVTestUtils.getValidTestStrandA(type1),
                    "chr1", 2001, SVTestUtils.getValidTestStrandB(type1), type1, null, Collections.emptyList(),
                    SVTestUtils.getLength(1000, 2001, type1), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                    Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
            for (final GATKSVVCFConstants.StructuralVariantAnnotationType type2 : GATKSVVCFConstants.StructuralVariantAnnotationType.values()) {
                final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, SVTestUtils.getValidTestStrandA(type2),
                        "chr1", 2001, SVTestUtils.getValidTestStrandB(type2), type2, null, Collections.emptyList(),
                        SVTestUtils.getLength(1000, 2001, type2), Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                // Should only cluster together if same type, except CNVs
                if ((type1 == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV && call2.isSimpleCNV()) ||
                        (type2 == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV && call1.isSimpleCNV())) {
                    Assert.assertTrue(engine.getLinkage().areClusterable(call1, call2));
                } else {
                    Assert.assertEquals(engine.getLinkage().areClusterable(call1, call2), type1 == type2);
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryStrands() {
        final List<Boolean> bools = Lists.newArrayList(Boolean.TRUE, Boolean.FALSE);
        for (final Boolean strand1A : bools) {
            for (final Boolean strand1B : bools) {
                final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, strand1A,
                        "chr1", 2001, strand1B, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                        null, Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                for (final Boolean strand2A : bools) {
                    for (final Boolean strand2B : bools) {
                        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, strand2A,
                                "chr1", 2001, strand2B, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                                null, Lists.newArrayList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                        // Should only cluster if strands match
                        Assert.assertEquals(engine.getLinkage().areClusterable(call1, call2), strand1A == strand2A && strand1B == strand2B);
                    }
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryContigs() {
        final List<String> contigs = Lists.newArrayList("chr1", "chrX");
        for (int i = 0; i < contigs.size(); i++) {
            final String contig1A = contigs.get(i);
            for (int j = i; j < contigs.size(); j++) {
                final String contig1B = contigs.get(j);
                final SVCallRecord call1 = new SVCallRecord("call1", contig1A, 1000, true,
                        contig1B, 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                        null, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                        Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                for (int k = 0; k < contigs.size(); k++) {
                    final String contig2A = contigs.get(k);
                    for (int m = k; m < contigs.size(); m++) {
                        final String contig2B = contigs.get(m);
                        final SVCallRecord call2 = new SVCallRecord("call2", contig2A, 1000, true,
                                contig2B, 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, Collections.emptyList(),
                                null, SVTestUtils.PESR_ONLY_ALGORITHM_LIST,
                                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                        // Should only cluster if contigs match
                        Assert.assertEquals(engine.getLinkage().areClusterable(call1, call2), contig1A.equals(contig2A) && contig1B.equals(contig2B));
                    }
                }
            }
        }
    }

    @Test
    public void testClusterTogetherVaryAlgorithms() {
        final List<List<String>> algorithmsList = Lists.newArrayList(
                Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Arrays.asList(GATKSVVCFConstants.DEPTH_ALGORITHM, SVTestUtils.PESR_ALGORITHM),
                SVTestUtils.PESR_ONLY_ALGORITHM_LIST
        );
        for (final List<String> algorithms1 : algorithmsList) {
            final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, true,
                    "chr1", 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                    1002, algorithms1, Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
            for (final List<String> algorithms2 : algorithmsList) {
                final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1000, true,
                        "chr1", 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                        1002, algorithms2, Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
                // All combinations should cluster
                Assert.assertTrue(engine.getLinkage().areClusterable(call1, call2));
            }
        }
    }

    @Test
    public void testClusterTogetherCNVs() {
        final SVCallRecord cnv1 = SVTestUtils.makeRecord("cnv1", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                Collections.emptyList());

        final SVCallRecord cnv2 = SVTestUtils.makeRecord("cnv2", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.CNV,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP),
                Collections.emptyList());

        final SVCallRecord del1 = SVTestUtils.makeRecord("del1", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DEL),
                Collections.emptyList());

        final SVCallRecord dup1 = SVTestUtils.makeRecord("dup1", "chr1", 1001, null,
                "chr1", 2001, null, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP,
                null, SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST,
                Lists.newArrayList(Allele.REF_N, Allele.SV_SIMPLE_DUP),
                Collections.emptyList());

        Assert.assertTrue(engine.getLinkage().areClusterable(cnv1, cnv2));
        Assert.assertTrue(engine.getLinkage().areClusterable(cnv1, del1));
        Assert.assertTrue(engine.getLinkage().areClusterable(cnv1, dup1));
        Assert.assertFalse(engine.getLinkage().areClusterable(del1, dup1));
    }

    @DataProvider(name = "testClusterTogetherIntervaledComplexData")
    public Object[][] testClusterTogetherIntervaledComplexData() {
        return new Object[][]{
                // exact match
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                            new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                            new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        true
                },
                // match within parameters
                {"chr1", 1100, "chr1", 1900,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                            new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1150, 1550)),
                            new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1650, 1800))
                        ),
                        true
                },
                // different contigs
                {"chr2", 1000, "chr2", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // different coordinates
                {"chr1", 1600, "chr1", 2400,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // different subtypes
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINVdel,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // different number of intervals
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500))
                        ),
                        false
                },
                // different cpx intervals
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 800, 1100)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
                // second cpx interval type is different
                {"chr1", 1000, "chr1", 2000,
                        GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, new SimpleInterval("chr1", 1100, 1500)),
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 1600, 1900))
                        ),
                        false
                },
        };
    }

    @Test(dataProvider= "testClusterTogetherIntervaledComplexData")
    public void testClusterTogetherIntervaledComplex(final String contigA, final int posA, final String contigB, final int posB,
                                                     final GATKSVVCFConstants.ComplexVariantSubtype subtype,
                                                     final List<SVCallRecord.ComplexEventInterval> cpxIntervals, final boolean expected) {
        final SVCallRecord cpx1 = new SVCallRecord("cpx1", "chr1", 1000, null,
                "chr1", 2000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.delINV,
                Arrays.asList(SVCallRecord.ComplexEventInterval.decode("DEL_chr1:1100-1500", SVTestUtils.hg38Dict), SVCallRecord.ComplexEventInterval.decode("INV_chr1:1600-1900", SVTestUtils.hg38Dict)),
                1000, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final Integer length2 = inferLength(contigA, posA, contigB, posB);
        final SVCallRecord cpx2 = new SVCallRecord("cpx2", contigA, posA, null,
                contigB, posB, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                subtype,
                cpxIntervals,
                length2, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertEquals(engine.getLinkage().areClusterable(cpx1, cpx2), expected);
    }

    @DataProvider(name = "testClusterTogetherInsertedComplexData")
    public Object[][] testClusterTogetherInsertedComplexData() {
        return new Object[][]{
                // exact match
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6000, 8000))
                        ),
                        true
                },
                // close match
                {"chr1", 1010, "chr1", 1010,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6200, 7990))
                        ),
                        true
                },
                // not match by coordinates
                {"chr1", 2000, "chr1", 3000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6200, 7990))
                        ),
                        false
                },
                // not match by subtype
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6000, 8000))
                        ),
                        false
                },
                // not match by cpx interval
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 9000, 11000))
                        ),
                        false
                },
                // different cpx interval type
                {"chr1", 1000, "chr1", 1000,
                        GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                        Arrays.asList(
                                new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.INV, new SimpleInterval("chr1", 6000, 8000))
                        ),
                        false
                },
        };
    }

    @Test(dataProvider= "testClusterTogetherInsertedComplexData")
    public void testClusterTogetherInsertedComplex(final String contigA, final int posA, final String contigB, final int posB,
                                                   final GATKSVVCFConstants.ComplexVariantSubtype subtype,
                                                   final List<SVCallRecord.ComplexEventInterval> cpxIntervals, final boolean expected) {
        final SVCallRecord cpx1 = new SVCallRecord("cpx1", "chr1", 1000, null,
                "chr1", 1000, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                GATKSVVCFConstants.ComplexVariantSubtype.dDUP,
                Arrays.asList(new SVCallRecord.ComplexEventInterval(GATKSVVCFConstants.StructuralVariantAnnotationType.DUP, new SimpleInterval("chr1", 6000, 8000))),
                2000, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final Integer length2 = cpxIntervals.get(0).getInterval().size();
        final SVCallRecord cpx2 = new SVCallRecord("cpx2", contigA, posA, null,
                contigB, posB, null,
                GATKSVVCFConstants.StructuralVariantAnnotationType.CPX,
                subtype,
                cpxIntervals,
                length2, Collections.singletonList(SVTestUtils.PESR_ALGORITHM),
                Lists.newArrayList(Allele.REF_N, SVTestUtils.CPX_ALLELE),
                Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        Assert.assertEquals(engine.getLinkage().areClusterable(cpx1, cpx2), expected);
    }

    @Test
    public void testClusterTogetherVaryParameters() {
        final SVClusterEngine testEngine1 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", 1000, true,
                "chr1", 2001, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1002, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2 = new SVCallRecord("call2", "chr1", 1100, true,
                "chr1", 2101, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                1002, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM), Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        // Cluster with default parameters
        Assert.assertTrue(testEngine1.getLinkage().areClusterable(call1, call2));
        final ClusteringParameters exactMatchParameters = ClusteringParameters.createDepthParameters(1.0, 0, 0, 1.0);
        final CanonicalSVLinkage<SVCallRecord> exactMatchLinkage = SVTestUtils.getNewDefaultLinkage();
        exactMatchLinkage.setDepthOnlyParams(exactMatchParameters);
        // Do not cluster requiring exact overlap
        Assert.assertFalse(exactMatchLinkage.areClusterable(call1, call2));
    }


    @DataProvider(name = "testAddVaryPositionsProvider")
    public Object[][] testAddVaryPositionsProvider() {
        return new Object[][]{
                {10000, 20001, 12000, 22001, 14000, 24001, SINGLE_LINKAGE, 1},
                {10000, 20001, 12000, 22001, 14001, 24002, SINGLE_LINKAGE, 2},
                {10000, 20001, 12001, 22002, 14002, 24003, SINGLE_LINKAGE, 3},
                {10000, 20001, 11000, 21001, 12000, 22001, MAX_CLIQUE, 1},
                {10000, 20001, 12000, 22001, 14000, 24001, MAX_CLIQUE, 2},
                {10000, 20001, 12001, 22002, 14002, 24003, MAX_CLIQUE, 3}
        };
    }

    @Test(dataProvider= "testAddVaryPositionsProvider")
    public void testAddVaryPositions(final int positionA1, final int positionB1,
                                     final int positionA2, final int positionB2,
                                     final int positionA3, final int positionB3,
                                     final SVClusterEngine.CLUSTERING_TYPE type,
                                     final int result) {
        final SVClusterEngine engine;
        if (type == SINGLE_LINKAGE) {
            engine = SVTestUtils.getNewDefaultSingleLinkageEngine();
        } else if (type == MAX_CLIQUE) {
            engine = SVTestUtils.getNewDefaultMaxCliqueEngine();
        } else {
            throw new TestException("Unimplemented clustering type " + type.name());
        }
        final SVCallRecord call1 = new SVCallRecord("call1", "chr1", positionA1, true,
                "chr1", positionB1, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                positionB1 - positionA1 + 1, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call2 = new SVCallRecord("call1", "chr1", positionA2, true,
                "chr1", positionB2, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                positionB2 - positionA2 + 1, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final SVCallRecord call3 = new SVCallRecord("call1", "chr1", positionA3, true,
                "chr1", positionB3, false, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, Collections.emptyList(),
                positionB3 - positionA3 + 1, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        engine.add(call1);
        engine.add(call2);
        engine.add(call3);
        Assert.assertEquals(engine.forceFlush().size(), result);
    }

    @Test
    public void testAdd() {
        //single-sample merge case, ignoring sample sets
        final SVClusterEngine temp1 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        Assert.assertTrue(temp1.isEmpty());
        temp1.add(SVTestUtils.call1);
        Assert.assertFalse(temp1.isEmpty());
        //force new cluster by adding a non-overlapping event
        temp1.add(SVTestUtils.call3);
        final List<SVCallRecord> output1 = temp1.forceFlush(); //flushes all clusters
        Assert.assertTrue(temp1.isEmpty());
        Assert.assertEquals(output1.size(), 2);
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call1, output1.get(0));
        SVTestUtils.assertEqualsExceptMembershipAndGT(SVTestUtils.call3, output1.get(1));

        final SVClusterEngine temp2 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        temp2.add(SVTestUtils.call1);
        temp2.add(SVTestUtils.overlapsCall1);
        //force new cluster by adding a call on another contig
        temp2.add(SVTestUtils.call4_chr10);
        final List<SVCallRecord> output2 = temp2.forceFlush();
        Assert.assertEquals(output2.size(), 2);
        //median of two items ends up being the second item here
        Assert.assertEquals(output2.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output2.get(0).getPositionB(), SVTestUtils.call1.getPositionB());
        SVTestUtils.assertEqualsExceptMembershipAndGT(output2.get(1), SVTestUtils.call4_chr10);

        //checking insensitivity to sample set overlap
        final SVClusterEngine temp3 = SVTestUtils.getNewDefaultSingleLinkageEngine();
        temp3.add(SVTestUtils.call1);
        temp3.add(SVTestUtils.sameBoundsSampleMismatch);
        final List<SVCallRecord> output3 = temp3.forceFlush();
        Assert.assertEquals(output3.size(), 1);
        Assert.assertEquals(output3.get(0).getPositionA(), SVTestUtils.call1.getPositionA());
        Assert.assertEquals(output3.get(0).getPositionB(), SVTestUtils.call1.getPositionB());
        Assert.assertEquals(output3.get(0).getPositionA(), SVTestUtils.sameBoundsSampleMismatch.getPositionA());
        Assert.assertEquals(output3.get(0).getPositionB(), SVTestUtils.sameBoundsSampleMismatch.getPositionB());
    }

    @Test
    public void testAddMaxCliqueLarge() {
        final int numRecords = 100;
        final SVClusterEngine engine = SVTestUtils.getNewDefaultMaxCliqueEngine();
        final int length = 5000;
        for (int i = 0; i < numRecords; i++) {
            final int start = 1000 + 10 * i;
            final int end = start + length - 1;
            engine.add(SVTestUtils.newPESRCallRecordWithIntervalAndType(start, end, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL));
        }
        final List<SVCallRecord> result = engine.forceFlush();
        Assert.assertEquals(result.size(), 50);
        for (final SVCallRecord resultRecord : result) {
            Assert.assertTrue(resultRecord.getAttributes().containsKey(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY));
            final int numMembers = VariantContextGetters.attributeToList(resultRecord.getAttributes().get(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY)).size();
            Assert.assertEquals(numMembers, 51);
        }
    }

    @DataProvider(name = "testGetCarrierSamplesBiallelicData")
    public Object[][] testGetCarrierSamplesBiallelicData() {
        return new Object[][]{
                // Empty
                {0, null, null, new int[]{}, new String[]{}},
                {0, null, Allele.SV_SIMPLE_DEL, new int[]{}, new String[]{}},
                {0, Allele.REF_N, null, new int[]{}, new String[]{}},
                {0, Allele.NO_CALL, null, new int[]{}, new String[]{}},
                {0, Allele.NO_CALL, Allele.SV_SIMPLE_DEL, new int[]{}, new String[]{}},
                {0, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{}, new String[]{}},
                // Ploidy 0
                {0, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{0}, new String[]{}},
                // Ploidy 1, no carrier
                {1, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{1, 1}, new String[]{}},
                // Ploidy 1, with carrier
                {1, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{1, 0}, new String[]{"1"}},
                // Ploidy 2, no carrier
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 2, 2}, new String[]{}},
                // Ploidy 2, with carriers
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 1, 1}, new String[]{"1", "2"}},
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 1, 0}, new String[]{"1", "2"}},
                {2, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{2, 1, 3}, new String[]{"1"}}, // third is not consistent with DEL
                // DUP
                {2, Allele.REF_N, Allele.SV_SIMPLE_DUP, new int[]{2, 3, 2}, new String[]{"1"}},
                {2, Allele.REF_N, Allele.SV_SIMPLE_DUP, new int[]{2, 1, 3}, new String[]{"2"}}, // second is not consistent with DUP
                // Ploidy 3
                {3, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{3, 3, 3}, new String[]{}},
                {3, Allele.REF_N, Allele.SV_SIMPLE_DEL, new int[]{0, 1, 2, 3}, new String[]{"0", "1", "2"}}
        };
    }

    @Test(dataProvider= "testGetCarrierSamplesBiallelicData")
    public void testGetCarrierSamplesBiallelic(final int ploidy, final Allele refAllele, final Allele altAllele,
                                               final int[] copyNumbers, final String[] expectedSamples) {
        final Set<String> expectedResult = Stream.of(expectedSamples).collect(Collectors.toSet());
        final List<Allele> alleles = new ArrayList<>(2);
        if (refAllele != null) {
            alleles.add(refAllele);
        }
        final GATKSVVCFConstants.StructuralVariantAnnotationType svtype;
        if (altAllele != null) {
            if (altAllele.equals(Allele.SV_SIMPLE_DEL)) {
                svtype = GATKSVVCFConstants.StructuralVariantAnnotationType.DEL;
            } else if (altAllele.equals(Allele.SV_SIMPLE_DUP)) {
                svtype = GATKSVVCFConstants.StructuralVariantAnnotationType.DUP;
            } else {
                throw new TestException("Invalid alt allele: " + altAllele.getDisplayString());
            }
            alleles.add(altAllele);
        } else {
            svtype = GATKSVVCFConstants.StructuralVariantAnnotationType.DEL; // default
        }

        // Create genotypes with copy number attribute (and no GT)
        final List<Genotype> genotypesWithCopyNumber = IntStream.range(0, copyNumbers.length)
                .mapToObj(i -> new GenotypeBuilder(String.valueOf(i))
                        .attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumbers[i])
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy)
                        .alleles(SVTestUtils.buildHomAlleleListWithPloidy(Allele.NO_CALL, ploidy))
                        .make())
                .collect(Collectors.toList());
        final SVCallRecord recordWithCopyNumber = new SVCallRecord("", "chr1", 1000, SVTestUtils.getValidTestStrandA(svtype),
                "chr1", 1999, SVTestUtils.getValidTestStrandB(svtype), svtype, null, Collections.emptyList(),
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                alleles, GenotypesContext.copy(genotypesWithCopyNumber), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final Set<String> resultWithCopyNumber =  recordWithCopyNumber.getCarrierSampleSet();

        Assert.assertEquals(resultWithCopyNumber, expectedResult);

        // Create genotypes with GT (and no copy number attribute)
        final List<Genotype> genotypesWithGenotype = IntStream.range(0, copyNumbers.length)
                .mapToObj(i -> new GenotypeBuilder(String.valueOf(i))
                        .alleles(SVTestUtils.buildBiallelicListWithPloidy(altAllele, refAllele, copyNumbers[i], ploidy))
                        .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ploidy)
                        .make())
                .collect(Collectors.toList());
        final SVCallRecord recordWithGenotype = new SVCallRecord("", "chr1", 1000, SVTestUtils.getValidTestStrandA(svtype),
                "chr1", 1999, SVTestUtils.getValidTestStrandB(svtype), svtype, null, Collections.emptyList(),
                1000, Collections.singletonList(GATKSVVCFConstants.DEPTH_ALGORITHM),
                alleles, GenotypesContext.copy(genotypesWithGenotype), Collections.emptyMap(), Collections.emptySet(), null, SVTestUtils.hg38Dict);
        final Set<String> resultWithGenotype = recordWithGenotype.getCarrierSampleSet();

        Assert.assertEquals(resultWithGenotype, expectedResult);
    }

    @Test
    public void testLargeRandom() {
        final Random rand = new Random(42);
        final List<SVCallRecord> records = new ArrayList<>(1000);
        for (int i = 0; i < 1000; i++) {
            final int pos1 = rand.nextInt(9000) + 1;
            final int pos2 = rand.nextInt(10000) + 1;
            records.add(SVTestUtils.newPESRCallRecordWithIntervalAndType(Math.min(pos1, pos2), Math.max(pos1, pos2), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL));
        }
        final SVClusterEngine engine = SVTestUtils.getNewDefaultMaxCliqueEngine();
        records.stream().sorted(SVCallRecordUtils.getCallComparator(SVTestUtils.hg38Dict)).forEach(engine::add);
        final List<SVCallRecord> output = engine.forceFlush();
        Assert.assertEquals(output.size(), 2926);
    }
}