package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.StructuralVariantAnnotationType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Checks that {@link SVClusterEngine} emits exactly the clusters of the frozen pre-rework {@link LegacySVClusterEngine}:
 * same clusters, same member order within each cluster, same emission order, and same per-call flush boundaries,
 * on randomized dense record streams. Also checks that {@link CanonicalSVLinkage#getMaxLinkableStartingPosition} is
 * never violated, since the engine relies on it to skip pairwise tests.
 */
public class SVClusterEngineEquivalenceTest {

    private static final List<String> CONTIGS = Arrays.asList("chr1", "chr2", "chr3");
    private static final List<String> DEPTH = SVTestUtils.DEPTH_ONLY_ALGORITHM_LIST;
    private static final List<String> PESR = SVTestUtils.PESR_ONLY_ALGORITHM_LIST;
    private static final List<String> BOTH = Arrays.asList(SVTestUtils.PESR_ALGORITHM, GATKSVVCFConstants.DEPTH_ALGORITHM);

    /**
     * Random coordinate-sorted records, dense enough to form large chained clusters: a mix of small and large
     * DEL/DUP/CNV/INV/INS/BND/CPX with depth-only, PE/SR, and mixed evidence, tiny lengths that exercise rounding in
     * the overlap bound, and stacks of identical calls.
     */
    private static List<SVCallRecord> randomRecords(final long seed, final int numRecords, final int regionSize) {
        final Random random = new Random(seed);
        final List<SVCallRecord> records = new ArrayList<>(numRecords);
        for (int i = 0; i < numRecords; i++) {
            final String contig = CONTIGS.get(random.nextInt(CONTIGS.size()));
            final int start = 1 + random.nextInt(regionSize);
            final int stackSize = random.nextInt(20) == 0 ? 2 + random.nextInt(10) : 1;
            final SVCallRecord record = randomRecord(random, "r" + i, contig, start);
            for (int j = 0; j < stackSize && records.size() < numRecords; j++) {
                records.add(j == 0 ? record : copyWithId(record, "r" + i + "_" + j));
            }
        }
        final Map<String, Integer> contigOrder = new HashMap<>();
        for (int i = 0; i < CONTIGS.size(); i++) {
            contigOrder.put(CONTIGS.get(i), i);
        }
        // Stable sort: ties keep generation order, as an input merge would
        records.sort(Comparator.comparing((SVCallRecord r) -> contigOrder.get(r.getContigA()))
                .thenComparingInt(SVCallRecord::getPositionA));
        return records;
    }

    private static SVCallRecord randomRecord(final Random random, final String id, final String contig, final int start) {
        final int lengthClass = random.nextInt(10);
        final int length;
        if (lengthClass < 2) {
            length = 1 + random.nextInt(10);
        } else if (lengthClass < 8) {
            length = 50 + random.nextInt(5000);
        } else {
            length = 50_000 + random.nextInt(3_000_000);
        }
        final List<String> algorithms;
        final int algorithmClass = random.nextInt(3);
        algorithms = algorithmClass == 0 ? DEPTH : (algorithmClass == 1 ? PESR : BOTH);
        final StructuralVariantAnnotationType[] types = {StructuralVariantAnnotationType.DEL,
                StructuralVariantAnnotationType.DUP, StructuralVariantAnnotationType.CNV, StructuralVariantAnnotationType.INV,
                StructuralVariantAnnotationType.INS, StructuralVariantAnnotationType.BND, StructuralVariantAnnotationType.CPX};
        final StructuralVariantAnnotationType type = types[random.nextInt(types.length)];
        switch (type) {
            case DEL:
                return newRecord(id, contig, start, true, contig, start + length - 1, false, type, null, length, algorithms);
            case DUP:
                return newRecord(id, contig, start, false, contig, start + length - 1, true, type, null, length, algorithms);
            case CNV:
                return newRecord(id, contig, start, null, contig, start + length - 1, null, type, null, length, DEPTH);
            case INV: {
                final boolean strand = random.nextBoolean();
                return newRecord(id, contig, start, strand, contig, start + length - 1, strand, type, null, length, PESR);
            }
            case INS:
                return newRecord(id, contig, start, true, contig, start, false, type, null,
                        random.nextInt(4) == 0 ? null : length, PESR);
            case BND: {
                // Contig B may not precede contig A
                final int contigIndex = CONTIGS.indexOf(contig);
                final String contigB = random.nextBoolean() ? contig
                        : CONTIGS.get(contigIndex + random.nextInt(CONTIGS.size() - contigIndex));
                final int positionB = contigB.equals(contig) ? start + random.nextInt(length) : 1 + random.nextInt(100_000);
                return newRecord(id, contig, start, random.nextBoolean(), contigB, positionB, random.nextBoolean(),
                        type, null, null, PESR);
            }
            default: {
                // CPX with and without complex intervals; empty interval lists link at any distance
                final List<SVCallRecord.ComplexEventInterval> intervals = random.nextBoolean() ? Collections.emptyList()
                        : Collections.singletonList(new SVCallRecord.ComplexEventInterval(StructuralVariantAnnotationType.DUP,
                        new SimpleInterval(contig, start, start + length - 1)));
                return newRecord(id, contig, start, null, contig, start + length - 1, null, type, intervals, length, PESR);
            }
        }
    }

    private static SVCallRecord newRecord(final String id, final String contigA, final int positionA, final Boolean strandA,
                                          final String contigB, final int positionB, final Boolean strandB,
                                          final StructuralVariantAnnotationType type,
                                          final List<SVCallRecord.ComplexEventInterval> cpxIntervals,
                                          final Integer length, final List<String> algorithms) {
        return new SVCallRecord(id, contigA, positionA, strandA, contigB, positionB, strandB, type,
                type == StructuralVariantAnnotationType.CPX ? GATKSVVCFConstants.ComplexVariantSubtype.dDUP : null,
                cpxIntervals == null ? Collections.emptyList() : cpxIntervals, length, Collections.emptyList(),
                algorithms, Collections.emptyList(), Collections.emptyList(), Collections.emptyMap(),
                Collections.emptySet(), null, SVTestUtils.hg38Dict);
    }

    private static SVCallRecord copyWithId(final SVCallRecord r, final String id) {
        return new SVCallRecord(id, r.getContigA(), r.getPositionA(), r.getStrandA(), r.getContigB(), r.getPositionB(),
                r.getStrandB(), r.getType(), r.getComplexSubtype(), r.getComplexEventIntervals(), r.getLength(),
                r.getEvidence(), r.getAlgorithms(), r.getAlleles(), r.getGenotypes(), r.getAttributes(), r.getFilters(),
                r.getLog10PError(), SVTestUtils.hg38Dict);
    }

    private static CanonicalSVLinkage<SVCallRecord> linkage(final double overlap, final int pesrWindow,
                                                           final int depthWindow, final boolean overlapAndProximity,
                                                           final boolean clusterDelWithDup) {
        final CanonicalSVLinkage<SVCallRecord> linkage = new CanonicalSVLinkage<>(SVTestUtils.hg38Dict, clusterDelWithDup);
        linkage.setDepthOnlyParams(new ClusteringParameters(overlap, 0, depthWindow, 0, overlapAndProximity,
                (a, b) -> a.isDepthOnly() && b.isDepthOnly()));
        linkage.setMixedParams(new ClusteringParameters(overlap, 0, 1000, 0, overlapAndProximity,
                (a, b) -> a.isDepthOnly() != b.isDepthOnly()));
        linkage.setEvidenceParams(new ClusteringParameters(Math.min(overlap, 0.5), 0.5, pesrWindow, 0, overlapAndProximity,
                (a, b) -> !a.isDepthOnly() && !b.isDepthOnly()));
        return linkage;
    }

    @DataProvider(name = "linkages")
    public Object[][] linkages() {
        return new Object[][]{
                // overlap, PE/SR window, depth window, overlap-and-proximity, DEL/DUP clustering
                {0.8, 500, 10_000_000, true, false},   // tool defaults
                {0.8, 500, 10_000_000, true, true},
                {0.5, 5000, 100_000, true, false},
                {0.1, 500, 10_000_000, true, false},
                {0.0, 500, 1000, true, false},          // zero overlap threshold: proximity alone decides
                {0.8, 500, 1000, false, false},         // overlap-or-proximity mode
        };
    }

    @Test(dataProvider = "linkages")
    public void testLinkageBoundIsNeverViolated(final double overlap, final int pesrWindow, final int depthWindow,
                                                final boolean overlapAndProximity, final boolean clusterDelWithDup) {
        final CanonicalSVLinkage<SVCallRecord> linkage = linkage(overlap, pesrWindow, depthWindow, overlapAndProximity, clusterDelWithDup);
        final List<SVCallRecord> records = randomRecords(17, 3000, 200_000);
        for (int i = 0; i < records.size(); i++) {
            final SVCallRecord a = records.get(i);
            final int bound = linkage.getMaxLinkableStartingPosition(a);
            for (int j = i + 1; j < records.size(); j++) {
                final SVCallRecord b = records.get(j);
                if (!b.getContigA().equals(a.getContigA())) {
                    break;
                }
                if (b.getPositionA() > bound) {
                    Assert.assertFalse(linkage.areClusterable(a, b).getResult(), a.getId() + " links " + b.getId() + " past bound " + bound);
                    Assert.assertFalse(linkage.areClusterable(b, a).getResult(), b.getId() + " links " + a.getId() + " past bound " + bound);
                }
            }
        }
    }

    @Test(dataProvider = "linkages")
    public void testMatchesLegacyEngine(final double overlap, final int pesrWindow, final int depthWindow,
                                        final boolean overlapAndProximity, final boolean clusterDelWithDup) {
        for (final SVClusterEngine.CLUSTERING_TYPE type : SVClusterEngine.CLUSTERING_TYPE.values()) {
            for (long seed = 0; seed < 4; seed++) {
                final CanonicalSVLinkage<SVCallRecord> linkage = linkage(overlap, pesrWindow, depthWindow, overlapAndProximity, clusterDelWithDup);
                // Max-clique enumerates cliques, so keep its streams sparser
                final int regionSize = type == SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE ? 2_000_000 : 300_000;
                assertSameClusters(type, linkage, randomRecords(seed, 2000, regionSize));
            }
        }
    }

    /**
     * A linkage that does not override {@link SVClusterLinkage#getMaxLinkableStartingPosition} disables pruning; the
     * remaining bookkeeping changes must still reproduce the legacy output.
     */
    @Test
    public void testMatchesLegacyEngineWithoutPruning() {
        final CanonicalSVLinkage<SVCallRecord> canonical = linkage(0.8, 500, 10_000_000, true, false);
        final SVClusterLinkage<SVCallRecord> unpruned = new SVClusterLinkage<SVCallRecord>() {
            @Override
            public LinkageResult areClusterable(final SVCallRecord a, final SVCallRecord b) {
                return canonical.areClusterable(a, b);
            }

            @Override
            public int getMaxClusterableStartingPosition(final SVCallRecord item) {
                return canonical.getMaxClusterableStartingPosition(item);
            }
        };
        for (final SVClusterEngine.CLUSTERING_TYPE type : SVClusterEngine.CLUSTERING_TYPE.values()) {
            assertSameClusters(type, unpruned, randomRecords(5, 2000, 300_000));
        }
    }

    /**
     * The same large event called by depth and by PE/SR in many batches: depth-only members keep the chained cluster
     * open for 0.2 x length, while PE/SR members stop being linkable after the 1 kb mixed window. Those PE/SR members
     * should no longer be tested against every new item.
     */
    @Test
    public void testPrunesPairwiseTests() {
        final Random random = new Random(11);
        final List<SVCallRecord> records = new ArrayList<>();
        for (int i = 0; i < 3000; i++) {
            final int start = 1 + random.nextInt(2_000_000);
            final int length = 2_400_000 + random.nextInt(200_000);
            records.add(newRecord("dup" + i, "chr1", start, false, "chr1", start + length - 1, true,
                    StructuralVariantAnnotationType.DUP, null, length, random.nextBoolean() ? DEPTH : PESR));
        }
        records.sort(Comparator.comparingInt(SVCallRecord::getPositionA));
        final long[] tests = countPairwiseTests(records);
        Assert.assertTrue(tests[1] * 2 < tests[0], "pairwise tests: " + tests[1] + " vs legacy " + tests[0]);
    }

    /**
     * On a generic stream the engines must agree, and pruning must never add pairwise tests.
     */
    @Test
    public void testNeverAddsPairwiseTests() {
        final long[] tests = countPairwiseTests(randomRecords(3, 4000, 300_000));
        Assert.assertTrue(tests[1] <= tests[0], "pairwise tests: " + tests[1] + " vs legacy " + tests[0]);
    }

    /**
     * Returns {legacy, reworked} pairwise test counts for single-linkage over the given records, after checking that
     * both engines emit the same clusters.
     */
    private static long[] countPairwiseTests(final List<SVCallRecord> records) {
        final CanonicalSVLinkage<SVCallRecord> canonical = linkage(0.8, 500, 10_000_000, true, false);
        final AtomicLong count = new AtomicLong();
        final SVClusterLinkage<SVCallRecord> counting = new SVClusterLinkage<SVCallRecord>() {
            @Override
            public LinkageResult areClusterable(final SVCallRecord a, final SVCallRecord b) {
                count.incrementAndGet();
                return canonical.areClusterable(a, b);
            }

            @Override
            public int getMaxClusterableStartingPosition(final SVCallRecord item) {
                return canonical.getMaxClusterableStartingPosition(item);
            }

            @Override
            public int getMaxLinkableStartingPosition(final SVCallRecord item) {
                return canonical.getMaxLinkableStartingPosition(item);
            }
        };
        final List<List<String>> legacyClusters = new ArrayList<>();
        final List<List<String>> clusters = new ArrayList<>();
        final LegacySVClusterEngine legacy = new LegacySVClusterEngine(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE,
                recordingCollapser(legacyClusters, LegacySVClusterEngine.OutputCluster::getItems), counting, SVTestUtils.hg38Dict);
        records.forEach(legacy::addAndFlush);
        legacy.flush();
        final long legacyTests = count.getAndSet(0);
        final SVClusterEngine engine = new SVClusterEngine(SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE,
                recordingCollapser(clusters, SVClusterEngine.OutputCluster::getItems), counting, SVTestUtils.hg38Dict);
        records.forEach(engine::addAndFlush);
        engine.flush();
        Assert.assertEquals(clusters, legacyClusters);
        return new long[]{legacyTests, count.get()};
    }

    private static void assertSameClusters(final SVClusterEngine.CLUSTERING_TYPE type,
                                           final SVClusterLinkage<SVCallRecord> linkage,
                                           final List<SVCallRecord> records) {
        final List<List<String>> legacyClusters = new ArrayList<>();
        final List<List<String>> clusters = new ArrayList<>();
        final LegacySVClusterEngine legacy = new LegacySVClusterEngine(type,
                recordingCollapser(legacyClusters, LegacySVClusterEngine.OutputCluster::getItems), linkage, SVTestUtils.hg38Dict);
        final SVClusterEngine engine = new SVClusterEngine(type,
                recordingCollapser(clusters, SVClusterEngine.OutputCluster::getItems), linkage, SVTestUtils.hg38Dict);
        int numMultiItemClusters = 0;
        for (final SVCallRecord record : records) {
            final int legacyEmitted = legacy.addAndFlush(record).size();
            final int emitted = engine.addAndFlush(record).size();
            Assert.assertEquals(emitted, legacyEmitted, "clusters emitted when adding " + record.getId());
        }
        Assert.assertEquals(engine.flush().size(), legacy.flush().size());
        Assert.assertEquals(clusters, legacyClusters);
        for (final List<String> cluster : clusters) {
            if (cluster.size() > 1) {
                numMultiItemClusters++;
            }
        }
        // Guard against vacuous streams
        Assert.assertTrue(numMultiItemClusters > 10, "only " + numMultiItemClusters + " multi-item clusters");
    }

    private static <C> Function<C, SVCallRecord> recordingCollapser(final List<List<String>> clusters,
                                                                    final Function<C, List<SVCallRecord>> getItems) {
        return cluster -> {
            final List<SVCallRecord> items = getItems.apply(cluster);
            clusters.add(items.stream().map(SVCallRecord::getId).collect(Collectors.toList()));
            return items.get(0);
        };
    }
}
