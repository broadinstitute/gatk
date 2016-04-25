package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class SharedVertexSequenceSplitterUnitTest extends BaseTest {
    private final static boolean PRINT_GRAPHS = false;

    @DataProvider(name = "PrefixSuffixData")
    public Object[][] makePrefixSuffixData() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{Arrays.asList("A", "C"), 0, 0});
        tests.add(new Object[]{Arrays.asList("C", "C"), 1, 0});
        tests.add(new Object[]{Arrays.asList("ACT", "AGT"), 1, 1});
        tests.add(new Object[]{Arrays.asList("ACCT", "AGT"), 1, 1});
        tests.add(new Object[]{Arrays.asList("ACT", "ACT"), 3, 0});
        tests.add(new Object[]{Arrays.asList("ACTA", "ACT"), 3, 0});
        tests.add(new Object[]{Arrays.asList("ACTA", "ACTG"), 3, 0});
        tests.add(new Object[]{Arrays.asList("ACTA", "ACTGA"), 3, 1});
        tests.add(new Object[]{Arrays.asList("GCTGA", "ACTGA"), 0, 4});

        tests.add(new Object[]{Arrays.asList("A", "C", "A"), 0, 0});
        tests.add(new Object[]{Arrays.asList("A", "A", "A"), 1, 0});
        tests.add(new Object[]{Arrays.asList("A", "AA", "A"), 1, 0});
        tests.add(new Object[]{Arrays.asList("A", "ACA", "A"), 1, 0});
        tests.add(new Object[]{Arrays.asList("ACT", "ACAT", "ACT"), 2, 1});
        tests.add(new Object[]{Arrays.asList("ACT", "ACAT", "ACGT"), 2, 1});
        tests.add(new Object[]{Arrays.asList("AAAT", "AAA", "CAAA"), 0, 0});
        tests.add(new Object[]{Arrays.asList("AACTTT", "AAGTTT", "AAGCTTT"), 2, 3});
        tests.add(new Object[]{Arrays.asList("AAA", "AAA", "CAAA"), 0, 3});
        tests.add(new Object[]{Arrays.asList("AAA", "AAA", "AAA"), 3, 0});

        tests.add(new Object[]{Arrays.asList("AC", "ACA", "AC"), 2, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PrefixSuffixData")
    public void testPrefixSuffix(final List<String> strings, int expectedPrefixLen, int expectedSuffixLen) {
        final List<byte[]> bytes = new ArrayList<>();
        int min = Integer.MAX_VALUE;
        for ( final String s : strings ) {
            bytes.add(s.getBytes());
            min = Math.min(min, s.length());
        }

        final int actualPrefixLen = GraphUtils.commonMaximumPrefixLength(bytes);
        Assert.assertEquals(actualPrefixLen, expectedPrefixLen, "Failed prefix test");

        final int actualSuffixLen = GraphUtils.commonMaximumSuffixLength(bytes, min - actualPrefixLen);
        Assert.assertEquals(actualSuffixLen, expectedSuffixLen, "Failed suffix test");
    }

    @Test(dataProvider = "PrefixSuffixData")
    public void testPrefixSuffixVertices(final List<String> strings, int expectedPrefixLen, int expectedSuffixLen) {
        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : strings ) {
            v.add(new SeqVertex(s));
        }

        final String expectedPrefix = strings.get(0).substring(0, expectedPrefixLen);
        final String expectedSuffix = strings.get(0).substring(strings.get(0).length() - expectedSuffixLen);

        final Pair<SeqVertex, SeqVertex> result = SharedVertexSequenceSplitter.commonPrefixAndSuffixOfVertices(v);
        Assert.assertEquals(result.getLeft().getSequenceString(), expectedPrefix, "Failed suffix test");
        Assert.assertEquals(result.getRight().getSequenceString(), expectedSuffix, "Failed suffix test");

        Assert.assertEquals(result.getLeft().isEmpty(), expectedPrefix.isEmpty());
        Assert.assertEquals(result.getRight().isEmpty(), expectedSuffix.isEmpty());
    }

    @Test(dataProvider = "PrefixSuffixData")
    public void testSplitter(final List<String> strings, int expectedPrefixLen, int expectedSuffixLen) {
        final SeqGraph graph = new SeqGraph(11);

        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : strings ) {
            v.add(new SeqVertex(s));
        }

        graph.addVertices(v.toArray(new SeqVertex[v.size()]));

        final String expectedPrefix = strings.get(0).substring(0, expectedPrefixLen);
        final String expectedSuffix = strings.get(0).substring(strings.get(0).length() - expectedSuffixLen);

        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(graph, v);
        splitter.split();

        Assert.assertEquals(splitter.getPrefixV().getSequenceString(), expectedPrefix);
        Assert.assertEquals(splitter.getSuffixV().getSequenceString(), expectedSuffix);

        Assert.assertTrue(splitter.getSplitGraph().outDegreeOf(splitter.getPrefixV()) <= strings.size());
        Assert.assertEquals(splitter.getSplitGraph().inDegreeOf(splitter.getPrefixV()), 0);

        Assert.assertTrue(splitter.getSplitGraph().inDegreeOf(splitter.getSuffixV()) <= strings.size());
        Assert.assertEquals(splitter.getSplitGraph().outDegreeOf(splitter.getSuffixV()), 0);

        for ( final SeqVertex mid : splitter.getNewMiddles()) {
            Assert.assertNotNull(splitter.getSplitGraph().getEdge(splitter.getPrefixV(), mid));
            Assert.assertNotNull(splitter.getSplitGraph().getEdge(mid, splitter.getSuffixV()));
        }
    }

    @DataProvider(name = "CompleteCycleData")
    public Object[][] makeCompleteCycleData() {
        List<Object[]> tests = new ArrayList<>();

        for ( final boolean hasTop : Arrays.asList(true, false) ) {
            for ( final boolean hasBot : Arrays.asList(true, false) ) {
                if ( ! hasTop && ! hasBot ) continue;
                tests.add(new Object[]{Arrays.asList("A", "A"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "AC"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "CA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("AC", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("AT", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("ATA", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("ATAA", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("ATAACA", "ACA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CCCAAA", "AAA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CCCAAAAAA", "AAA"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CCCAAAAAA", "CCCAAA"), hasTop, hasBot});

                tests.add(new Object[]{Arrays.asList("A", "A", "A"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "A", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("A", "C", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("AC", "C", "C"), hasTop, hasBot});
                tests.add(new Object[]{Arrays.asList("CA", "C", "C"), hasTop, hasBot});
                // all merged
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "AGA"), hasTop, hasBot});
                // prefix and suffix
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "ACA"), hasTop, hasBot});
                // 2 -> prefix, leave C
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "AGAC"), hasTop, hasBot});
                // 2 -> prefix, leave CCC
                tests.add(new Object[]{Arrays.asList("AGA", "AGA", "AGACCC"), hasTop, hasBot});
                // 2 -> suffix, leave A/T
                tests.add(new Object[]{Arrays.asList("TAGA", "TAGA", "AAGA"), hasTop, hasBot});
                // 2 -> suffix, leave T, delete 1
                tests.add(new Object[]{Arrays.asList("TAGA", "TAGA", "AGA"), hasTop, hasBot});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CompleteCycleData")
    public void testSplitterCompleteCycle(final List<String> strings, final boolean hasTop, final boolean hasBot) {
        final SeqGraph graph = new SeqGraph(11);

        int edgeWeight = 1;
        final SeqVertex top = hasTop ? new SeqVertex("AAAAAAAA") : null;
        final SeqVertex bot = hasBot ? new SeqVertex("GGGGGGGG") : null;
        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : strings ) {
            v.add(new SeqVertex(s));
        }
        graph.addVertices(v.toArray(new SeqVertex[v.size()]));
        final SeqVertex first = v.get(0);

        if ( hasTop ) {
            graph.addVertex(top);
            for ( final SeqVertex vi : v )
                graph.addEdge(top, vi, new BaseEdge(vi == first, edgeWeight++));
        }

        if ( hasBot ) {
            graph.addVertex(bot);
            for ( final SeqVertex vi : v )
                graph.addEdge(vi, bot, new BaseEdge(vi == first, edgeWeight++));
        }

        final Set<String> haplotypes = new HashSet<>();
        final KBestHaplotypeFinder originalPaths = new KBestHaplotypeFinder(graph.clone(),graph.getSources(),graph.getSinks());
        for ( final KBestHaplotype path : originalPaths )
            haplotypes.add(new String(path.bases()));

        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(graph, v);
        splitter.split();
        if ( PRINT_GRAPHS ) graph.printGraph(new File(Utils.join("_", strings) + "_" + hasTop + "_" + hasBot + ".original.dot"), 0);
        if ( PRINT_GRAPHS ) splitter.getSplitGraph().printGraph(new File(Utils.join("_", strings) + "_" + hasTop + "_" + hasBot + ".split.dot"), 0);
        splitter.updateGraph(top, bot);
        if ( PRINT_GRAPHS ) graph.printGraph(new File(Utils.join("_", strings) + "_" + hasTop + "_" + hasBot + ".updated.dot"), 0);

        final KBestHaplotypeFinder splitPaths = new KBestHaplotypeFinder(graph,graph.getSources(),graph.getSinks());
        for ( final KBestHaplotype path : splitPaths ) {
            final String h = new String(path.bases());
            Assert.assertTrue(haplotypes.contains(h), "Failed to find haplotype " + h);
        }


        final List<byte[]> sortedOriginalPaths = new ArrayList<>(originalPaths.size());
        for (final KBestHaplotype kbh : unique(originalPaths))
            sortedOriginalPaths.add(kbh.bases());
        Collections.sort(sortedOriginalPaths, BaseUtils.BASES_COMPARATOR);
        final List<byte[]> sortedSplitPaths = new ArrayList<>(splitPaths.size());
        for (final KBestHaplotype kbh : unique(splitPaths))
            sortedSplitPaths.add(kbh.bases());
        Collections.sort(sortedSplitPaths, BaseUtils.BASES_COMPARATOR);

        Assert.assertEquals(sortedSplitPaths, sortedOriginalPaths, Utils.join("_", strings) + "_" + hasTop + "_" + hasBot);
    }

    /**
     * Returns a unique list of haplotypes solutions.
     * <p>
     *     The result will not contain more than one haplotype with the same base sequence. The solution of the best
     *     score is returned.
     * </p>
     * <p>
     *     This makes sense when there are more than one possible path through the graph to create the same haplotype.
     * </p>
     * <p>
     *     The resulting list is sorted by the score with more likely haplotype search results first.
     * </p>
     *
     * @return never {@code null}, perhaps an empty list.
     */
    public static List<KBestHaplotype> unique(final KBestHaplotypeFinder kbhf) {
        final int requiredCapacity = kbhf.size();
        final Set<Haplotype> haplotypes = new HashSet<>(requiredCapacity);
        final List<KBestHaplotype> result = new ArrayList<>(requiredCapacity);
        for (final KBestHaplotype kbh : kbhf) {
            if (haplotypes.add(kbh.haplotype())) {
                result.add(kbh);
            }
        }
        return result;
    }


    @DataProvider(name = "MeetsMinSequenceData")
    public Object[][] makeMeetsMinSequenceData() {
        final List<Object[]> tests = new ArrayList<>();

        final boolean prefixBiased = true;
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 0, true, true});
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 1, prefixBiased, ! prefixBiased});
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 2, prefixBiased, ! prefixBiased});
        tests.add(new Object[]{Arrays.asList("AC", "AC"), 3, false, false});
        tests.add(new Object[]{Arrays.asList("A", "AC"), 1, true, false});
        tests.add(new Object[]{Arrays.asList("A", "AC"), 2, false, false});
        tests.add(new Object[]{Arrays.asList("AT", "AC"), 1, true, false});
        tests.add(new Object[]{Arrays.asList("AAT", "AAC"), 1, true, false});
        tests.add(new Object[]{Arrays.asList("AAT", "AAC"), 2, true, false});
        tests.add(new Object[]{Arrays.asList("AAT", "AAC"), 3, false, false});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 1, true, true});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 2, true, true});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 3, false, true});
        tests.add(new Object[]{Arrays.asList("AATCCC", "AACCCC"), 4, false, false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MeetsMinSequenceData")
    public void testSplitterCompleteCycle(final List<String> mids, final int minSeqLength, final boolean prefixMeets, final boolean suffixMeets) {
        final SeqGraph graph = new SeqGraph(11);

        final SeqVertex top = new SeqVertex("AAAAAAAA");
        final SeqVertex bot = new SeqVertex("GGGGGGGG");
        final List<SeqVertex> v = new ArrayList<>();
        for ( final String s : mids ) { v.add(new SeqVertex(s)); }
        graph.addVertices(v.toArray(new SeqVertex[v.size()]));
        graph.addVertices(top, bot);
        for ( final SeqVertex vi : v ) { graph.addEdge(top, vi); graph.addEdge(vi, bot); }

        final SharedVertexSequenceSplitter splitter = new SharedVertexSequenceSplitter(graph, v);
        Assert.assertEquals(splitter.meetsMinMergableSequenceForPrefix(minSeqLength), prefixMeets, "Prefix failed");
        Assert.assertEquals(splitter.meetsMinMergableSequenceForSuffix(minSeqLength), suffixMeets, "Suffix failed");
        Assert.assertEquals(splitter.meetsMinMergableSequenceForEitherPrefixOrSuffix(minSeqLength), suffixMeets || prefixMeets, "Either prefix or suffix failed");
    }
}
