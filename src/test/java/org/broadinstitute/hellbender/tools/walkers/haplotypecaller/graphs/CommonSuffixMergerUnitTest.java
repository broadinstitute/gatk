package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class CommonSuffixMergerUnitTest extends BaseTest {
    private static final boolean PRINT_GRAPHS = false;

    @DataProvider(name = "CompleteCycleData")
    public Object[][] makeCompleteCycleData() {
        return makeSplitMergeData(-1);
    }

    public static class SplitMergeData {
        final SeqGraph graph;
        final SeqVertex v;
        final String commonSuffix;

        public SplitMergeData(SeqGraph graph, SeqVertex v, String commonSuffix) {
            this.graph = graph;
            this.v = v;
            this.commonSuffix = commonSuffix;
        }

        @Override
        public String toString() {
            return "SplitMergeData{" +
                    "graph=" + graph +
                    ", v=" + v +
                    ", commonSuffix='" + commonSuffix + '\'' +
                    '}';
        }
    }

    public static Object[][] makeSplitMergeData(final int maxTests) {
        List<Object[]> tests = new ArrayList<>();

        final List<String> bases = Arrays.asList("A", "C", "G", "T");
        for ( final String commonSuffix : Arrays.asList("", "A", "AT") ) {
            for ( final int nBots : Arrays.asList(0, 1, 2) ) {
                for ( final int nMids : Arrays.asList(1, 2, 3) ) {
                    for ( int nTops = 0; nTops < nMids; nTops++ ) {
                        for ( int nTopConnections = 1; nTopConnections <= nMids; nTopConnections++ ) {
                            int multi = 1;
                            final SeqGraph graph = new SeqGraph(11);
                            final SeqVertex v = new SeqVertex("GGGG");
                            graph.addVertex(v);

                            final LinkedList<SeqVertex> tops = new LinkedList<>();
                            final LinkedList<SeqVertex> mids = new LinkedList<>();

                            for ( int i = 0; i < nMids; i++) {
                                final SeqVertex mid = new SeqVertex(bases.get(i) + commonSuffix);
                                graph.addVertex(mid);
                                graph.addEdge(mid, v, new BaseEdge(i == 0, multi++));
                                mids.add(mid);

                                tops.add(new SeqVertex(bases.get(i)));
                            }

                            graph.addVertices(tops);
                            for ( final SeqVertex t : tops ) {
                                for ( int i = 0; i < nTopConnections; i++ ) {
                                    graph.addEdge(t, mids.get(i), new BaseEdge(i == 0, multi++));
                                }
                            }

                            for ( int i = 0; i < nBots; i++ ) {
                                final SeqVertex bot = new SeqVertex(bases.get(i));
                                graph.addVertex(bot);
                                graph.addEdge(v, bot, new BaseEdge(i == 0, multi++));

                            }

                            tests.add(new Object[]{new SplitMergeData(graph, v, commonSuffix)});
                        }
                    }
                }
            }
        }

        final List<Object[]> toUse = maxTests == -1 ? tests : tests.subList(0, Math.min(tests.size(), maxTests));
        return toUse.toArray(new Object[][]{});
    }

    /**
     * Compares KBestHaplotype solutions, first by the haplotype base sequence and the by their score.
     */
    private static final Comparator<KBestHaplotype> KBESTHAPLOTYPE_COMPARATOR = new Comparator<KBestHaplotype>() {

        /**
         * Compares KBestHaplotype solutions, first by the haplotype base sequence and the by their score.
         *
         * @return {@inheritDoc}
         */
        @Override
        public int compare(final KBestHaplotype o1,final KBestHaplotype o2) {
            final int baseCmp = o1.haplotype().getBaseString().compareTo(o2.haplotype().getBaseString());
            if (baseCmp != 0)
                return baseCmp;
            return - Double.compare(o1.score(), o2.score());
        }
    };


    public static void assertSameHaplotypes(final String name, final SeqGraph actual, final SeqGraph original) {
        final KBestHaplotypeFinder originalKBestHaplotypes = new KBestHaplotypeFinder(original,original.getSources(),original.getSinks());
        final KBestHaplotypeFinder actualKBestHaplotypes = new KBestHaplotypeFinder(actual,actual.getSources(),actual.getSinks());
        final List<KBestHaplotype> sortedOriginalKBestHaplotypes = new ArrayList<>(originalKBestHaplotypes);
        Collections.sort(sortedOriginalKBestHaplotypes, KBESTHAPLOTYPE_COMPARATOR);
        final List<KBestHaplotype> sortedActualKBestHaplotypes = new ArrayList<>(actualKBestHaplotypes);
        Collections.sort(sortedActualKBestHaplotypes, KBESTHAPLOTYPE_COMPARATOR);
        try {
            final Set<String> haplotypes = new HashSet<>();

            for (final KBestHaplotype kbh : originalKBestHaplotypes)
                haplotypes.add(new String(kbh.bases()));

            for ( final KBestHaplotype kbh : actualKBestHaplotypes ) {
                final String h = new String(kbh.bases());
                Assert.assertTrue(haplotypes.contains(h), "Failed to find haplotype " + h);
            }

            Assert.assertEquals(sortedActualKBestHaplotypes, sortedOriginalKBestHaplotypes);
        } catch ( AssertionError e ) {
            if ( PRINT_GRAPHS ) original.printGraph(new File(String.format("%s.original.dot", name, actual.vertexSet().size())), 0);
            if ( PRINT_GRAPHS ) actual.printGraph(new File(String.format("%s.actual.dot", name, actual.vertexSet().size())), 0);
            try {
                if ( PRINT_GRAPHS ) originalKBestHaplotypes.printDOTFile(String.format("%s.original.finder.dot", name));
                if ( PRINT_GRAPHS ) actualKBestHaplotypes.printDOTFile(String.format("%s.actual.finder.dot", name));
            } catch (IOException e2) {
                // do nothing.
            }
            throw e;
        }
    }

    @Test(dataProvider = "CompleteCycleData")
    public void testMerging(final SplitMergeData data) {
        final SeqGraph original = data.graph.clone();
        SharedSequenceMerger.merge(data.graph, data.v);
        assertSameHaplotypes(String.format("suffixMerge.%s.%d", data.commonSuffix, data.graph.vertexSet().size()), data.graph, original);
    }

    @Test
    public void testDoesntMergeSourceNodes() {
        final SeqGraph g = new SeqGraph(11);
        final SeqVertex v1 = new SeqVertex("A");
        final SeqVertex v2 = new SeqVertex("A");
        final SeqVertex v3 = new SeqVertex("A");
        final SeqVertex top = new SeqVertex("T");
        final SeqVertex b = new SeqVertex("C");
        g.addVertices(top, v1, v2, v3, top, b);
        g.addEdges(top, v1, b);
        g.addEdges(v2, b); // v2 doesn't have previous node, cannot be merged
        g.addEdges(top, v3, b);
        Assert.assertFalse(SharedSequenceMerger.merge(g, b), "Shouldn't be able to merge shared vertices, when one is a source");
    }
}
