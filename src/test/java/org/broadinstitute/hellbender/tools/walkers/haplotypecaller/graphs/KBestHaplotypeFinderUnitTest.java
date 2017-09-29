package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.base.Strings;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

public final class KBestHaplotypeFinderUnitTest extends BaseTest {

    @Test
    public void testScore(){
        final SeqGraph g = new SeqGraph(3);
        final SeqVertex v1 = new SeqVertex("A");
        final SeqVertex v2 = new SeqVertex("C");
        final SeqVertex v3 = new SeqVertex("G");
        final SeqVertex v4 = new SeqVertex("T");
        final SeqVertex v5 = new SeqVertex("A");
        g.addVertex(v1);   //source
        g.addVertex(v2);
        g.addVertex(v3);
        g.addVertex(v4);
        g.addVertex(v5);
        g.addEdge(v1, v2);
        g.addEdge(v2, v3);
        g.addEdge(v2, v4);
        g.addEdge(v2, v5);
        final KBestHaplotypeFinder finder = new KBestHaplotypeFinder(g);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 3);
        final double score = finder.score("ACG".getBytes());
        Assert.assertEquals(score, -0.47712125471966244);
        final double scoreH = finder.score(new Haplotype("ACG".getBytes()));
        Assert.assertEquals(scoreH, -0.47712125471966244);
    }

    @Test
    public void testCycleRemove(){
        final SeqGraph g = new SeqGraph(3);
        final SeqVertex v1 = new SeqVertex("a");
        final SeqVertex v2 = new SeqVertex("b");
        final SeqVertex v3 = new SeqVertex("c");
        final SeqVertex v4 = new SeqVertex("d");
        g.addVertex(v1);   //source
        g.addVertex(v2);
        g.addVertex(v3);
        g.addVertex(v4);  //sink
        g.addEdge(v1, v2);
        g.addEdge(v2, v3);
        g.addEdge(v3, v2); //cycle
        g.addEdge(v3, v4);
        final KBestHaplotypeFinder finder = new KBestHaplotypeFinder(g);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
    }

    @Test
    public void testCyclePrint() throws FileNotFoundException {
        final SeqGraph g = new SeqGraph(3);
        final SeqVertex v1 = new SeqVertex("a");
        final SeqVertex v2 = new SeqVertex("b");
        final SeqVertex v3 = new SeqVertex("c");
        final SeqVertex v4 = new SeqVertex("d");
        g.addVertex(v1);   //source
        g.addVertex(v2);
        g.addVertex(v3);
        g.addVertex(v4);  //sink
        g.addEdge(v1, v2);
        g.addEdge(v2, v3);
        g.addEdge(v3, v2); //cycle
        g.addEdge(v3, v4);
        final KBestHaplotypeFinder finder = new KBestHaplotypeFinder(g);
        finder.printDOT(new PrintWriter(System.out));     //check not blowing up
        final File tmp = createTempFile("foo", ".dot");
        finder.printDOT(tmp);
        final String fname="fred.dot" ;
        finder.printDOTFile(fname);
        Assert.assertTrue(tmp.exists());
        Assert.assertTrue(new File(fname).exists());
        new File(fname).deleteOnExit();
    }

    @Test
    public void testNoSourceOrSink(){
        final SeqGraph g = new SeqGraph(3);
        final SeqVertex v1 = new SeqVertex("a");
        final SeqVertex v2 = new SeqVertex("b");
        g.addVertex(v1);   //source
        g.addVertex(v2);
        g.addEdge(v1, v2);
        final KBestHaplotypeFinder finder = new KBestHaplotypeFinder(g);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);

        final KBestHaplotypeFinder finder2 = new KBestHaplotypeFinder(g, Collections.emptySet(), Collections.singleton(v2));
        Assert.assertEquals(finder2.sources.size(), 0);
        Assert.assertEquals(finder2.sinks.size(), 1);

        final KBestHaplotypeFinder finder3 = new KBestHaplotypeFinder(g, Collections.singleton(v1), Collections.emptySet());
        Assert.assertEquals(finder3.sources.size(), 1);
        Assert.assertEquals(finder3.sinks.size(), 0);
        Assert.assertEquals(finder3.size(), 0);
    }

    @Test
    public void testDeadNode(){
        final SeqGraph g = new SeqGraph(3);
        final SeqVertex v1 = new SeqVertex("a");
        final SeqVertex v2 = new SeqVertex("b");
        final SeqVertex v3 = new SeqVertex("c");
        final SeqVertex v4 = new SeqVertex("d");
        final SeqVertex v5 = new SeqVertex("e");
        g.addVertex(v1);   //source
        g.addVertex(v2);
        g.addVertex(v3);
        g.addVertex(v4);  //sink
        g.addVertex(v5);  //sink
        g.addEdge(v1, v2);
        g.addEdge(v2, v3);
        g.addEdge(v3, v2);//cycle
        g.addEdge(v2, v5);
        g.addEdge(v1, v4);
        final KBestHaplotypeFinder finder1 = new KBestHaplotypeFinder(g);
        Assert.assertEquals(finder1.sources.size(), 1);
        Assert.assertEquals(finder1.sinks.size(), 2);

        final KBestHaplotypeFinder finder2 = new KBestHaplotypeFinder(g, v1, v4); //v5 is a dead node (can't reach the sink v4)
        Assert.assertEquals(finder2.sources.size(), 1);
        Assert.assertEquals(finder2.sinks.size(), 1);
    }

    @DataProvider(name = "BasicPathFindingData")
    public Object[][] makeBasicPathFindingData() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int nStartNodes : Arrays.asList(1, 2, 3) ) {
            for ( final int nBranchesPerBubble : Arrays.asList(2, 3) ) {
                for ( final int nEndNodes : Arrays.asList(1, 2, 3) ) {
                    tests.add(new Object[]{nStartNodes, nBranchesPerBubble, nEndNodes});
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    private static int weight = 1;
    final Set<SeqVertex> createVertices(final SeqGraph graph, final int n, final SeqVertex source, final SeqVertex target) {
        final List<String> seqs = Arrays.asList("A", "C", "G", "T");
        final Set<SeqVertex> vertices = new LinkedHashSet<>();
        for ( int i = 0; i < n; i++ ) {
            final SeqVertex v = new SeqVertex(seqs.get(i));
            graph.addVertex(v);
            vertices.add(v);
            if ( source != null ) graph.addEdge(source, v, new BaseEdge(false, weight++));
            if ( target != null ) graph.addEdge(v, target, new BaseEdge(false, weight++));
        }
        return vertices;
    }

    @Test(dataProvider = "BasicPathFindingData")
    public void testBasicPathFinding(final int nStartNodes, final int nBranchesPerBubble, final int nEndNodes) {
        final SeqGraph graph = new SeqGraph(11);

        final SeqVertex middleTop = new SeqVertex("GTAC");
        final SeqVertex middleBottom = new SeqVertex("ACTG");
        graph.addVertices(middleTop, middleBottom);
        final Set<SeqVertex> starts = createVertices(graph, nStartNodes, null, middleTop);
        @SuppressWarnings("unused")
        final Set<SeqVertex> bubbles = createVertices(graph, nBranchesPerBubble, middleTop, middleBottom);
        final Set<SeqVertex> ends = createVertices(graph, nEndNodes, middleBottom, null);

        // enumerate all possible paths
        final KBestHaplotypeFinder paths = new KBestHaplotypeFinder(graph, starts, ends);

        final int expectedNumOfPaths = nStartNodes * nBranchesPerBubble * nEndNodes;
        Assert.assertEquals(paths.size(), expectedNumOfPaths, "Didn't find the expected number of paths");

        double lastScore = 0;
        for ( final KBestHaplotype kbh : paths ) {
            final Path<SeqVertex,BaseEdge> path = kbh.path();
            Assert.assertTrue(kbh.score() <= lastScore, "Paths out of order.   Path " + path + " has score " + path.getScore() + " above previous " + lastScore);
            lastScore = kbh.score();
        }


        double lastScoreIter = 0;
        final Iterator<KBestHaplotype> iter = paths.iterator(paths.size());
        while ( iter.hasNext() ) {
            final KBestHaplotype kbh = iter.next();
            final Path<SeqVertex,BaseEdge> path = kbh.path();
            Assert.assertTrue(kbh.score() <= lastScoreIter, "Paths out of order.   Path " + path + " has score " + path.getScore() + " above previous " + lastScore);
            lastScoreIter = kbh.score();
        }

        // get the best path, and make sure it's the same as our optimal path overall
        final Path<SeqVertex,BaseEdge> best = paths.get(0).path();
        final List<KBestHaplotype> justOne = new KBestHaplotypeFinder(graph,starts, ends).subList(0,1);
        Assert.assertEquals(justOne.size(), 1);

        Assert.assertTrue(justOne.get(0).path().pathsAreTheSame(best), "Best path from complete enumerate " + best + " not the same as from k = 1 search " + justOne.get(0));
    }


    @DataProvider(name = "BasicBubbleDataProvider")
    public Object[][] makeBasicBubbleDataProvider() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int refBubbleLength : Arrays.asList(1, 5, 10) ) {
            for ( final int altBubbleLength : Arrays.asList(1, 5, 10) ) {
                tests.add(new Object[]{refBubbleLength, altBubbleLength});
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "BasicBubbleDataProvider")
    public void testBasicBubbleData(final int refBubbleLength, final int altBubbleLength) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(3);
        final String preRef = "ATGG";
        final String postRef = "GGGGC";

        SeqVertex v = new SeqVertex(preRef);
        SeqVertex v2Ref = new SeqVertex(Strings.repeat("A", refBubbleLength));
        SeqVertex v2Alt = new SeqVertex(Strings.repeat("A", altBubbleLength-1) + "T");
        SeqVertex v3 = new SeqVertex(postRef);

        graph.addVertex(v);
        graph.addVertex(v2Ref);
        graph.addVertex(v2Alt);
        graph.addVertex(v3);
        graph.addEdge(v, v2Ref, new BaseEdge(true, 10));
        graph.addEdge(v2Ref, v3, new BaseEdge(true, 10));
        graph.addEdge(v, v2Alt, new BaseEdge(false, 5));
        graph.addEdge(v2Alt, v3, new BaseEdge(false, 5));

        // Construct the test path
        Path<SeqVertex,BaseEdge> path = new Path<>(v, graph);
        path = new Path<>(path, graph.getEdge(v, v2Alt));
        path = new Path<>(path, graph.getEdge(v2Alt, v3));

        // Construct the actual cigar string implied by the test path
        Cigar expectedCigar = new Cigar();
        expectedCigar.add(new CigarElement(preRef.length(), CigarOperator.M));
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength, CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength, CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(postRef.length(), CigarOperator.M));

        final String ref = preRef + v2Ref.getSequenceString() + postRef;
        Assert.assertEquals(path.calculateCigar(ref.getBytes()).toString(), AlignmentUtils.consolidateCigar(expectedCigar).toString(), "Cigar string mismatch");
    }

    @DataProvider(name = "GetBasesData")
    public Object[][] makeGetBasesData() {
        List<Object[]> tests = new ArrayList<>();

        final List<String> frags = Arrays.asList("ACT", "GAC", "CAT");

        for ( int n = 1; n <= frags.size(); n++ ) {
            for ( final List<String> comb : Utils.makePermutations(frags, n, false) ) {
                tests.add(new Object[]{comb});
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GetBasesData")
    public void testGetBases(final List<String> frags) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(3);

        SeqVertex prev = null;
        for (final String s : frags) {
            SeqVertex v = new SeqVertex(s);
            graph.addVertex(v);
            if ( prev != null )
                graph.addEdge(prev, v);
            prev = v;
        }

        // enumerate all possible paths
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph);
        Assert.assertEquals(paths.size(), 1);
        final Path<SeqVertex,BaseEdge> path = paths.get(0).path();
        Assert.assertEquals(new String(path.getBases()), Utils.join("", frags), "Path doesn't have the expected sequence");
    }

    @DataProvider(name = "TripleBubbleDataProvider")
    public Object[][] makeTripleBubbleDataProvider() {
        final List<Object[]> tests = new ArrayList<>();
        for ( final int refBubbleLength : Arrays.asList(1, 5, 10) ) {
            for ( final int altBubbleLength : Arrays.asList(1, 5, 10) ) {
                for ( final boolean offRefEnding : Arrays.asList(true, false) ) {
                    for ( final boolean offRefBeginning : Arrays.asList(false) ) {
                        tests.add(new Object[]{refBubbleLength, altBubbleLength, offRefBeginning, offRefEnding});
                    }
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TripleBubbleDataProvider")
    public void testTripleBubbleData(final int refBubbleLength, final int altBubbleLength, final boolean offRefBeginning, final boolean offRefEnding) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(11);
        final String preAltOption = "ATCGATCGATCGATCGATCG";
        final String postAltOption = "CCCC";
        final String preRef = "ATGG";
        final String postRef = "GGCCG";
        final String midRef1 = "TTCCT";
        final String midRef2 = "CCCAAAAAAAAAAAA";

        SeqVertex preV = new SeqVertex(preAltOption);
        SeqVertex v = new SeqVertex(preRef);
        SeqVertex v2Ref = new SeqVertex(Strings.repeat("A", refBubbleLength));
        SeqVertex v2Alt = new SeqVertex(Strings.repeat("A", altBubbleLength - 1) + "T");
        SeqVertex v4Ref = new SeqVertex(Strings.repeat("C", refBubbleLength));
        SeqVertex v4Alt = new SeqVertex(Strings.repeat("C", altBubbleLength - 1) + "T");
        SeqVertex v6Ref = new SeqVertex(Strings.repeat("G", refBubbleLength));
        SeqVertex v6Alt = new SeqVertex(Strings.repeat("G", altBubbleLength - 1) + "T");
        SeqVertex v3 = new SeqVertex(midRef1);
        SeqVertex v5 = new SeqVertex(midRef2);
        SeqVertex v7 = new SeqVertex(postRef);
        SeqVertex postV = new SeqVertex(postAltOption);

        final String ref = preRef + v2Ref.getSequenceString() + midRef1 + v4Ref.getSequenceString() + midRef2 + v6Ref.getSequenceString() + postRef;

        graph.addVertex(preV);
        graph.addVertex(v);
        graph.addVertex(v2Ref);
        graph.addVertex(v2Alt);
        graph.addVertex(v3);
        graph.addVertex(v4Ref);
        graph.addVertex(v4Alt);
        graph.addVertex(v5);
        graph.addVertex(v6Ref);
        graph.addVertex(v6Alt);
        graph.addVertex(v7);
        graph.addVertex(postV);
        graph.addEdge(preV, v, new BaseEdge(false, 1));
        graph.addEdge(v, v2Ref, new BaseEdge(true, 10));
        graph.addEdge(v2Ref, v3, new BaseEdge(true, 10));
        graph.addEdge(v, v2Alt, new BaseEdge(false, 5));
        graph.addEdge(v2Alt, v3, new BaseEdge(false, 5));
        graph.addEdge(v3, v4Ref, new BaseEdge(true, 10));
        graph.addEdge(v4Ref, v5, new BaseEdge(true, 10));
        graph.addEdge(v3, v4Alt, new BaseEdge(false, 5));
        graph.addEdge(v4Alt, v5, new BaseEdge(false, 5));
        graph.addEdge(v5, v6Ref, new BaseEdge(true, 11));
        graph.addEdge(v6Ref, v7, new BaseEdge(true, 11));
        graph.addEdge(v5, v6Alt, new BaseEdge(false, 55));
        graph.addEdge(v6Alt, v7, new BaseEdge(false, 55));
        graph.addEdge(v7, postV, new BaseEdge(false, 1));

        // Construct the test path
        Path<SeqVertex,BaseEdge> path = new Path<>( (offRefBeginning ? preV : v), graph);
        if( offRefBeginning )
            path = new Path<>(path, graph.getEdge(preV, v));
        path = new Path<>(path, graph.getEdge(v, v2Alt));
        path = new Path<>(path, graph.getEdge(v2Alt, v3));
        path = new Path<>(path, graph.getEdge(v3, v4Ref));
        path = new Path<>(path, graph.getEdge(v4Ref, v5));
        path = new Path<>(path, graph.getEdge(v5, v6Alt));
        path = new Path<>(path, graph.getEdge(v6Alt, v7));
        if( offRefEnding )
            path = new Path<>(path, graph.getEdge(v7,postV));

        // Construct the actual cigar string implied by the test path
        Cigar expectedCigar = new Cigar();
        if( offRefBeginning ) {
            expectedCigar.add(new CigarElement(preAltOption.length(), CigarOperator.I));
        }
        expectedCigar.add(new CigarElement(preRef.length(), CigarOperator.M));
        // first bubble
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength, CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength, CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(midRef1.length(), CigarOperator.M));
        // second bubble is ref path
        expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        expectedCigar.add(new CigarElement(midRef2.length(), CigarOperator.M));
        // third bubble
        if( refBubbleLength > altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength - altBubbleLength, CigarOperator.D));
            expectedCigar.add(new CigarElement(altBubbleLength, CigarOperator.M));
        } else if ( refBubbleLength < altBubbleLength ) {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
            expectedCigar.add(new CigarElement(altBubbleLength - refBubbleLength, CigarOperator.I));
        } else {
            expectedCigar.add(new CigarElement(refBubbleLength, CigarOperator.M));
        }
        expectedCigar.add(new CigarElement(postRef.length(), CigarOperator.M));
        if( offRefEnding ) {
            expectedCigar.add(new CigarElement(postAltOption.length(), CigarOperator.I));
        }

        Assert.assertEquals(path.calculateCigar(ref.getBytes()).toString(),
                AlignmentUtils.consolidateCigar(expectedCigar).toString(),
                "Cigar string mismatch: ref = " + ref + " alt " + new String(path.getBases()));
    }

    @Test
    public void testIntraNodeInsertionDeletion() {
        // Construct the assembly graph
        final SeqGraph graph = new SeqGraph(11);
        final SeqVertex top = new SeqVertex("T");
        final SeqVertex bot = new SeqVertex("T");
        final SeqVertex alt = new SeqVertex("AAACCCCC");
        final SeqVertex ref = new SeqVertex("CCCCCGGG");

        graph.addVertices(top, bot, alt, ref);
        graph.addEdges(() -> new BaseEdge(true, 1), top, ref, bot);
        graph.addEdges(() -> new BaseEdge(false, 1), top, alt, bot);

        @SuppressWarnings("all")
        final KBestHaplotypeFinder bestPathFinder = new KBestHaplotypeFinder(graph,top,bot);

        Assert.assertEquals(bestPathFinder.size(), 2);

        final Path<SeqVertex,BaseEdge> refPath = bestPathFinder.get(0).path();
        final Path<SeqVertex,BaseEdge> altPath = bestPathFinder.get(1).path();

        final String refString = top.getSequenceString() + ref.getSequenceString() + bot.getSequenceString();
        Assert.assertEquals(refPath.calculateCigar(refString.getBytes()).toString(), "10M");
        Assert.assertEquals(altPath.calculateCigar(refString.getBytes()).toString(), "1M3I5M3D1M");
    }

    @Test
    public void testHardSWPath() {
        // Construct the assembly graph
        final SeqGraph graph = new SeqGraph(11);
        final SeqVertex top = new SeqVertex( "NNN" );
        final SeqVertex bot = new SeqVertex( "NNN" );
        final SeqVertex alt = new SeqVertex( "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA" );
        final SeqVertex ref = new SeqVertex( "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA" );
        graph.addVertices(top, bot, alt, ref);
        graph.addEdges(() -> new BaseEdge(true, 1), top, ref, bot);
        graph.addEdges(() -> new BaseEdge(false, 1), top, alt, bot);

        @SuppressWarnings("all")
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph, top, bot);

        Assert.assertEquals(paths.size(), 2);

        final Path<SeqVertex,BaseEdge> refPath = paths.get(0).path();
        final Path<SeqVertex,BaseEdge> altPath = paths.get(1).path();

        final String refString = top.getSequenceString() + ref.getSequenceString() + bot.getSequenceString();

        logger.warn("RefPath : " + refPath + " cigar " + refPath.calculateCigar(refString.getBytes()));
        logger.warn("AltPath : " + altPath + " cigar " + altPath.calculateCigar(refString.getBytes()));

        Assert.assertEquals(refPath.calculateCigar(refString.getBytes()).toString(), "51M");
        Assert.assertEquals(altPath.calculateCigar(refString.getBytes()).toString(), "3M6I48M");
    }

    // -----------------------------------------------------------------
    //
    // Systematic tests to ensure that we get the correct SW result for
    // a variety of variants in the ref vs alt bubble
    //
    // -----------------------------------------------------------------

    @DataProvider(name = "SystematicRefAltSWTestData")
    public Object[][] makeSystematicRefAltSWTestData() {
        final List<Object[]> tests = new ArrayList<>();

        final List<List<String>> allDiffs = Arrays.asList(
                Arrays.asList("G", "C", "1M"),
                Arrays.asList("G", "", "1D"),
                Arrays.asList("", "C", "1I"),
                Arrays.asList("AAA", "CGT", "3M"),
                Arrays.asList("TAT", "CAC", "3M"),
                Arrays.asList("GCTG", "GTCG", "4M"),
                Arrays.asList("AAAAA", "", "5D"),
                Arrays.asList("", "AAAAA", "5I"),
                Arrays.asList("AAAAACC", "CCGGGGGG", "5D2M6I")
        );

        for ( final String prefix : Arrays.asList("", "X", "XXXXXXXXXXXXX")) {
            for ( final String end : Arrays.asList("", "X", "XXXXXXXXXXXXX")) {
                for ( final List<String> diffs : allDiffs )
                    tests.add(new Object[]{prefix, end, diffs.get(0), diffs.get(1), diffs.get(2)});
            }
        }

        return tests.toArray(new Object[][]{});
    }


    /**
     * Convenience constructor for testing that creates a path through vertices in graph
     */
    private static <T extends BaseVertex, E extends BaseEdge> Path<T,E> makePath(final List<T> vertices, final BaseGraph<T, E> graph) {
        Path<T,E> path = new Path<>(vertices.get(0), graph);
        for ( int i = 1; i < vertices.size(); i++ ) {
            path = new Path<>(path, graph.getEdge(path.getLastVertex(), vertices.get(i)));
        }
        return path;
    }

    @Test(dataProvider = "SystematicRefAltSWTestData")
    public void testRefAltSW(final String prefix, final String end, final String refMid, final String altMid, final String midCigar) {
        // Construct the assembly graph
        SeqGraph graph = new SeqGraph(11);

        final int padSize = 0;
        SeqVertex top = new SeqVertex(Strings.repeat("N", padSize));
        SeqVertex ref = new SeqVertex(prefix + refMid + end);
        SeqVertex alt = new SeqVertex(prefix + altMid + end);
        SeqVertex bot = new SeqVertex(Strings.repeat("N", padSize));

        graph.addVertices(top, ref, alt, bot);
        graph.addEdges(() -> new BaseEdge(true, 1), top, ref, bot);
        graph.addEdges(() -> new BaseEdge(false, 1), top, alt, bot);

        // Construct the test path
        Path<SeqVertex,BaseEdge> path = makePath(Arrays.asList(top, alt, bot), graph);

        Cigar expected = new Cigar();
        expected.add(new CigarElement(padSize, CigarOperator.M));
        if ( ! prefix.equals("") ) expected.add(new CigarElement(prefix.length(), CigarOperator.M));
        for ( final CigarElement elt : TextCigarCodec.decode(midCigar).getCigarElements() ) expected.add(elt);
        if ( ! end.equals("") ) expected.add(new CigarElement(end.length(), CigarOperator.M));
        expected.add(new CigarElement(padSize, CigarOperator.M));
        expected = AlignmentUtils.consolidateCigar(expected);

        final String refString = top.getSequenceString() + ref.getSequenceString() + bot.getSequenceString();
        final Cigar pathCigar = path.calculateCigar(refString.getBytes());

        logger.warn("diffs: " + ref + " vs. " + alt + " cigar " + midCigar);
        logger.warn("Path " + path + " with cigar " + pathCigar);
        logger.warn("Expected cigar " + expected);

        Assert.assertEquals(pathCigar, expected, "Cigar mismatch: ref = " + refString + " vs alt = " + new String(path.getBases()));
    }

    @Test
    public void testLeftAlignCigarSequentially() {
        String preRefString = "GATCGATCGATC";
        String postRefString = "TTT";
        String refString = "ATCGAGGAGAGCGCCCCG";
        String indelString1 = "X";
        String indelString2 = "YZ";
        int refIndel1 = 10;
        int refIndel2 = 12;

        for ( final int indelSize1 : Arrays.asList(1, 2, 3, 4) ) {
            for ( final int indelOp1 : Arrays.asList(1, -1) ) {
                for ( final int indelSize2 : Arrays.asList(1, 2, 3, 4) ) {
                    for ( final int indelOp2 : Arrays.asList(1, -1) ) {

                        Cigar expectedCigar = new Cigar();
                        expectedCigar.add(new CigarElement(refString.length(), CigarOperator.M));
                        expectedCigar.add(new CigarElement(indelSize1, (indelOp1 > 0 ? CigarOperator.I : CigarOperator.D)));
                        expectedCigar.add(new CigarElement((indelOp1 < 0 ? refIndel1 - indelSize1 : refIndel1), CigarOperator.M));
                        expectedCigar.add(new CigarElement(refString.length(), CigarOperator.M));
                        expectedCigar.add(new CigarElement(indelSize2 * 2, (indelOp2 > 0 ? CigarOperator.I : CigarOperator.D)));
                        expectedCigar.add(new CigarElement((indelOp2 < 0 ? (refIndel2 - indelSize2) * 2 : refIndel2 * 2), CigarOperator.M));
                        expectedCigar.add(new CigarElement(refString.length(), CigarOperator.M));

                        Cigar givenCigar = new Cigar();
                        givenCigar.add(new CigarElement(refString.length() + refIndel1/2, CigarOperator.M));
                        givenCigar.add(new CigarElement(indelSize1, (indelOp1 > 0 ? CigarOperator.I : CigarOperator.D)));
                        givenCigar.add(new CigarElement((indelOp1 < 0 ? (refIndel1/2 - indelSize1) : refIndel1/2) + refString.length() + refIndel2/2 * 2, CigarOperator.M));
                        givenCigar.add(new CigarElement(indelSize2 * 2, (indelOp2 > 0 ? CigarOperator.I : CigarOperator.D)));
                        givenCigar.add(new CigarElement((indelOp2 < 0 ? (refIndel2/2 - indelSize2) * 2 : refIndel2/2 * 2) + refString.length(), CigarOperator.M));

                        String theRef = preRefString + refString + Strings.repeat(indelString1, refIndel1) + refString + Strings.repeat(indelString2, refIndel2) + refString + postRefString;
                        String theRead = refString + Strings.repeat(indelString1, refIndel1 + indelOp1 * indelSize1) + refString + Strings.repeat(indelString2, refIndel2 + indelOp2 * indelSize2) + refString;

                        Cigar calculatedCigar = CigarUtils.leftAlignCigarSequentially(AlignmentUtils.consolidateCigar(givenCigar), theRef.getBytes(), theRead.getBytes(), preRefString.length(), 0);
                        Assert.assertEquals(AlignmentUtils.consolidateCigar(calculatedCigar).toString(), AlignmentUtils.consolidateCigar(expectedCigar).toString(), "Cigar strings do not match!");
                    }
                }
            }
        }
    }

    @Test(enabled = true)
    public void testLeftAlignCigarSequentiallyAdjacentID() {
        final String ref = "GTCTCTCTCTCTCTCTCTATATATATATATATATTT";
        final String hap = "GTCTCTCTCTCTCTCTCTCTCTATATATATATATTT";
        final Cigar originalCigar = TextCigarCodec.decode("18M4I12M4D2M");

        final Cigar result = CigarUtils.leftAlignCigarSequentially(originalCigar, ref.getBytes(), hap.getBytes(), 0, 0);
        logger.warn("Result is " + result);
        Assert.assertEquals(originalCigar.getReferenceLength(), result.getReferenceLength(), "Reference lengths are different");
    }
}
