package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public class JunctionTreeLinkedDeBruijnGraphUnitTest extends BaseTest {

    @DataProvider (name = "loopingReferences")
    public static Object[][] loopingReferences() {
        return new Object[][]{
                new Object[]{"GACACACAGTCA", 3}, //ACA and CAC are repeated
                new Object[]{"GACACTCACGCACGG", 3}, //CAC repeated twice, showing that the loops are handled somewhat sanely
                new Object[]{"GACATCGACGG", 3}, //GAC repeated twice, starting with GAC, can it recover the reference if it starts on a repeated kmer
                new Object[]{"GACATCACATC", 3} //final kmer ATC is repeated, can we recover the reference for repating final kmer
        };
    }

    @Test (dataProvider = "loopingReferences")
    public void testRecoveryOfLoopingReferenceSequences(final String ref, final int kmerSize) {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.buildGraphIfNecessary();
        List<MultiDeBruijnVertex> refVertexes = assembler.getReferencePath(ReadThreadingGraph.TraversalDirection.downwards);
        final StringBuilder builder = new StringBuilder(refVertexes.get(0).getSequenceString());
        refVertexes.stream().skip(1).forEach(v -> builder.append(v.getSuffixString()));
        Assert.assertEquals(builder.toString(), ref);
    }

    @Test
    // This test is intended to assert that a potentially dangerous infinite loop is closed
    public void testDanglingEndRecoveryInfiniteHeadLoop() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        final String ref = "AACTGGGCTAGAGAGCGTT";
        final String loopingDanglingHead = "TTCGAAGGTCGAAG"; // "TCGAA" is repeated, causing dangling head recovery to fail

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(ref), false);
        // Coverage = 3, thus above the prune factor of 2
        assembler.addSequence("loopingUnconecetedHead", getBytes(loopingDanglingHead), false);
        assembler.addSequence("loopingUnconecetedHead", getBytes(loopingDanglingHead), false);
        assembler.addSequence("loopingUnconecetedHead", getBytes(loopingDanglingHead), false);

        // This graph should have generated 3 junction trees (one at GAAAA, one at TCGGG, and one at AATCG)
        assembler.buildGraphIfNecessary();

        //Try dangling head recovery (this might have infinitely looped if we are not careful)
        assembler.recoverDanglingHeads(2, 0, true, SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA));

        Assert.assertTrue(true, "If we reached here than the above code didn't infinitely loop");
    }

    @Test
    // As above, pruning the dangling end recovery path can cause the infinite loop recovery to fall over
    public void testDanglingEndRecoveryInfiniteHeadLoopWithPruning() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        final String ref = "AACTGGGCTAGAGAGCGTT";
        final String loopingDanglingHead = "TTCGAAGGTCGAA"; // "TCGAA" is repeated, causing dangling head recovery to fail
        final String loopingDanglingHeadShortened = "TTCGAAGGTC"; // "TCGAA" is repeated, causing dangling head recovery to fail


        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(ref), false);
        assembler.addSequence("loopingUnconecetedHead", getBytes(loopingDanglingHead), false);
        assembler.addSequence("loopingUnconecetedHead", getBytes(loopingDanglingHead), false);
        // Now the some of the kimers in the infinite loop should have < 3 coverage, causing pruning to kick in
        assembler.addSequence("loopingUnconecetedHead", getBytes(loopingDanglingHeadShortened), false);

        // This graph should have generated 3 junction trees (one at GAAAA, one at TCGGG, and one at AATCG)
        assembler.buildGraphIfNecessary();

        //Try dangling head recovery (this might have infinitely looped if we are not careful)
        assembler.recoverDanglingHeads(3, 0, true, SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA));

        Assert.assertTrue(true, "If we reached here than the above code didn't infinitely loop");
    }

    @Test
    public void testJunctionTreePaperTestExample() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "ACTGATTTCGATGCGATGCGATGCCACGGTGG"; // a loop of length 3 in the middle

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(ref), false);

        // This graph should have generated 3 junction trees (one at GAAAA, one at TCGGG, and one at AATCG)
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();
        Assert.assertTrue(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("TCGAT"))).isPresent());
        Assert.assertTrue(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("GCGAT"))).isPresent());

        Assert.assertEquals(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("TCGAT"))).get().getPathsPresentAsBaseChoiceStrings().get(0), "GGC_");
    }

    @Test
    public void testGraphPruningSanityCheck() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "TTTTTAAAAACCCCCGGGGGATATATCGCGCG"; // the first site has an interesting graph structure and the second site is used to ensurethe graph is intersting

        String altRead1 = "TTTTTGAAAACCCTCGGGGGATAAATCGCGCG"; // 3 SNPs
        String altRead2 = "TTTTTGAAAACCCTCGGGGG";
        String altRead3 = "TTTTTGAAAACCCTCG";
        String altRead4 = "TTTTTGAAAA";
        String altRead5 = "TTTTTG";

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(altRead1), false);
        assembler.addSequence("anonymous", getBytes(altRead2), false);
        assembler.addSequence("anonymous", getBytes(altRead3), false);
        assembler.addSequence("anonymous", getBytes(altRead4), false);
        assembler.addSequence("anonymous", getBytes(altRead5), false);

        // This graph should have generated 3 junction trees (one at GAAAA, one at TCGGG, and one at AATCG)
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();
        Assert.assertTrue(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("GAAAA"))).isPresent());
        Assert.assertTrue(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("TCGGG"))).isPresent());
        Assert.assertTrue(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("AATCG"))).isPresent());

        Assert.assertEquals(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("GAAAA"))).get().getPathsPresentAsBaseChoiceStrings().get(0), "TA_");

        // Now we prune requiring at least 2 bases of support (which should kill all but the first graph and then only the first two nodes of the first graph)
        assembler.pruneJunctionTrees(2);

        Assert.assertTrue(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("GAAAA"))).isPresent());
        Assert.assertFalse(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("TCGGG"))).isPresent());
        Assert.assertFalse(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("AATCG"))).isPresent());

        Assert.assertEquals(assembler.getJunctionTreeForNode(assembler.findKmer( new Kmer("GAAAA"))).get().getPathsPresentAsBaseChoiceStrings().get(0), "T");
    }

    @Test
    // The simplest case where we expect two uninformative junction trees to be constructed
    public void testSimpleJunctionTreeExample() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "GGGGAAATTTCCCGGGG";

        // A simple SNP het
        String refRead1 = "GGGAAATTTCC";
        String refRead2 = "GAAATTTCCCGGG";
        String altRead1 = "GGAAATGTC";
        String altRead2 = "GGAAATGTCCCGGG";

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(refRead1), false);
        assembler.addSequence("anonymous", getBytes(refRead2), false);
        assembler.addSequence("anonymous", getBytes(altRead1), false);
        assembler.addSequence("anonymous", getBytes(altRead2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2);

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree1 = junctionTrees.get(assembler.findKmer(new Kmer("GTCCC")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree2 = junctionTrees.get(assembler.findKmer(new Kmer("TTCCC")));
        Assert.assertNotNull(tree1); // Make sure we actually have tree data for that site
        Assert.assertNotNull(tree2); // Make sure we actually have tree data for that site

        List<String> tree1Choices = tree1.getPathsPresentAsBaseChoiceStrings();
        List<String> tree2Choices = tree2.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(tree1Choices, Collections.singletonList(""));
        Assert.assertEquals(tree2Choices, Collections.singletonList(""));

        // Since the remaining graphs are empty prune them and assert as much
        junctionTrees = assembler.getReadThreadingJunctionTrees(true);
        Assert.assertEquals(junctionTrees.size(), 0);
    }

    @Test
    // This is a test enforcing that the behavior around nodes are both outDegree > 1 while also having downstream children with inDegree > 1.
    // We are asserting that the JunctionTree generated by this case lives on the node itself
    public void testEdgeCaseInvolvingHighInDegreeAndOutDegreeCars() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(4);
        String ref = "AAAACAC"+"CCGA"+"ATGTGGGG"+"A"+"GGGTT"; // the first site has an interesting graph structure and the second site is used to ensurethe graph is intersting

        // A simple snip het
        String refRead1 = "AAAACAC"+"CCGA"+"ATGTGGGG"+"A"+"GGGTT";
        String read1 = "AAAACAC"+"CCGA"+"CTGTGGGG"+"C"+"GGGTT"; // CGAC merges with the below read
        String read2 = "AAAACAC"+"TCGA"+"CTGTGGGG"+"C"+"GGGTT"; // CGAC merges with the above read

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(refRead1), false);
        assembler.addSequence("anonymous", getBytes(read1), false);
        assembler.addSequence("anonymous", getBytes(read2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 6);

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree1 = junctionTrees.get(assembler.findKmer(new Kmer("CCGA")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree2 = junctionTrees.get(assembler.findKmer(new Kmer("TCGA")));
        Assert.assertNotNull(tree1); // Make sure we actually have tree data for that site
        Assert.assertNotNull(tree2); // Make sure we actually have tree data for that site

        List<String> tree1Choices = tree1.getPathsPresentAsBaseChoiceStrings();
        List<String> tree2Choices = tree2.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(new HashSet<>(tree1Choices), new HashSet<>(Arrays.asList("AA_","CC_"))); // The AA and CC correspond to taking A and C out of the current nodes.
        Assert.assertEquals(tree2Choices, Collections.singletonList("C_"));
    }


    @Test
    // This is a test enforcing that the behavior around nodes are both outDegree > 1 while also having downstream children with inDegree > 1.
    // We are asserting that the JunctionTree generated by this case lives on the node itself
    public void testGraphRecoveryWhenKmerIsRepeated() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "AAATCTTCGGGGGGGGGGGGGGTTTCTGGG"; // the first site has an interesting graph structure and the second site is used to ensurethe graph is intersting

        // A simple snip het
        String refRead = "AAATCTTCGGGGGGGGGGGGGGTTTCTGGG";

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(refRead), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2);

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree1 = junctionTrees.get(assembler.findKmer(new Kmer("CGGGG")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree2 = junctionTrees.get(assembler.findKmer(new Kmer("GGGGG")));
        Assert.assertNotNull(tree1); // Make sure we actually have tree data for that site
        Assert.assertNotNull(tree2); // Make sure we actually have tree data for that site


        List<String> tree1Choices = tree1.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(tree1Choices, Collections.singletonList("GGGGGGGGGT_"));

        List<String> tree2Choices = tree2.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(new HashSet<>(tree2Choices), new HashSet<>(Arrays.asList(
                 "GGGGGGGGGT_"
                ,"GGGGGGGGT_"
                ,"GGGGGGGT_"
                ,"GGGGGGT_"
                ,"GGGGGT_"
                ,"GGGGT_"
                ,"GGGT_"
                ,"GGT_"
                ,"GT_"
                ,"T_")));
    }

    @Test
    // The simplest case where we expect two uninformative junction trees to be constructed
    public void testPruneGraph() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "GGGGAAATTTCCCGGGG";

        // A simple snip het
        String refRead1 = "GGGAAATTTCC";
        String refRead2 = "GAAATTTCCCGGG";
        String altRead1 = "GGAAATGTC";
        String altRead2 = "GGAAATGTCCCGGG";

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(refRead1), false);
        assembler.addSequence("anonymous", getBytes(refRead2), false);
        assembler.addSequence("anonymous", getBytes(altRead1), false);
        assembler.addSequence("anonymous", getBytes(altRead2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        assembler.pruneJunctionTrees(0);

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 0);
    }

    @Test
    public void testSimpleJunctionTreeIncludeRefInJunctionTreeNoStopEnd() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "GGGGAAATTTCCCGGGG"; // Ref starts 1 base before/after any of the reads

        // A simple snip het
        String altARead1 = "GGGAAATATCC"; // Replaces a T with an A
        String altARead2 = "GAAATATCCCGGG";
        String altTRead1 = "GGAAATGTC"; // Replaces a T with a G
        String altTRead2 = "GGAAATGTCCCGGG";

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(altTRead1), false);
        assembler.addSequence("anonymous", getBytes(altTRead2), false);
        assembler.addSequence("anonymous", getBytes(altARead1), false);
        assembler.addSequence("anonymous", getBytes(altARead2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2);

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree1 = junctionTrees.get(assembler.findKmer(new Kmer("ATCCC")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree2 = junctionTrees.get(assembler.findKmer(new Kmer("GTCCC")));
        //NOTE there is no "TTCCC" edge here because it only existed as a reference path
        Assert.assertNotNull(tree1); // Make sure we actually have tree data for that site
        Assert.assertNotNull(tree2); // Make sure we actually have tree data for that site

        List<String> tree1Choices = tree1.getPathsPresentAsBaseChoiceStrings();
        List<String> tree2Choices = tree2.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(tree1Choices, Collections.singletonList(""));
        Assert.assertEquals(tree2Choices, Collections.singletonList(""));
    }

    @Test
    // NOTE this test is less useful now that the starting JT has been removed
    public void testSimpleJunctionTreeIncludeRefInJunctionTreeNoWithStartEndChanges() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "GGGAAATTTCCCGGG"; // Reads span to the start/stop of the ref, thus we have to worry about symbolic alleles

        // A simple snip het
        String altARead1 = "GGGAAATATCC"; // Replaces a T with an A
        String altARead2 = "GAAATATCCCGGG";
        String altTRead1 = "GGAAATGTC"; // Replaces a T with a G
        String altTRead2 = "GGAAATGTCCCGGG";

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(altTRead1), false);
        assembler.addSequence("anonymous", getBytes(altTRead2), false);
        assembler.addSequence("anonymous", getBytes(altARead1), false);
        assembler.addSequence("anonymous", getBytes(altARead2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2);

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree starttree = junctionTrees.get(assembler.findKmer(new Kmer("GGGAA")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree1 = junctionTrees.get(assembler.findKmer(new Kmer("ATCCC")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree2 = junctionTrees.get(assembler.findKmer(new Kmer("GTCCC")));
        //NOTE there is no "TTCCC" edge here because it only existed as a reference path
        Assert.assertNotNull(tree1); // Make sure we actually have tree data for that site
        Assert.assertNotNull(tree2); // Make sure we actually have tree data for that site

        List<String> tree1Choices = tree1.getPathsPresentAsBaseChoiceStrings();
        List<String> tree2Choices = tree2.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(tree1Choices, Collections.singletonList("_"));
        Assert.assertEquals(tree2Choices, Collections.singletonList("_"));
    }

    @Test
    public void testSimpleJunctionTreeIncludeRefInJunctionTreeTwoSites() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "GGAAAT"+"T"+"TCCGGC"+"T"+"CGTTTAG"; //Two variant sites in close proximity

        // A simple snip het
        String altAARead1 = "GAAAT"+"A"+"TCCGGC"+"A"+"CGTTTA"; // Replaces a T with an A, then a T with a A
        String altAARead2 = "GAAAT"+"A"+"TCCGGC"+"A"+"CGTTTA"; // Replaces a T with an A, then a T with a A
        String altTCRead1 = "GAAAT"+"T"+"TCCGGC"+"C"+"CGTTTA"; // Keeps the T, then replaces a T with a C
        String altTCRead2 = "GAAAT"+"T"+"TCCGGC"+"C"+"CGTTTA"; // Keeps the T, then replaces a T with a C

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(altAARead1), false);
        assembler.addSequence("anonymous", getBytes(altAARead2), false);
        assembler.addSequence("anonymous", getBytes(altTCRead1), false);
        assembler.addSequence("anonymous", getBytes(altTCRead2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(true);
        Assert.assertEquals(junctionTrees.size(), 2);

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree1 = junctionTrees.get(assembler.findKmer(new Kmer("ATCCG")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree tree2 = junctionTrees.get(assembler.findKmer(new Kmer("TTCCG")));
        //NOTE there is no "TTCCC" edge here because it only existed as a reference path
        Assert.assertNotNull(tree1); // Make sure we actually have tree data for that site
        Assert.assertNotNull(tree2); // Make sure we actually have tree data for that site

        List<String> tree1Choices = tree1.getPathsPresentAsBaseChoiceStrings();
        List<String> tree2Choices = tree2.getPathsPresentAsBaseChoiceStrings();
        // Perfectly phased variants (site 2 has 2 mismatches to the reference)
        Assert.assertEquals(tree1Choices, Collections.singletonList("A"));
        Assert.assertEquals(tree2Choices, Collections.singletonList("C"));
    }

    @Test
    public void testSimpleJuncionTreeLoopingReference() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "ATGGTTAGGGGAAATTTAAATTTAAAGCGCCCCCG"; // AAATT, AATTT, ATTTA, TTTAA, and TTAAA are all repeated in a loop

        assembler.addSequence("reference", getBytes(ref), true);
        for (int i = 0; i + 20 < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(ref.substring(i, i + 20)), false);
        }
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2); //

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree outsideTree = junctionTrees.get(assembler.findKmer(new Kmer("GAAAT")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree insideTree = junctionTrees.get(assembler.findKmer(new Kmer("TAAAT")));
        Assert.assertNotNull(outsideTree); // Make sure we actually have tree data for that site
        Assert.assertNotNull(insideTree); // Make sure we actually have tree data for that site

        List<String> insideChoices = insideTree.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(insideChoices.size(), 1);
        Assert.assertTrue(insideChoices.contains("G"));

        List<String> outsideChoices = outsideTree.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(outsideChoices.size(), 1);
        Assert.assertTrue(outsideChoices.contains("TG"));
    }

    @Test
    public void testSimpleJuncionTreeThriceLoopingReference() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "ATGGTTAGGGGAAATTTAAATTTAAATTTAAAGCGCCCCCG"; // AAATT, AATTT, ATTTA, TTTAA, and TTAAA are all repeated in a loop

        assembler.addSequence("reference", getBytes(ref), true);
        for (int i = 0; i + 23 < ref.length(); i++) {
            // 20 bases should be exactly enough to recover the whole path through the loop (from GGAAAT to TTAAAG)
            assembler.addSequence("anonymous", getBytes(ref.substring(i, i + 23)), false);
        }
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2); //

        JunctionTreeLinkedDeBruijnGraph.ThreadingTree outsideTree = junctionTrees.get(assembler.findKmer(new Kmer("GAAAT")));
        JunctionTreeLinkedDeBruijnGraph.ThreadingTree insideTree = junctionTrees.get(assembler.findKmer(new Kmer("TAAAT")));
        Assert.assertNotNull(outsideTree); // Make sure we actually have tree data for that site
        Assert.assertNotNull(insideTree); // Make sure we actually have tree data for that site

        List<String> insideChoices = insideTree.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(insideChoices.size(), 2);
        Assert.assertTrue(insideChoices.contains("TG"));
        Assert.assertTrue(insideChoices.contains("G"));

        List<String> outsideChoices = outsideTree.getPathsPresentAsBaseChoiceStrings();
        Assert.assertEquals(outsideChoices.size(), 1);
        Assert.assertTrue(outsideChoices.contains("TTG"));
    }

    private static final boolean DEBUG = false;

    public static byte[] getBytes(final String alignment) {
        return alignment.replace("-","").getBytes();
    }

    private void assertNonUniqueKmersInvolveLoops(final JunctionTreeLinkedDeBruijnGraph assembler, String... nonUniques) {
        final Set<String> actual = new HashSet<>();
        assembler.buildGraphIfNecessary();
        for (String kmer : nonUniques) {
            MultiDeBruijnVertex vertex = assembler.findKmer(new Kmer(kmer));
            Assert.assertTrue(assembler.incomingEdgesOf(vertex).size() > 1 || assembler.outgoingEdgesOf(vertex).size() > 1);
        }
    }

    @Test
    public void testSimpleHaplotypeRethreading() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(11);
        final String ref   = "CATGCACTTTAAAACTTGCCTTTTTAACAAGACTTCCAGATG";
        final String alt   = "CATGCACTTTAAAACTTGCCGTTTTAACAAGACTTCCAGATG";
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(alt), false);
        assembler.buildGraphIfNecessary();
        Assert.assertNotEquals(ref.length() - 11 + 1, assembler.vertexSet().size(), "the number of vertex in the graph is the same as if there was no alternative sequence");
        Assert.assertEquals(ref.length() - 11 + 1 + 11, assembler.vertexSet().size(), "the number of vertex in the graph is not the same as if there is an alternative sequence");
        MultiDeBruijnVertex startAlt = assembler.findKmer(new Kmer(alt.getBytes(),20,11));
        Assert.assertNotNull(startAlt);
    }

    @Test
    public void testJunctionTreeErrorCorrection() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(11);
        final String ref            = "AAAAAAAAAAACCCCCC"+"G"+"CCCCCTTTTTTTTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        final String supportedAlt   = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCTTTTTTTTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"A"+"ATATATATAATAT";
        // This path has an undersupported edge that gets pruned, we want to assert that we can recover the proper junction tree weights regardless
        final String unSupportedAlt = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCTTTTTTTTTTTGGGGGGG"+"C"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        // Only provide a single instance of the path so that its middle C variant gets pruned
        assembler.addSequence("anonymous", getBytes(unSupportedAlt), false);
        assembler.buildGraphIfNecessary();
        new LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(2).pruneLowWeightChains(assembler);
        assembler.generateJunctionTrees();

        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> bestPaths =
                new JunctionTreeKBestHaplotypeFinder<>(assembler,assembler.getReferenceSourceVertex(),assembler.getReferenceSinkVertex(), 0, false).findBestHaplotypes(10);

        Assert.assertEquals(bestPaths.size(),2);
        Assert.assertEquals(new String(bestPaths.get(0).getBases()), supportedAlt);
        Assert.assertEquals(new String(bestPaths.get(1).getBases()), "AAAAAAAAAAACCCCCC"+"T"+"CCCCCTTTTTTTTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT");
    }

    @Test
    public void testJunctionTreeErrorCorrectionUnrecoverableWithInsertion() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(11);
        final String ref            = "AAAAAAAAAAACCCCCC"+"G"+"CCCCCTTTTTTTTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        final String supportedAlt   = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCTTTTTTTTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"A"+"ATATATATAATAT";
        // This path has an undersupported edge that gets pruned, we want to assert that we can recover the proper junction tree weights regardless
        final String unSupportedAlt = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCTTTTTTTTTTTGGGGGGG"+"CAAT"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        // Only provide a single instance of the path so that its middle C variant gets pruned
        assembler.addSequence("anonymous", getBytes(unSupportedAlt), false);
        assembler.buildGraphIfNecessary();
        new LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(2).pruneLowWeightChains(assembler);
        assembler.generateJunctionTrees();

        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> bestPaths =
                new JunctionTreeKBestHaplotypeFinder<>(assembler,assembler.getReferenceSourceVertex(),assembler.getReferenceSinkVertex(), 0, false).findBestHaplotypes(10);

        // We only see the supported alt path because the unsupported alt path is never recovered
        Assert.assertEquals(bestPaths.size(),1);
        Assert.assertEquals(new String(bestPaths.get(0).getBases()), supportedAlt);
    }

    @Test
    public void testJunctionTreeErrorCorrectionKmerSizeDistanceAllowance() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(11);
        final String ref            = "AAAAAAAAAAACCCCCC"+"G"+"CCCCCCTTTTTT"+"A"+"TTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        final String supportedAlt   = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCCTTTTTT"+"A"+"TTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"A"+"ATATATATAATAT";
        // This path has two unsupported edges which should be pruned, however they are more than kmer size apart so the junction tree code should still work to recover the connectivity
        final String unSupportedAlt = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCCTTTTTT"+"G"+"TTTTTGGGGGGG"+"T"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        // Only provide a single instance of the path so that its middle C variant gets pruned
        assembler.addSequence("anonymous", getBytes(unSupportedAlt), false);
        assembler.buildGraphIfNecessary();
        new LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(2).pruneLowWeightChains(assembler);
        assembler.generateJunctionTrees();

        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> bestPaths =
                new JunctionTreeKBestHaplotypeFinder<>(assembler,assembler.getReferenceSourceVertex(),assembler.getReferenceSinkVertex(), 0, false).findBestHaplotypes(10);

        Assert.assertEquals(bestPaths.size(),2);
        Assert.assertEquals(new String(bestPaths.get(0).getBases()), supportedAlt);
        Assert.assertEquals(new String(bestPaths.get(1).getBases()), "AAAAAAAAAAACCCCCC"+"T"+"CCCCCCTTTTTT"+"A"+"TTTTTGGGGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT");
    }

    @Test
    public void testJuncitonTreeErrorCorrectionFailingIfWithinKmerSizeForConsistencySake() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(11);
        final String ref            = "AAAAAAAAAAACCCCCC"+"G"+"CCCCCCTTTTTT"+"A"+"TTTTGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        final String supportedAlt   = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCCTTTTTT"+"A"+"TTTTGGGG"+"A"+"GGGGTGTGTGTGTGCCCGTGTGT"+"A"+"ATATATATAATAT";
        // This path has two unsupported edges which should be pruned, however they are less than a kmer size apart and should result in the path being unable to find
        final String unSupportedAlt = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCCTTTTTT"+"G"+"TTTTGGGG"+"T"+"GGGGTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt), false);
        // Only provide a single instance of the path so that its middle C variant gets pruned
        assembler.addSequence("anonymous", getBytes(unSupportedAlt), false);
        assembler.buildGraphIfNecessary();
        new LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(2).pruneLowWeightChains(assembler);
        assembler.generateJunctionTrees();

        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> bestPaths =
                new JunctionTreeKBestHaplotypeFinder<>(assembler,assembler.getReferenceSourceVertex(),assembler.getReferenceSinkVertex(), 0, false).findBestHaplotypes(10);

        // We only see the supported alt path because the unsupported alt path is never recovered
        Assert.assertEquals(bestPaths.size(),1);
        Assert.assertEquals(new String(bestPaths.get(0).getBases()), supportedAlt);
    }

    @Test
    // test asserting that we don't make spurious edges when we are trying to recover error bases
    public void testJunctionTreeErrorCorrectionNotAddingTreesWhenOverTentativeBases() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(11);
        final String ref             = "AAAAAAAAAAACCCCCC"+"G"+"CCCCCCTTTTTT"+"A"+"TTGG"+"A"+"GGG"+"C"+"GTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        final String supportedAlt   = "AAAAAAAAAAACCCCCC"+"T"+"CCCCCCTTTTTT"+"A"+"TTGG"+"C"+"GGG"+"C"+"GTGTGTGTGTGCCCGTGTGT"+"A"+"ATATATATAATAT";

        // This path shold be pruned, but if we aren't careful might end up pairing the first SNP with the middle SNP at this site
        final String unSupportedAltWithError = "AAAAAAAAAAACCCCCC"+"G"+"CCCCCCTTTTTT"+"C"+"TTGG"+"C"+"GGG"+"A"+"GTGTGTGTGTGCCCGTGTGT"+"C"+"ATATATATAATAT";
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(ref), false);
        assembler.addSequence("anonymous", getBytes(ref), false);
        assembler.addSequence("anonymous", getBytes(supportedAlt),  false);
        assembler.addSequence("anonymous", getBytes(supportedAlt),  false);
        assembler.addSequence("anonymous", getBytes(unSupportedAltWithError), 1, false);
        assembler.buildGraphIfNecessary();
        new LowWeightChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(2).pruneLowWeightChains(assembler);
        assembler.generateJunctionTrees();

        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> bestPaths =
                new JunctionTreeKBestHaplotypeFinder<>(assembler,assembler.getReferenceSourceVertex(),assembler.getReferenceSinkVertex(), 0, false).setJunctionTreeEvidenceWeightThreshold(1).findBestHaplotypes(10);

        // We only see the supported alt path because the unsupported alt path is never recovered
        // If we saw 3 that means the undersupported pruned path ended up in the junciton trees and created a path that shouldn't be there
        Assert.assertEquals(bestPaths.size(),2);
        Assert.assertEquals(new String(bestPaths.get(0).getBases()), ref);
        Assert.assertEquals(new String(bestPaths.get(1).getBases()), supportedAlt);
    }

    @Test(enabled = ! DEBUG)
    public void testNonUniqueMiddle() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(3);
        final String ref   = "GACACACAGTCA";
        final String read1 = "GACAC---GTCA";
        final String read2 =   "CAC---GTCA";
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
        assertNonUniqueKmersInvolveLoops(assembler, "ACA", "CAC");
    }

    @Test(enabled = ! DEBUG)
    public void testReadsCreateNonUnique() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(3);
        final String ref   = "GCAC--GTCA"; // CAC is unique
        final String read1 = "GCACACGTCA"; // makes CAC non unique because it has a duplication
        final String read2 =    "CACGTCA"; // shouldn't be allowed to match CAC as start
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
//        assembler.convertToSequenceGraph().printGraph(new File("test.dot"), 0);

        assertNonUniqueKmersInvolveLoops(assembler, "CAC");
        //assertSingleBubble(assembler, ref, "CAAAATCGGG");
    }

    @Test(enabled = ! DEBUG)
    // NOTE that now we do still count the reference as edge multiplicity
    public void testCountingOfStartEdges() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(3);
        final String ref   = "NNNGTCAAA"; // ref has some bases before start
        final String read1 =    "GTCAAA"; // starts at first non N base

        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.buildGraphIfNecessary();
//        assembler.printGraph(new File("test.dot"), 0);

        for ( final MultiSampleEdge edge : assembler.edgeSet() ) {
            final MultiDeBruijnVertex source = assembler.getEdgeSource(edge);
            final MultiDeBruijnVertex target = assembler.getEdgeTarget(edge);
            final boolean headerVertex = source.getSuffix() == 'N' || target.getSuffix() == 'N';
            if ( headerVertex ) {
                Assert.assertEquals(edge.getMultiplicity(), 1, "Bases in the unique reference header should have multiplicity of 0");
            } else {
                Assert.assertEquals(edge.getMultiplicity(), 2, "Should have multiplicity of 1 for any edge outside the ref header but got " + edge + " " + source + " -> " + target);
            }
        }
    }

    @Test(enabled = !DEBUG)
    public void testCountingOfStartEdgesWithMultiplePrefixes() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(3);
        assembler.setIncreaseCountsThroughBranches(true);
        final String ref   = "NNNGTCAXX"; // ref has some bases before start
        final String alt1  = "NNNCTCAXX"; // alt1 has SNP right after N
        final String read  =     "TCAXX"; // starts right after SNP, but merges right before branch

        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(alt1), false);
        assembler.addSequence(getBytes(read), false);
        assembler.buildGraphIfNecessary();
        assembler.printGraph(createTempFile("test",".dot"), 0);

        final List<String> oneCountVertices = Arrays.asList("NNN", "NNC", "NCT", "NNG", "NGT");
        final List<String> twoCountVertices = Arrays.asList("CTC", "TCA", "GTC");

        for ( final MultiSampleEdge edge : assembler.edgeSet() ) {
            final MultiDeBruijnVertex source = assembler.getEdgeSource(edge);
            final MultiDeBruijnVertex target = assembler.getEdgeTarget(edge);
            final int expected;
//            if (source.getSequenceString().equals("GTC") && target.getSequenceString().equals("TCA")) {
//                expected = 1;
//            } else {
                expected = oneCountVertices.contains(target.getSequenceString()) ? 1 : ((twoCountVertices.contains(target.getSequenceString()) ? 2 : 3));
//            }
            Assert.assertEquals(edge.getMultiplicity(), expected, "Bases at edge " + edge + " from " + source + " to " + target + " has bad multiplicity");
        }
    }

    @Test(enabled = !DEBUG)
    public void testCyclesInGraph() {

        // b37 20:12655200-12655850
        final String ref = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTATACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        // SNP at 20:12655528 creates a cycle for small kmers
        final String alt = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTGTACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        final List<GATKRead> reads = new ArrayList<>();
        for ( int index = 0; index < alt.length() - 100; index += 20 ) {
            reads.add(ArtificialReadUtils.createArtificialRead(Arrays.copyOfRange(alt.getBytes(), index, index + 100), Utils.dupBytes((byte) 30, 100), 100 + "M"));
        }

        // test that there are cycles detected for small kmer
        final JunctionTreeLinkedDeBruijnGraph rtgraph25 = new JunctionTreeLinkedDeBruijnGraph(25);
        rtgraph25.addSequence("ref", ref.getBytes(), true);
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        for ( final GATKRead read : reads ) {
            rtgraph25.addRead(read, header);
        }
        rtgraph25.buildGraphIfNecessary();
        Assert.assertTrue(rtgraph25.hasCycles());

        // test that there are no cycles detected for large kmer
        final JunctionTreeLinkedDeBruijnGraph rtgraph75 = new JunctionTreeLinkedDeBruijnGraph(75);
        rtgraph75.addSequence("ref", ref.getBytes(), true);
        for ( final GATKRead read : reads ) {
            rtgraph75.addRead(read, header);
        }
        rtgraph75.buildGraphIfNecessary();
        Assert.assertFalse(rtgraph75.hasCycles());
    }

    @Test(enabled = !DEBUG)
    // Test showing that if a read gets completely clipped in the ReadThreadingAssembler, that the assembly will not crash
    public void testEmptyReadBeingAddedToGraph() {

        // b37 20:12655200-12655850
        final String ref = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTATACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        // SNP at 20:12655528 creates a cycle for small kmers
        final String alt = "CAATTGTCATAGAGAGTGACAAATGTTTCAAAAGCTTATTGACCCCAAGGTGCAGCGGTGCACATTAGAGGGCACCTAAGACAGCCTACAGGGGTCAGAAAAGATGTCTCAGAGGGACTCACACCTGAGCTGAGTTGTGAAGGAAGAGCAGGATAGAATGAGCCAAAGATAAAGACTCCAGGCAAAAGCAAATGAGCCTGAGGGAAACTGGAGCCAAGGCAAGAGCAGCAGAAAAGAGCAAAGCCAGCCGGTGGTCAAGGTGGGCTACTGTGTATGCAGAATGAGGAAGCTGGCCAAGTAGACATGTTTCAGATGATGAACATCCTGTGTACTAGATGCATTGGAACTTTTTTCATCCCCTCAACTCCACCAAGCCTCTGTCCACTCTTGGTACCTCTCTCCAAGTAGACATATTTCAGATCATGAACATCCTGTGTACTAGATGCATTGGAAATTTTTTCATCCCCTCAACTCCACCCAGCCTCTGTCCACACTTGGTACCTCTCTCTATTCATATCTCTGGCCTCAAGGAGGGTATTTGGCATTAGTAAATAAATTCCAGAGATACTAAAGTCAGATTTTCTAAGACTGGGTGAATGACTCCATGGAAGAAGTGAAAAAGAGGAAGTTGTAATAGGGAGACCTCTTCGG";

        final List<GATKRead> reads = new ArrayList<>();
        for ( int index = 0; index < alt.length() - 100; index += 20 ) {
            reads.add(ArtificialReadUtils.createArtificialRead(Arrays.copyOfRange(alt.getBytes(), index, index + 100), Utils.dupBytes((byte) 30, 100), 100 + "M"));
        }
        reads.add(ReadUtils.emptyRead(reads.get(0)));

        // test that there are cycles detected for small kmer
        final JunctionTreeLinkedDeBruijnGraph rtgraph25 = new JunctionTreeLinkedDeBruijnGraph(25);
        rtgraph25.addSequence("ref", ref.getBytes(), true);
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        for ( final GATKRead read : reads ) {
            rtgraph25.addRead(read, header);
        }
        rtgraph25.buildGraphIfNecessary();
    }

    @Test(enabled = false)
    // TODO determine what is causing null bases here
    public void testNsInReadsAreNotUsedForGraph() {

        final int length = 100;
        final byte[] ref = Utils.dupBytes((byte) 'A', length);

        final JunctionTreeLinkedDeBruijnGraph rtgraph = new JunctionTreeLinkedDeBruijnGraph(25);
        rtgraph.addSequence("ref", ref, true);

        // add reads with Ns at any position
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        for ( int i = 0; i < length; i++ ) {
            final byte[] bases = ref.clone();
            bases[i] = 'N';
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, Utils.dupBytes((byte) 30, length), length + "M");
            rtgraph.addRead(read, header);
        }
        rtgraph.buildGraphIfNecessary();

        //final SeqGraph graph = rtgraph.toSequenceGraph();
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths = new GraphBasedKBestHaplotypeFinder<>(rtgraph, rtgraph.getReferenceSourceVertex(), rtgraph.getReferenceSinkVertex()).findBestHaplotypes();
        Assert.assertEquals(paths.size(), 1);
    }

// TODO -- update to use determineKmerSizeAndNonUniques directly
//    @DataProvider(name = "KmerSizeData")
//    public Object[][] makeKmerSizeDataProvider() {
//        List<Object[]> tests = new ArrayList<Object[]>();
//
//        // this functionality can be adapted to provide input data for whatever you might want in your data
//        tests.add(new Object[]{3, 3, 3, Arrays.asList("ACG"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 3, Arrays.asList("CAGACG"), Arrays.asList()});
//
//        tests.add(new Object[]{3, 3, 3, Arrays.asList("AAAAC"), Arrays.asList("AAA")});
//        tests.add(new Object[]{3, 4, 4, Arrays.asList("AAAAC"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 4, Arrays.asList("AAAAC"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 3, Arrays.asList("CAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 4, 4, Arrays.asList("CAAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 4, Arrays.asList("CAAAA"), Arrays.asList()});
//        tests.add(new Object[]{3, 5, 5, Arrays.asList("ACGAAAAACG"), Arrays.asList()});
//
//        for ( int maxSize = 3; maxSize < 20; maxSize++ ) {
//            for ( int dupSize = 3; dupSize < 20; dupSize++ ) {
//                final int expectedSize = Math.min(maxSize, dupSize);
//                final String dup = Utils.dupString("C", dupSize);
//                final List<String> nonUnique = dupSize > maxSize ? Arrays.asList(Utils.dupString("C", maxSize)) : Collections.<String>emptyList();
//                tests.add(new Object[]{3, maxSize, expectedSize, Arrays.asList("ACGT", "A" + dup + "GT"), nonUnique});
//                tests.add(new Object[]{3, maxSize, expectedSize, Arrays.asList("A" + dup + "GT", "ACGT"), nonUnique});
//            }
//        }
//
//        return tests.toArray(new Object[][]{});
//    }
//
//    /**
//     * Example testng test using MyDataProvider
//     */
//    @Test(dataProvider = "KmerSizeData")
//    public void testDynamicKmerSizing(final int min, final int max, final int expectKmer, final List<String> seqs, final List<String> expectedNonUniques) {
//        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(min, max);
//        for ( String seq : seqs ) assembler.addSequence(seq.getBytes(), false);
//        assembler.buildGraphIfNecessary();
//        Assert.assertEquals(assembler.getKmerSize(), expectKmer);
//        assertNonUniqueKmersInvolveLoops(assembler, expectedNonUniques.toArray(new String[]{}));
//    }

    @DataProvider(name = "DanglingTails")
    public Object[][] makeDanglingTailsData() {
        List<Object[]> tests = new ArrayList<>();

        // add 1M to the expected CIGAR because it includes the previous (common) base too
        tests.add(new Object[]{"AAAAAAAAAA", "CAAA", "5M", true, 3});                  // incomplete haplotype
        tests.add(new Object[]{"AAAAAAAAAA", "CAAAAAAAAAA", "1M1I10M", true, 10});     // insertion
        tests.add(new Object[]{"CCAAAAAAAAAA", "AAAAAAAAAA", "1M2D10M", true, 10});    // deletion
        tests.add(new Object[]{"AAAAAAAA", "CAAAAAAA", "9M", true, 7});                // 1 snp
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2});                // several snps
        tests.add(new Object[]{"AAAAA", "C", "1M4D1M", false, -1});                    // funky SW alignment
        tests.add(new Object[]{"AAAAA", "CA", "1M3D2M", false, 1});                    // very little data
        tests.add(new Object[]{"AAAAAAA", "CAAAAAC", "8M", true, -1});                 // ends in mismatch
        tests.add(new Object[]{"AAAAAA", "CGAAAACGAA", "1M2I4M2I2M", false, 0});       // alignment is too complex
        tests.add(new Object[]{"AAAAA", "XXXXX", "1M5I", false, -1});                  // insertion

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingTails")
    public void testDanglingTails(final String refEnd,
                                  final String altEnd,
                                  final String cigar,
                                  final boolean cigarIsGood,
                                  final int mergePointDistanceFromSink) {

        final int kmerSize = 15;

        // construct the haplotypes
        final String commonPrefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        final String ref = commonPrefix + refEnd;
        final String alt = commonPrefix + altEnd;

        // create the graph and populate it
        final JunctionTreeLinkedDeBruijnGraph rtgraph = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        rtgraph.addRead(read, header);
        rtgraph.buildGraphIfNecessary();

        // confirm that we have just a single dangling tail
        MultiDeBruijnVertex altSink = null;
        for ( final MultiDeBruijnVertex v : rtgraph.vertexSet() ) {
            if ( rtgraph.isSink(v) && !rtgraph.isReferenceNode(v) ) {
                Assert.assertTrue(altSink == null, "We found more than one non-reference sink");
                altSink = v;
            }
        }

        Assert.assertTrue(altSink != null, "We did not find a non-reference sink");

        // confirm that the SW alignment agrees with our expectations
        final JunctionTreeLinkedDeBruijnGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstDownwardsReferencePath(altSink, 0, 4, false, SmithWatermanJavaAligner
                .getInstance());

        if ( result == null ) {
            Assert.assertFalse(cigarIsGood);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the goodness of the cigar agrees with our expectations
        Assert.assertEquals(JunctionTreeLinkedDeBruijnGraph.cigarIsOkayToMerge(result.cigar, false, true), cigarIsGood);

        // confirm that the tail merging works as expected
        if ( cigarIsGood ) {
            final int mergeResult = rtgraph.mergeDanglingTail(result);
            Assert.assertTrue(mergeResult == 1 || mergePointDistanceFromSink == -1);

            // confirm that we created the appropriate edge
            if ( mergePointDistanceFromSink >= 0 ) {
                MultiDeBruijnVertex v = altSink;
                for ( int i = 0; i < mergePointDistanceFromSink; i++ ) {
                    if ( rtgraph.inDegreeOf(v) != 1 )
                        Assert.fail("Encountered vertex with multiple edges");
                    v = rtgraph.getEdgeSource(rtgraph.incomingEdgeOf(v));
                }
                Assert.assertTrue(rtgraph.outDegreeOf(v) > 1);
            }
        }
    }

    @Test (enabled = false)
    //TODO this test will need to be reenabled when dangling head recovery is handled with the JTBesthaplotypeFinder
    public void testForkedDanglingEnds() {

        final int kmerSize = 15;

        // construct the haplotypes
        final String commonPrefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        final String refEnd = "GCTAGCTAATCG";
        final String altEnd1 = "ACTAGCTAATCG";
        final String altEnd2 = "ACTAGATAATCG";
        final String ref = commonPrefix + refEnd;
        final String alt1 = commonPrefix + altEnd1;
        final String alt2 = commonPrefix + altEnd2;

        // create the graph and populate it
        final JunctionTreeLinkedDeBruijnGraph rtgraph = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(alt1.getBytes(), Utils.dupBytes((byte) 30, alt1.length()), alt1.length() + "M");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(alt2.getBytes(), Utils.dupBytes((byte) 30, alt2.length()), alt2.length() + "M");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        rtgraph.addRead(read1, header);
        rtgraph.addRead(read2, header);
        rtgraph.buildGraphIfNecessary();

        Assert.assertEquals(rtgraph.getSinks().size(), 3);

        for (final MultiDeBruijnVertex altSink : rtgraph.getSinks()) {
            if (rtgraph.isReferenceNode(altSink)) {
                continue;
            }

            // confirm that the SW alignment agrees with our expectations
            final JunctionTreeLinkedDeBruijnGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstDownwardsReferencePath(altSink,
                    0, 4, true, SmithWatermanJavaAligner.getInstance());
            Assert.assertNotNull(result);
            Assert.assertTrue(JunctionTreeLinkedDeBruijnGraph.cigarIsOkayToMerge(result.cigar, false, true));

            // confirm that the tail merging works as expected
            final int mergeResult = rtgraph.mergeDanglingTail(result);
            Assert.assertTrue(mergeResult == 1);
        }

        final List<String> paths = new JunctionTreeKBestHaplotypeFinder<>(rtgraph).findBestHaplotypes().stream()
                .map(kBestHaplotype -> new String(kBestHaplotype.getBases()))
                .distinct()
                .sorted()
                .collect(Collectors.toList());

        Assert.assertEquals(paths.size(), 3);
        Assert.assertEquals(paths, Arrays.asList(ref, alt1, alt2).stream().sorted().collect(Collectors.toList()));
    }

    @Test
    public void testWholeTailIsInsertion() {
        final JunctionTreeLinkedDeBruijnGraph rtgraph = new JunctionTreeLinkedDeBruijnGraph(10);
        final JunctionTreeLinkedDeBruijnGraph.DanglingChainMergeHelper result = new JunctionTreeLinkedDeBruijnGraph.DanglingChainMergeHelper(null, null, "AXXXXX".getBytes(), "AAAAAA".getBytes(), TextCigarCodec.decode("5I1M"));
        final int mergeResult = rtgraph.mergeDanglingTail(result);
        Assert.assertEquals(mergeResult, 0);
    }

    @Test
    public void testGetBasesForPath() {

        final int kmerSize = 4;
        final String testString = "AATGGGGCAATACTA";

        final JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        graph.addSequence(testString.getBytes(), true);
        graph.buildGraphIfNecessary();

        final List<MultiDeBruijnVertex> vertexes = new ArrayList<>();
        MultiDeBruijnVertex v = graph.getReferenceSourceVertex();
        while ( v != null ) {
            vertexes.add(v);
            v = graph.getNextReferenceVertex(v);
        }

        final String resultForTails = new String(graph.getBasesForPath(vertexes, false));
        Assert.assertEquals(resultForTails, testString.substring(kmerSize - 1));
        final String resultForHeads = new String(graph.getBasesForPath(vertexes, true));
        Assert.assertEquals(resultForHeads, "GTAAGGGCAATACTA");  // because the source node will be reversed
    }

    @DataProvider(name = "DanglingHeads")
    public Object[][] makeDanglingHeadsData() {
        List<Object[]> tests = new ArrayList<>();

        // add 1M to the expected CIGAR because it includes the last (common) base too
        tests.add(new Object[]{"XXTXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true});        // 1 snp
        tests.add(new Object[]{"XXXAACCGGTTACGT", "XAAACCGGTTACGT", "7M", false});         // 1 snp
        tests.add(new Object[]{"XXTXXXXAACCGGTTACGT", "XAACGGTTACGT", "4M1D4M", false});   // deletion
        tests.add(new Object[]{"XXTXXXXAACCGGTTACGT", "AYYCGGTTACGT", "8M", true});        // 2 snps
        tests.add(new Object[]{"XXTXXXXAACCGGTTACGTAA", "AYCYGGTTACGTAA", "9M", true});    // 2 snps
        tests.add(new Object[]{"XXXTXXXAACCGGTTACGT", "AYCGGTTACGT", "7M", true});         // very little data
        tests.add(new Object[]{"XXTXXXXAACCGGTTACGT", "YCCGGTTACGT", "6M", true});         // begins in mismatch

        //TODO these tests are diffiuclt to evaluate, as GraphBasedKBestHaplotypeFinder currently cannot handle looping sequence which poses a problem for testing explicitly looping refrence behavior
//        tests.add(new Object[]{"XXTXXXXAACCGGTTACGTTTACGT", "AAYCGGTTACGT", "8M", true});        // dangling head hanging off of a non-unique reference kmer
//        tests.add(new Object[]{"XXTXXXXAACCGGTTACGTCGGTTACGT", "AAYCGGTTACGT", "8M", true});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingHeads")
    public void testDanglingHeads(final String ref,
                                  final String alt,
                                  final String cigar,
                                  final boolean shouldBeMerged) {

        final int kmerSize = 5;

        // create the graph and populate it
        final JunctionTreeLinkedDeBruijnGraph rtgraph = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        rtgraph.addRead(read, header);
        rtgraph.setMaxMismatchesInDanglingHead(10);
        rtgraph.buildGraphIfNecessary();

        // confirm that we have just a single dangling head
        MultiDeBruijnVertex altSource = null;
        for ( final MultiDeBruijnVertex v : rtgraph.vertexSet() ) {
            if ( rtgraph.isSource(v) && !rtgraph.isReferenceNode(v) ) {
                Assert.assertTrue(altSource == null, "We found more than one non-reference source");
                altSource = v;
            }
        }

        Assert.assertTrue(altSource != null, "We did not find a non-reference source");

        // confirm that the SW alignment agrees with our expectations
        final JunctionTreeLinkedDeBruijnGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstUpwardsReferencePath(altSource, 0, 1, false, SmithWatermanJavaAligner
                .getInstance());

        if ( result == null ) {
            Assert.assertFalse(shouldBeMerged);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the tail merging works as expected
        final int mergeResult = rtgraph.mergeDanglingHeadLegacy(result);
        Assert.assertTrue(mergeResult > 0 || !shouldBeMerged);

        // confirm that we created the appropriate bubble in the graph only if expected
        rtgraph.cleanNonRefPaths();
//        final SeqGraph seqGraph = rtgraph.toSequenceGraph();
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths = new GraphBasedKBestHaplotypeFinder<>(rtgraph, rtgraph.getReferenceSourceVertex(), rtgraph.getReferenceSinkVertex()).findBestHaplotypes();
        Assert.assertEquals(paths.size(), shouldBeMerged ? 2 : 1);
    }

    @Test(dataProvider = "DanglingHeads", enabled = false)
    //TODO this test will need to be enabled when we change the dangling end recovery logic, as of right now definitionally we are
    //TODO incapable of recovering dangling heads and that will have to change eventually
    public void testDanglingHeadsWithJTBestHaplotypeFinder(final String ref,
                                  final String alt,
                                  final String cigar,
                                  final boolean shouldBeMerged) {

        final int kmerSize = 5;

        // create the graph and populate it
        final JunctionTreeLinkedDeBruijnGraph rtgraph = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        rtgraph.addRead(read, header);
        rtgraph.setMaxMismatchesInDanglingHead(10);
        rtgraph.buildGraphIfNecessary();

        // confirm that we have just a single dangling head
        MultiDeBruijnVertex altSource = null;
        for ( final MultiDeBruijnVertex v : rtgraph.vertexSet() ) {
            if ( rtgraph.isSource(v) && !rtgraph.isReferenceNode(v) ) {
                Assert.assertTrue(altSource == null, "We found more than one non-reference source");
                altSource = v;
            }
        }

        Assert.assertTrue(altSource != null, "We did not find a non-reference source");

        // confirm that the SW alignment agrees with our expectations
        final JunctionTreeLinkedDeBruijnGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstUpwardsReferencePath(altSource, 0, 1, false, SmithWatermanJavaAligner
                .getInstance());

        if ( result == null ) {
            Assert.assertFalse(shouldBeMerged);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the tail merging works as expected
        final int mergeResult = rtgraph.mergeDanglingHead(result);
        Assert.assertTrue(mergeResult > 0 || !shouldBeMerged);

        // confirm that we created the appropriate bubble in the graph only if expected
        rtgraph.cleanNonRefPaths();
//        final SeqGraph seqGraph = rtgraph.toSequenceGraph();
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths = new JunctionTreeKBestHaplotypeFinder<>(rtgraph, rtgraph.getReferenceSourceVertex(), rtgraph.getReferenceSinkVertex()).findBestHaplotypes();
        Assert.assertEquals(paths.size(), shouldBeMerged ? 2 : 1);
    }

    /**
     * Find the ending position of the longest uniquely matching
     * run of bases of kmer in seq.
     *
     * for example, if seq = ACGT and kmer is NAC, this function returns 1,2 as we have the following
     * match:
     *
     *  0123
     * .ACGT
     * NAC..
     *
     * @param seq a non-null sequence of bytes
     * @param kmer a non-null kmer
     * @return the ending position and length where kmer matches uniquely in sequence, or null if no
     *         unique longest match can be found
     */
    private static Pair<Integer, Integer> findLongestUniqueSuffixMatch(final byte[] seq, final byte[] kmer) {
        int longestPos = -1;
        int length = 0;
        boolean foundDup = false;

        for ( int i = 0; i < seq.length; i++ ) {
            final int matchSize = JunctionTreeLinkedDeBruijnGraph.longestSuffixMatch(seq, kmer, i);
            if ( matchSize > length ) {
                longestPos = i;
                length = matchSize;
                foundDup = false;
            } else if ( matchSize == length ) {
                foundDup = true;
            }
        }

        return foundDup ? null : Pair.of(longestPos, length);
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "findLongestUniqueMatchData")
    public void testfindLongestUniqueMatch(final String seq, final String kmer, final int start, final int length) {
        // adaptor this code to do whatever testing you want given the arguments start and size
        final Pair<Integer, Integer> actual = findLongestUniqueSuffixMatch(seq.getBytes(), kmer.getBytes());
        if ( start == -1 )
            Assert.assertNull(actual);
        else {
            Assert.assertNotNull(actual);
            Assert.assertEquals((int)actual.getLeft(), start);
            Assert.assertEquals((int)actual.getRight(), length);
        }
    }

    @DataProvider(name = "findLongestUniqueMatchData")
    public Object[][] makefindLongestUniqueMatchData() {
        List<Object[]> tests = new ArrayList<>();

        { // test all edge conditions
            final String ref = "ACGT";
            for ( int start = 0; start < ref.length(); start++ ) {
                for ( int end = start + 1; end <= ref.length(); end++ ) {
                    final String kmer = ref.substring(start, end);
                    tests.add(new Object[]{ref, kmer, end - 1, end - start});
                    tests.add(new Object[]{ref, "N" + kmer, end - 1, end - start});
                    tests.add(new Object[]{ref, "NN" + kmer, end - 1, end - start});
                    tests.add(new Object[]{ref, kmer + "N", -1, 0});
                    tests.add(new Object[]{ref, kmer + "NN", -1, 0});
                }
            }
        }

        { // multiple matches
            final String ref = "AACCGGTT";
            for ( final String alt : Arrays.asList("A", "C", "G", "T") )
                tests.add(new Object[]{ref, alt, -1, 0});
            tests.add(new Object[]{ref, "AA", 1, 2});
            tests.add(new Object[]{ref, "CC", 3, 2});
            tests.add(new Object[]{ref, "GG", 5, 2});
            tests.add(new Object[]{ref, "TT", 7, 2});
        }

        { // complex matches that have unique substrings of lots of parts of kmer in the ref
            final String ref = "ACGTACGTACGT";
            tests.add(new Object[]{ref, "ACGT", -1, 0});
            tests.add(new Object[]{ref, "TACGT", -1, 0});
            tests.add(new Object[]{ref, "GTACGT", -1, 0});
            tests.add(new Object[]{ref, "CGTACGT", -1, 0});
            tests.add(new Object[]{ref, "ACGTACGT", -1, 0});
            tests.add(new Object[]{ref, "TACGTACGT", 11, 9});
            tests.add(new Object[]{ref, "NTACGTACGT", 11, 9});
            tests.add(new Object[]{ref, "GTACGTACGT", 11, 10});
            tests.add(new Object[]{ref, "NGTACGTACGT", 11, 10});
            tests.add(new Object[]{ref, "CGTACGTACGT", 11, 11});
        }

        return tests.toArray(new Object[][]{});
    }
}