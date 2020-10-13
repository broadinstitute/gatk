package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.base.Strings;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.JunctionTreeLinkedDeBruijnGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class JunctionTreeKBestHaplotypeFinderUnitTest extends GATKBaseTest {
    private static final SWParameters DANGLING_END_SW_PARAMETERS = SmithWatermanAlignmentUtils.STANDARD_NGS;
    private static final SWParameters PATH_TO_REFERENCE_SW_PARAMETERS = SmithWatermanAlignmentUtils.NEW_SW_PARAMETERS;

    public static byte[] getBytes(final String alignment) {
        return alignment.replace("-","").getBytes();
    }

    @DataProvider (name = "loopingReferences")
    public static Object[][] loopingReferencesRecoverablehaplotypes() {
        return new Object[][]{
                new Object[]{"GACACACAGTCA", 3, 8, true}, //ACA and CAC are repeated, 9 is enough bases to span the path
                //TODO this fails due to new conservative edge pruning code  new Object[]{"GACACACAGTCA", 3, 5, false}, //ACA and CAC are repeated, 8 is not
                new Object[]{"GACACTCACGCACGG", 3, 6, true}, //CAC repeated thrice, showing that the loops are handled somewhat sanely
                new Object[]{"GACATCGACGG", 3, 11, false}, //GAC repeated twice, starting with GAC, can it recover the reference if it starts on a repeated kmer (should not be recoverable based on our decisions)
                new Object[]{"GACATCGACGG", 3, 6, false}, //With non-unique reference start we fall over for looping structures
                new Object[]{"GACATCGACGG", 3, 6, false}, //reads too short to be resolvable
                new Object[]{"GACATCACATC", 3, 8, true}, //final kmer ATC is repeated, can we recover the reference for repating final kmer

                // Some less complicated cases
                new Object[]{"ACTGTGGGGGGGGGGGGCCGCG", 5, 14, true}, //Just long enough to be resolvable (12 repeated G's with 1 base anchor on either side)
                new Object[]{"ACTGTGGGGGGGGGGGGCCGCG", 5, 12, false}, //Just too short to be resolvable

                new Object[]{"ACTGTTCTCTCTCTCTCCCGCG", 5, 14, true}, //Just long enough to be resolvable
                //TODO this fails due to new conservative edge pruning code new Object[]{"ACTGTTCTCTCTCTCTCCCGCG", 5, 9, false}, //Just too short to be resolvable
        };
    }

    @Test (dataProvider = "loopingReferences")
    public void testRecoveryOfLoopingReferenceSequences(final String ref, final int kmerSize, final int readlength, final boolean resolvable) {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        assembler.addSequence("anonymous", getBytes(ref), true);
        // Add "reads" to the graph
        for (int i = 0; i + readlength <= ref.length(); i ++) {
            assembler.addSequence("anonymous", Arrays.copyOfRange(getBytes(ref), i, i + readlength), false);
        }
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());

        // For all of these loops we expect to recover at least the reference haplotype
        Assert.assertTrue(bestPaths.contains(ref));

        if (resolvable) {
            Assert.assertEquals(bestPaths.size(), 1);
        } else {
            Assert.assertTrue(bestPaths.size() > 1);
        }
    }

    @Test (enabled = false) //TODO this test needs to be reconsidered, this is the "baby thrown out with the bathwater" for ignoring reference weight edges wherever possible,
    // todo                             this loop becomes unresolvable because there is no evidence for the path out of the loop so the chain stops expanding after a point
    // Asserting that for a very degenerate looping case (the only weight resides in the loop (which could happen due to non-unique kmers for unmerged dangling ends)
    // note that this particular case is recovered by virtue of the fact that the reference path itself has weight
    // This test asserts that ref path weight is still taken into accoun
    public void testDegenerateLoopingCase() {
        final String ref = "GACACTCACGCACGG";
        final int kmerSize = 3;
        final int readlength = 6;

        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        assembler.addSequence("anonymous", getBytes(ref), true);
        // Add "reads" to the graph
        for (int i = 0; i + readlength < ref.length(); i ++) {
            assembler.addSequence("anonymous", Arrays.copyOfRange(getBytes(ref), i, i + readlength), false);
        }
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());

        // For all of these loops we expect to recover at least the reference haplotype
        Assert.assertTrue(bestPaths.contains(ref));

        Assert.assertEquals(bestPaths.size(), 5);

    }

    @Test
    public void testRecoveryOfAPathDroppedByJunctionTrees() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "AAAACAC"+"C"+"ATGTGCGG"+"T"+"GGGTT"; // the first site has an interesting graph structure and the second site is used to ensure the graph isinterestingg

        // A simple snip het
        String read1 = "AAAACAC"+"T"+"CTGTGCGG"+"C"+"GGGTT"; // CGAC merges with the below read
        String read2 =                "TGTGCGG"+"A"+"GGGTT"; // doesn't show up in the first graph

        assembler.addSequence("reference", getBytes(ref), true);

        assembler.addSequence("anonymous", getBytes(read1), false);
        assembler.addSequence("anonymous", getBytes(read1), false);
        assembler.addSequence("anonymous", getBytes(read1), false);

        assembler.addSequence("anonymous", getBytes(read2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());

        Assert.assertEquals(bestPaths.size(), 2);
        Assert.assertEquals(bestPaths.get(0), read1);
        Assert.assertEquals(bestPaths.get(1), "AAAACAC"+"T"+"C" + read2); //asserting that the front of read 1 was appended to read 2 for haplotype construction
    }

    @Test
    // This test documents one of the known edge cases in the pivotal edge recovery code where sometimes an edge might be
    // dropped if there doesn't happen to be any result path that connects to the root vertex of this dorpped path
    public void testRecoveryOfAPathDroppedByJunctionTreeFailureDueToStillBeingUnreachable() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "AAAACAC"+"C"+"ATGTGCGG"+"TTTAGAGAG"+"GGGTT"; // the first site has an interesting graph structure and the second site is used to ensure the graph is interestingg

        // A simple snip het
        String read1 = "AAAACAC"+"T"+"CTGTGCGG"+"C"+"GGGTT"; // CGAC merges with the below read
        String read2 =                    "TAGAGTG"+"GGGTT"; // creates a branch off of the uncovered reference path

        assembler.addSequence("reference", getBytes(ref), true);

        assembler.addSequence("anonymous", getBytes(read1), false);
        assembler.addSequence("anonymous", getBytes(read1), false);
        assembler.addSequence("anonymous", getBytes(read1), false);

        assembler.addSequence("anonymous", getBytes(read2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());

        Assert.assertEquals(bestPaths.size(), 1);
        Assert.assertEquals(bestPaths.get(0), read1);
    }

    // TODO this test is disabled because currently we have NOT implemented, so now we are simply using the score of the output path as the
    // TODO heurisitic for stapling the front onto this path. this might be addressed in the future.
    @Test (enabled = false)
    public void testRecoveryOfDroppedPathChoosingMostLikePathDespiteThatPathHavingAWorseScore() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref =            "AAAACAC"+"C"+"ATGTGCGG"+"A"+"TTTAGAGAG"+"C"+"GGGTTCC"+"A"+"GAGAGATATA"+"C"+"GAGTTTTGTTT"; // the first site has an interesting graph structure and the second site is used to ensure the graph isinterestingg

        String completePath1 =  "AAAACAC"+"T"+"ATGTGCGG"+"A"+"TTTAGAGAG"+"A"+"GGGTTCC"+"T"+"GAGAGATATA"+"C"+"GAGTTTTGTTT";
        String completePath2 =  "AAAACAC"+"G"+"ATGTGCGG"+"A"+"TTTAGAGAG"+"A"+"GGGTTCC"+"A"+"GAGAGATATA"+"G"+"GAGTTTTGTTT";
        String incompletePath =                "TGTGCGG"+"C"+"TTTAGAGAG"+"A"+"GGGTTCC"+"A"+"GAGAGATATA"+"G"+"GAGTTTTGTTT"; // njote that this point is a copy of path 2 after its first C variant

        // Ensure that completePath1 will have the highest score by weighting it heavily
        assembler.addSequence("reference", getBytes(ref), true);

        assembler.addSequence("anonymous", getBytes(completePath1), false);
        assembler.addSequence("anonymous", getBytes(completePath1), false);
        assembler.addSequence("anonymous", getBytes(completePath1), false);
        assembler.addSequence("anonymous", getBytes(completePath1), false);
        assembler.addSequence("anonymous", getBytes(completePath1), false);

        // followed by lower weight complete path 2
        assembler.addSequence("anonymous", getBytes(completePath2), false);
        assembler.addSequence("anonymous", getBytes(completePath2), false);
        assembler.addSequence("anonymous", getBytes(completePath2), false);

        // And the incomplete path
        assembler.addSequence("anonymous", getBytes(incompletePath), false);
        assembler.addSequence("anonymous", getBytes(incompletePath), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());


        Assert.assertEquals(bestPaths.size(), 3);
        Assert.assertEquals(bestPaths.get(0), completePath1);
        Assert.assertEquals(bestPaths.get(1), completePath2);
        Assert.assertEquals(bestPaths.get(2), "AAAACAC"+"G"+"A" + incompletePath); // asserting that path 2 (the most similar path) was used for the head of incomplete path
    }

    @Test
    // This test demonstrates a potential source of infinite looping, specifically the chain extension
    public void testDegernerateLoopingDanglingEndInChainGatheringCode() {
        final String ref = "GACACTCACGCACGGCC";
        final String alt = "GACACTCACCCCCCCCC"; // the alt ends with a looping and un-escpaed string that is repetative, this will make a path that has no hope of escape (RUN!)
        final int kmerSize = 5;
        final int readlength = 8;

        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        assembler.addSequence("anonymous", getBytes(ref), true);
        // Add "reads" to the graph
        for (int i = 0; i + readlength < ref.length(); i ++) {
            assembler.addSequence("anonymous", Arrays.copyOfRange(getBytes(alt), i, i + readlength), false);
        }
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();
        assembler.pruneJunctionTrees(1);

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());

        Assert.assertTrue(true, "This case could have infinitely looped because the junction trees were all uninformative and pruned on the poly-C branch. There is no condition to stop because the junction tree was removed by pruning");
    }

    @Test
    // Asserting that for a very degenerate looping case (the only weight resides in the loop (which could happen due to non-unique kmers for unmerged dangling ends)
    // note that this particular case is recovered by virtue of the fact that the reference path itself has weight
    public void testDegernerateLoopingDanglingEnd() {
        final String ref = "GACACTCACGCACGGCC";
        final String alt = "GACACTCACCCCCCCCC"; // the alt ends with a looping and un-escpaed string that is repetative, this will make a path that has no hope of escape (RUN!)
        final int kmerSize = 5;
        final int readlength = 8;

        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(kmerSize);
        assembler.addSequence("anonymous", getBytes(ref), true);
        // Add "reads" to the graph
        for (int i = 0; i + readlength < ref.length(); i ++) {
            assembler.addSequence("anonymous", Arrays.copyOfRange(getBytes(alt), i, i + readlength), false);
        }
        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());

        // For this graph structure, there is no evidence whatsoever that a path follow the reference, unfortunately this means we can currently not recover the loop in the alt
        Assert.assertTrue(bestPaths.isEmpty());

    }

    @Test
    // Due to there being no JT covering the reference, we assert that we are only finding the read supported alt path despite the reference path being a reasonable start path
    public void testJunctionTreeRefusesToFollowReferencePathUnlessThereIsNoChoice() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(4);
        String ref = "AAAACAC"+"CCGA"+"ATGTGGGG"+"A"+"GGGTT"; // the first site has an interesting graph structure and the second site is used to ensure the graph isinterestingg

        // A simple snip het
        String read1 = "AAAACAC"+"TCGA"+"CTGTGGGG"+"C"+"GGGTT"; // CGAC merges with the below read
        String read2 = "AAAACAC"+"TCGA"+"CTGTGGGG"+"C"+"GGGTT"; // CGAC merges with the above read

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(read1), false);
        assembler.addSequence("anonymous", getBytes(read2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());

        Assert.assertEquals(bestPaths.size(), 1);
        Assert.assertEquals(new String(bestPaths.get(0).getBytes()), read1);
    }

    // TODO this test will soon be obsolete, we will change path termination to be based on loops soon enough (disabled because the number was bumped up to 7 rather than 4 which broke this test)
    @Test (enabled = false)
    // This test asserts the current behavior towards junction tree ending because of insufficient junction tree data, if we are past our point of junction tree evidence we will fail
    // to find paths, thus this test documents the limitations therein.
    public void testPathTerminationBasedOnDistanceFromLastHaplotype() {
        final JunctionTreeLinkedDeBruijnGraph assembler1 = new JunctionTreeLinkedDeBruijnGraph(6);
        final JunctionTreeLinkedDeBruijnGraph assembler2 = new JunctionTreeLinkedDeBruijnGraph(6);

        String ref =  "AAACAC"+"C"+"ATGGCGG"+"A"+"GGAGTT"+"T"+"GCTCGAA"+"G"+"GGCGTA"+"C"+"CCTACCT"; // the first site has an interesting graph structure and the second site is used to ensure the graph isinterestingg
        String alt1 = "AAACAC"+"A"+"ATGGCGG"+"T"+"GGAGTT"+"G"+"GCTCGAA"+"A"+"GGCGTA"+"G"+"CCTACCT";
        String alt2 = "AAACAC"+"A"+"ATGGCGG"+"T"+"GGAGTT"+"G"+"GCTCGAA"+"A"+"GGCGTA"+"C"+"CCTACCT"; //No expansion of the last path

        assembler1.addSequence("anonymous", getBytes(ref), true);
        assembler2.addSequence("anonymous", getBytes(ref), true);
        // Add short sequences of reads to generate the appropriate alt paths
        for (int i = 0; i + 6 < ref.length(); i ++) {
            assembler1.addSequence("anonymous", Arrays.copyOfRange(getBytes(alt1), i, i + 7), false);
            assembler2.addSequence("anonymous", Arrays.copyOfRange(getBytes(alt2), i, i + 7), false);
        }

        assembler1.buildGraphIfNecessary();
        assembler2.buildGraphIfNecessary();
        assembler1.generateJunctionTrees();
        assembler2.generateJunctionTrees();
        // Assert that the trees are all pointless
        assembler1.getReadThreadingJunctionTrees(false).values().forEach(tree -> Assert.assertTrue(tree.getRootNode().hasNoEvidence()));
        assembler2.getReadThreadingJunctionTrees(false).values().forEach(tree -> Assert.assertTrue(tree.getRootNode().hasNoEvidence()));

        // We had to close out our paths because none of the junction trees are informative and there are 4 decisions to make
        final List<String> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler1).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());
        Assert.assertEquals(bestPaths.size(), 0);

        // We didn't have to close out our paths because there were < 4 decisions made without junction trees, thus we do recover the path
        final List<String> bestPaths2 = new JunctionTreeKBestHaplotypeFinder<>(assembler2).setJunctionTreeEvidenceWeightThreshold(1)
                .findBestHaplotypes(5).stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toList());
        Assert.assertEquals(bestPaths2.size(), 1);
    }

    @Test
    // This is a test enforcing that the behavior around nodes are both outDegree > 1 while also having downstream children with inDegree > 1.
    // We are asserting that the JunctionTree generated by this case lives on the node itself
    public void testPerfectlyPhasedHaplotypeRecoveryRef() {
        int readlength = 20;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AAACAAG"+"G"+"TTGGGTTCG"+"A"+"GCGGGGTTC"+"T"+"CTCGAAGT"+"T"+"CTTGGTAATAT"+"A"+"GGGGGCCCC"; // Reference with 5 sites all separated by at least kmer size
        String alt1 = "AAACAAG"+"T"+"TTGGGTTCG"+"G"+"GCGGGGTTC"+"A"+"CTCGAAGT"+"C"+"CTTGGTAATAT"+"G"+"GGGGGCCCC"; // Alt with different values for all sites

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(alt1.substring(i, i + readlength)), false);
            assembler.addSequence("anonymous", getBytes(ref.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 10);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes();

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 2);
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that
        Assert.assertTrue(foundHaplotypes.contains(alt1));
        Assert.assertTrue(foundHaplotypes.contains(ref));
    }

    @Test
    // This test is intended to illustrate a problem we know causes probems for junction trees, namely that despite there being
    // evidence for a particular path, if there is insufficient evidence we simply end up with combinatorial expansion like before
    public void testInsufficientJunctionTreeDataCausingCombinatorialExpansion() {
        int readlength = 20;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AAACAAG"+"G"+"TTGGGTTCG"+"A"+"GCGGGGTTC"+"T"+"CTCGAAGT"+"T"+"CTTGGTAATAT"+"A"+"GGGGGCCCC"; // Reference with 5 sites all separated by at least kmer size
        String alt1 = "AAACAAG"+"T"+"TTGGGTTCG"+"G"+"GCGGGGTTC"+"A"+"CTCGAAGT"+"C"+"CTTGGTAATAT"+"G"+"GGGGGCCCC"; // Alt with different values for all sites

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);

        assembler.addSequence("anonymous", getBytes(alt1), false);
        assembler.addSequence("anonymous", getBytes(ref), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 10);

        JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes();

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 32); // 2^5 paths from combinatorial expansion
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that
        Assert.assertTrue(foundHaplotypes.contains(alt1));
        Assert.assertTrue(foundHaplotypes.contains(ref));

    }


    @Test
    // This is a test enforcing that the behavior around nodes are both outDegree > 1 while also having downstream children with inDegree > 1.
    // We are asserting that the JunctionTree generated by this case lives on the node itself
    public void testPerfectlyPhasedHaplotypeRecoveryTwoAlts() {
        int readlength = 20;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AAACAAG"+"G"+"TTGGGTTCG"+"A"+"GCGGGGTTC"+"T"+"CTCGAAGT"+"T"+"CTTGGTAATAT"+"A"+"GGGGGCCCC"; // Reference with 5 sites all separated by at least kmer size
        String alt1 = "AAACAAG"+"T"+"TTGGGTTCG"+"G"+"GCGGGGTTC"+"A"+"CTCGAAGT"+"C"+"CTTGGTAATAT"+"G"+"GGGGGCCCC"; // Alt with different values for all sites
        String alt2 = "AAACAAG"+"T"+"TTGGGTTCG"+"G"+"GCGGGGTTC"+"C"+"CTCGAAGT"+"C"+"CTTGGTAATAT"+"G"+"GGGGGCCCC"; // Alt with one different value from alt1

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(alt1.substring(i, i + readlength)), false);
            assembler.addSequence("anonymous", getBytes(alt2.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 6);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes();

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 2);
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that
        Assert.assertTrue(foundHaplotypes.contains(alt1));
        Assert.assertTrue(foundHaplotypes.contains(alt2));
    }

    @Test
    // Asserting that the reference end base is handled with care... This code has the tendency to cut early
    // TODO this test is disabled
    public void testNonUniqueRefStopPosition() {
        int readlength = 15;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AACTTGGGTGTGTGAAACCCGGGTTGTGTGTGAA"; // The sequence GTGTGTGAA is repeated

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i <= ref.length(); i+=2) {
            assembler.addSequence("anonymous", getBytes(ref.substring(i, i + readlength <= ref.length() ? i + readlength : ref.length())), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(3);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes(5);

        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that the correct reference haplotype actually existed
        Assert.assertTrue(foundHaplotypes.contains(ref));
        // Asserting that we did NOT find the haplotype with 0 loops of the repeat:
        Assert.assertFalse(foundHaplotypes.contains("AACTTGGGTGTGTG"));

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 1);
    }

    @Test
    public void testNonUniqueRefStopPositionNoStopJTToResolveIt() {
        int readlength = 15;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AACTTGGGTGTGTGAAACCCGGGTTGTGTGTGAA"; // The sequence GTGTGTGAA is repeated

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < ref.length() - 2; i++) { //the - 2 means that we don't span to the end
            assembler.addSequence("anonymous", getBytes(ref.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 1);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes(5);

        // We assert that we found all of the haplotypes present in reads and nothing else
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that the correct reference haplotype actually existed
        Assert.assertTrue(foundHaplotypes.contains(ref));
        // Asserting that we did NOT find the haplotype with 0 loops of the repeat:
        Assert.assertFalse(foundHaplotypes.contains("AACTTGGGTGTGTG"));
        //TODO maybe assert something about the lenght... unclear
    }

    @Test
    // Test asserting that the behavior is reasonable when there is read data past the reference end kmer
    public void testReferenceEndWithErrorSpanningPastEndBase() {
        int readlength = 15;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AACTGGGTT"  +  "GCGCGCGTTACCCGT"; // The sequence GTGTGTGAA is repeated
        String alt = "AACTGGGTT"+"T"+"GCGCGCGTTACCCGTTT"; // Alt contig with error spanning past the end of the reference sequence

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < alt.length(); i++) {
            assembler.addSequence("anonymous", getBytes(alt.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 1); //TODO this will become 2 once the change is implemented

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes(5);

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 1);
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that we did NOT find the reference haplotype itself
        Assert.assertFalse(foundHaplotypes.contains(ref));
        // Asserting that we did find the alt haplotype minus the ending bases
        Assert.assertTrue(foundHaplotypes.contains("AACTGGGTT"+"T"+"GCGCGCGTTACCCGT"));
    }

    @Test
    // Test asserting that the behavior is reasonable when there is read data past the reference end kmer
    public void testFullReferencePathRecoveryDespiteReadsNotReachingLastBase() {
        int readlength = 15;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AACTGGGTT"  +  "GCGCGCGTTACCCGT"; // The sequence GTGTGTGAA is repeated
        String alt = "AACTGGGTT"+"T"+"GCGCGCGTTACCC"; // Alt contig that doesn't reach the end of the reference

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < alt.length(); i++) {
            assembler.addSequence("anonymous", getBytes(alt.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 1);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes(5);

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 1);
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that we did NOT find the reference haplotype itself
        Assert.assertFalse(foundHaplotypes.contains(ref));
        // Asserting that we did find the alt haplotype minus the ending bases
        Assert.assertTrue(foundHaplotypes.contains("AACTGGGTT"+"T"+"GCGCGCGTTACCCGT"));
    }

    @Test
    // The read evidence at the loop is short (15 bases) and consequently doesn't span the entire loop.
    public void testOfLoopingReferenceReadsTooShortToRecoverIt() {
        int readlength = 15;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "AAACTTTCGCGGGCCCTTAAACCCGCCCTTAAACCCGCCCTTAAACCGCTGTAAGAAA"; // The sequence GCCCTTAAACCC (12 bases) is repeated

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(ref.substring(i, i + readlength)), false);
            assembler.addSequence("anonymous", getBytes(ref.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        Map<MultiDeBruijnVertex, JunctionTreeLinkedDeBruijnGraph.ThreadingTree> junctionTrees = assembler.getReadThreadingJunctionTrees(false);
        Assert.assertEquals(junctionTrees.size(), 2);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes(5); //NOTE we only ask for 5 haplotypes here since the graph might loop forever...

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 2);
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        // Asserting that the correct reference haplotype actually existed
        Assert.assertTrue(foundHaplotypes.contains(ref));
        // Asserting that we did NOT find the haplotype with 0 loops of the repeat:
        Assert.assertFalse(foundHaplotypes.contains("AAACTTTCGCGGGCCCTTAAACCGCTGTAAGAAA"));
    }

    // TODO this is a lingering and very serious problem that we aim to resolve shortly somehow, unclear as to exaclty how right now
    @Test (enabled = false)
    // This test illustrates a current known issue with the new algorithm, where old junction trees with subranchees that don't have a lot of data
    // are used in place of younger trees that might cointain evidence of new paths. This particular test shows that a variant might be dropped as a result.
    public void testOrphanedSubBranchDueToLackingOldJunctionTree() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(6);

        String ref        = "AAAACAC"+"T"+"ATGTGGGG"+"A"+"GGGTTAA"+"A"+"GTCTGAA";
        String haplotype1 = "AAAACAC"+"G"+"ATGTGGGG"+"T"+"GGGTTAA"+"A"+"GTCTGAA";
        String haplotype2 = "AAAACAC"+"G"+"ATGTGGGG"+"T"+"GGGTTAA"+"C"+"GTCTGAA";
        int longReadLength = 20;
        int shortReadLength = 10;

        // Add the reference a few times to get it into the junction trees
        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(ref), false);
        assembler.addSequence("anonymous", getBytes(ref), false);
        assembler.addSequence("anonymous", getBytes(ref), false);

        // Add haplotype 1 reads that are long enough to catch the last variant site from the first junction tree (22)
        for (int i = 0; i + longReadLength < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(haplotype1.substring(i, i + longReadLength)), false);
        }

        // Add haplotype 2 reads that are long enough to leapfrog the entire haplotype but insufficient to establish the G->T->C path on the first tree
        for (int i = 0; i + shortReadLength < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(haplotype2.substring(i, i + shortReadLength)), false);
        }

        assembler.generateJunctionTrees();

        final List<String> haplotypes = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(2)
                .findBestHaplotypes(5).stream().map(h -> new String(h.getBases())).collect(Collectors.toList());

        // Assert that the reference haplotype is recovered
        Assert.assertTrue(haplotypes.contains(ref));
        // Assert that haplotype1 is recovered
        Assert.assertTrue(haplotypes.contains(haplotype1));
        // Assert that haplotype2 is recovered
        // NOTE: this gets dropped because the oldest junction tree on the T->A path has data for haplotype1 but not haplotype2
        Assert.assertTrue(haplotypes.contains(haplotype2));
    }

    @Test
    // This is a test enforcing that the behavior around nodes are both outDegree > 1 while also having downstream children with inDegree > 1.
    // We are asserting that the JunctionTree generated by this case lives on the node itself
    public void testEdgeCaseInvolvingHighInDegreeAndOutDegreeChars() {
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

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        finder.setJunctionTreeEvidenceWeightThreshold(0);
        Assert.assertEquals(finder.sources.size(), 1);
        Assert.assertEquals(finder.sinks.size(), 1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes();

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 3);
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        Assert.assertTrue(foundHaplotypes.contains(refRead1));
        Assert.assertTrue(foundHaplotypes.contains(read1));
        Assert.assertTrue(foundHaplotypes.contains(read2));
    }


    @Test
    // This is a test enforcing that the behavior around nodes are both outDegree > 1 while also having downstream children with inDegree > 1.
    // We are asserting that the JunctionTree generated by this case lives on the node itself
    public void testGraphRecoveryWhenTreeContainsRepeatedKmers() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "AAATCTTCGGGGGGGGGGGGGGTTTCTGGG"; // the first site has an interesting graph structure and the second site is used to ensurethe graph is intersting

        // A simple snip het
        String refRead = "AAATCTTCGGGGGGGGGGGGGGTTTCTGGG";

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(refRead), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        finder.setJunctionTreeEvidenceWeightThreshold(0);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder.findBestHaplotypes();

        // We assert that we found all of the haplotypes present in reads and nothing else
        Assert.assertEquals(haplotypes.size(), 1);
        Set<String> foundHaplotypes = haplotypes.stream().map(haplotype -> new String(haplotype.getBases())).collect(Collectors.toSet());
        Assert.assertTrue(foundHaplotypes.contains(refRead));
    }

    @Test
    public void testSimpleJunctionTreeIncludeRefInJunctionTreeTwoSites() {
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(5);
        String ref = "GGGAAAT" + "T" + "TCCGGC" + "T" + "CGTTTA"; //Two variant sites in close proximity

        // A simple snip het
        String altAARead1 = "GGGAAAT" + "A" + "TCCGGC" + "A" + "CGTTTA"; // Replaces a T with an A, then a T with a A
        String altAARead2 = "GGGAAAT" + "A" + "TCCGGC" + "A" + "CGTTTA"; // Replaces a T with an A, then a T with a A
        String altTCRead1 = "GGGAAAT" + "T" + "TCCGGC" + "C" + "CGTTTA"; // Keeps the T, then replaces a T with a C
        String altTCRead2 = "GGGAAAT" + "T" + "TCCGGC" + "C" + "CGTTTA"; // Keeps the T, then replaces a T with a C

        assembler.addSequence("anonymous", getBytes(ref), true);
        assembler.addSequence("anonymous", getBytes(altAARead1), false);
        assembler.addSequence("anonymous", getBytes(altAARead2), false);
        assembler.addSequence("anonymous", getBytes(altTCRead1), false);
        assembler.addSequence("anonymous", getBytes(altTCRead2), false);

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();
        assembler.pruneJunctionTrees(0);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder1 = new JunctionTreeKBestHaplotypeFinder<>(assembler);
        Assert.assertEquals(finder1.sources.size(), 1);
        Assert.assertEquals(finder1.sinks.size(), 1);

        List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = finder1.findBestHaplotypes(10);
        System.out.println();
    }

    // Disabled until multi-sink/source edges are supported
    @Test (enabled = false)
    public void testDeadNode(){
        final JunctionTreeLinkedDeBruijnGraph g = new JunctionTreeLinkedDeBruijnGraph(3);
        final MultiDeBruijnVertex v1 = new MultiDeBruijnVertex("a".getBytes());
        final MultiDeBruijnVertex v2 = new MultiDeBruijnVertex("b".getBytes());
        final MultiDeBruijnVertex v3 = new MultiDeBruijnVertex("c".getBytes());
        final MultiDeBruijnVertex v4 = new MultiDeBruijnVertex("d".getBytes());
        final MultiDeBruijnVertex v5 = new MultiDeBruijnVertex("e".getBytes());
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
        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder1 = new JunctionTreeKBestHaplotypeFinder<>(g);
        Assert.assertEquals(finder1.sources.size(), 1);
        Assert.assertEquals(finder1.sinks.size(), 2);

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder2 = new JunctionTreeKBestHaplotypeFinder<>(g, v1, v4); //v5 is a dead node (can't reach the sink v4)
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
    final Set<MultiDeBruijnVertex> createVertices(final JunctionTreeLinkedDeBruijnGraph graph, final int n, final MultiDeBruijnVertex source, final MultiDeBruijnVertex target) {
        final List<String> seqs = Arrays.asList("A", "C", "G", "T");
        final Set<MultiDeBruijnVertex> vertices = new LinkedHashSet<>();
        for ( int i = 0; i < n; i++ ) {
            final MultiDeBruijnVertex v = new MultiDeBruijnVertex(seqs.get(i).getBytes());
            graph.addVertex(v);
            vertices.add(v);
            if ( source != null ) graph.addEdge(source, v, new MultiSampleEdge(false, weight++, 1));
            if ( target != null ) graph.addEdge(v, target, new MultiSampleEdge(false, weight++, 1));
        }
        return vertices;
    }


    @Test(dataProvider = "BasicPathFindingData")
    public void testBasicPathFindingNoJunctionTrees(final int nStartNodes, final int nBranchesPerBubble, final int nEndNodes) {
        final JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(11);

        final MultiDeBruijnVertex middleTop = new MultiDeBruijnVertex("GTAC".getBytes());
        final MultiDeBruijnVertex middleBottom = new MultiDeBruijnVertex("ACTG".getBytes());
        graph.addVertices(middleTop, middleBottom);
        final Set<MultiDeBruijnVertex> starts = createVertices(graph, nStartNodes, null, middleTop);
        @SuppressWarnings("unused")
        final Set<MultiDeBruijnVertex> bubbles = createVertices(graph, nBranchesPerBubble, middleTop, middleBottom);
        final Set<MultiDeBruijnVertex> ends = createVertices(graph, nEndNodes, middleBottom, null);

        final int expectedNumOfPaths = nStartNodes * nBranchesPerBubble * nEndNodes;
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> haplotypes = new JunctionTreeKBestHaplotypeFinder<>(graph, starts, ends, JunctionTreeKBestHaplotypeFinder.DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE, false).findBestHaplotypes();
        Assert.assertEquals(haplotypes.size(), expectedNumOfPaths);
        IntStream.range(1, haplotypes.size()).forEach(n -> Assert.assertTrue(haplotypes.get(n-1).score() >= haplotypes.get(n).score()));
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
        JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(4);
        final String preRef = "ATGG";
        final String postRef = "GCGGC";

        final String ref = preRef + Strings.repeat("A", refBubbleLength) + postRef;
        final String alt = preRef + Strings.repeat("A", altBubbleLength-1) + "T" + postRef;

        graph.addSequence("anonomyous", ref.getBytes(),  true);
        for (int i = 0; i < 5; i++) graph.addSequence("anonomyous", ref.getBytes(), 1, false);
        for (int i = 0; i < 10; i++) graph.addSequence("anonomyous", alt.getBytes(), 1, false);

        graph.buildGraphIfNecessary();
        graph.generateJunctionTrees();

        // Construct the test path
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(graph).findBestHaplotypes(5);
        Assert.assertEquals(bestPaths.size(), 2);
        KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge> path = bestPaths.get(0);

        // Construct the actual cigar string implied by the test path
        final CigarBuilder expectedCigar = new CigarBuilder();
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

        Assert.assertEquals(path.calculateCigar(ref.getBytes(), SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS).toString(), expectedCigar.make().toString(), "Cigar string mismatch");
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
    //TODO figure out what to do here as this test doesn't apply
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
        final CigarBuilder expectedCigar = new CigarBuilder();
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

        Assert.assertEquals(path.calculateCigar(ref.getBytes(), SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS).toString(),
                expectedCigar.make().toString(),
                "Cigar string mismatch: ref = " + ref + " alt " + new String(path.getBases()));
    }

    @Test (enabled = false)
    //TODO this test illustrates a problem with junction trees and dangling end recovery, namely the path that the JT points to
    //TODO is a dead end after dangling tail recovery. This needs to be resolved with either SmithWaterman or by coopting the threading code
    public void testIntraNodeInsertionDeletion() {
        // Construct the assembly graph
        final JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(5);
        final String ref = "TTTT" + "CCCCCGGG" + "TTT";
        final String alt = "TTTT" + "AAACCCCC" + "TTT";

        graph.addSequence("anonymous", getBytes(ref), true);
        graph.addSequence("anonymous", getBytes(ref), false);
        graph.addSequence("anonymous", getBytes(alt), false);
        graph.buildGraphIfNecessary();
        graph.recoverDanglingHeads(0, 0,  true, SmithWatermanJavaAligner.getInstance(), DANGLING_END_SW_PARAMETERS);
        graph.recoverDanglingTails(0, 0,  true, SmithWatermanJavaAligner.getInstance(), DANGLING_END_SW_PARAMETERS);
        graph.generateJunctionTrees();

        @SuppressWarnings("all")
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(graph).findBestHaplotypes();
        Assert.assertEquals(bestPaths.size(), 2);
        final Path<MultiDeBruijnVertex,MultiSampleEdge> refPath = bestPaths.get(0);
        final Path<MultiDeBruijnVertex,MultiSampleEdge> altPath = bestPaths.get(1);

        Assert.assertEquals(refPath.calculateCigar(ref.getBytes(), SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS).toString(), "15M");
        Assert.assertEquals(altPath.calculateCigar(ref.getBytes(), SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS).toString(), "4M3I5M3D3M");
    }

    @Test (enabled = false) //TODO this is disabled due to the k max paths per node optimization not being implemented yet
    /*
     This is a test of what can go awry if path pruning is based on the number of incoming edges to a given vertex.
     An illustration of the graph in this test:
               / ----- top \
              /(33)         \
     refStart -(32) -- mid -- midExt ---------------------------- refEnd
              \(34)        / (1)                               /
               \ ---- bot /- (33) botExt - (15) - botExtTop - /
                                        \ (18)               /
                                         \ ------ botExtBot /

      The expected best paths are (refStart->top->midExt->refEnd), and (refStart->mid->midExt->refEnd) because the bottom
      path is penalized for two extra forks despite it greedily looking the best at the start.

      Because the old behavior used to base pruning on the number of incoming edges < K, the edge (bot->midExt) would be created
      first, then (top -> midExt) but (mid -> midExt) would never be created because midExt already has 2 incoming edges.
      */

    public void testDegeneratePathPruningOptimizationCase() {
        // Construct an assembly graph demonstrating this issue
        final SeqGraph graph = new SeqGraph(11);
        final SeqVertex top = new SeqVertex("T");
        final SeqVertex mid = new SeqVertex("C");
        final SeqVertex midAndTopExt = new SeqVertex("GGG");
        final SeqVertex bot = new SeqVertex("G");
        final SeqVertex botExt = new SeqVertex("AAA");
        final SeqVertex botExtTop = new SeqVertex("A");
        final SeqVertex botExtBot = new SeqVertex("T");

        // Ref source and sink vertexes
        final SeqVertex refStart = new SeqVertex("CCCCCGGG");
        final SeqVertex refEnd = new SeqVertex("TTTT");

        graph.addVertices(top, bot, mid, midAndTopExt, bot, botExt, botExtTop, botExtBot, refStart, refEnd);
        // First "diamond" with 3 mostly equivalent cost paths
        graph.addEdges(() -> new BaseEdge(false, 34), refStart, bot, botExt);
        graph.addEdges(() -> new BaseEdge(true, 33), refStart, top, midAndTopExt);
        graph.addEdges(() -> new BaseEdge(false, 32), refStart, mid, midAndTopExt);

        // The best looking path reconnects with a very poor edge multiplicity to midAndTopExt (This is this is what casues the bug)
        graph.addEdges(() -> new BaseEdge(false, 1), bot, midAndTopExt);
        // There is another diamond at bot ext that will end up discounting that path from being in the k best
        graph.addEdges(() -> new BaseEdge(false, 15), botExt, botExtTop, refEnd);
        graph.addEdges(() -> new BaseEdge(false, 18), botExt, botExtBot, refEnd);

        // Wheras the path is smooth sailing from refEnd
        graph.addEdges(() -> new BaseEdge(true, 65), midAndTopExt, refEnd);

        @SuppressWarnings("all")
        final List<KBestHaplotype<SeqVertex, BaseEdge>> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(graph,refStart,refEnd).findBestHaplotypes(2);
        Assert.assertEquals(bestPaths.size(), 2);
        final Path<SeqVertex,BaseEdge> refPath = bestPaths.get(0);
        final Path<SeqVertex,BaseEdge> altPath = bestPaths.get(1);

        Assert.assertEquals(refPath.getVertices().toArray(new SeqVertex[0]), new SeqVertex[]{refStart, top, midAndTopExt, refEnd});
        Assert.assertEquals(altPath.getVertices().toArray(new SeqVertex[0]), new SeqVertex[]{refStart, mid, midAndTopExt, refEnd});
    }


    @Test
    @SuppressWarnings({"unchecked"})
    //TODO this test will make a good base for testing later realignment if the leading Ns are cut back down
    public void testHardSWPath() {
        // Construct the assembly graph
        final JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(11);
        String ref = "NNNNNNNNNNN"+"TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"+"NNN"; // Alt with one different value from alt1
        String alt = "NNNNNNNNNNN"+"ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"+"NNN"; // Alt with one different value from alt1

        // Generate some reads that do not span the entire active region
        graph.addSequence("anonymous", getBytes(ref), true);
        graph.addSequence("anonymous", getBytes(ref), false);
        graph.addSequence("anonymous", getBytes(alt), false);
        graph.buildGraphIfNecessary();
        graph.recoverDanglingHeads(0, 2, true, SmithWatermanJavaAligner.getInstance(), DANGLING_END_SW_PARAMETERS);
        graph.generateJunctionTrees();

        @SuppressWarnings("all")
        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> finder = new JunctionTreeKBestHaplotypeFinder<>(graph);
        finder.setJunctionTreeEvidenceWeightThreshold(1);
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths = finder.findBestHaplotypes();

        Assert.assertEquals(paths.size(), 2);

        final Path<MultiDeBruijnVertex, MultiSampleEdge> refPath = paths.get(0);
        final Path<MultiDeBruijnVertex, MultiSampleEdge> altPath = paths.get(1);

        logger.warn("RefPath : " + refPath + " cigar " + refPath.calculateCigar(ref.getBytes(),
                SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS));
        logger.warn("AltPath : " + altPath + " cigar " + altPath.calculateCigar(ref.getBytes(),
                SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS));

        Assert.assertEquals(refPath.calculateCigar(ref.getBytes(), SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS).toString(), "59M");
        Assert.assertEquals(altPath.calculateCigar(ref.getBytes(), SmithWatermanJavaAligner.getInstance(), PATH_TO_REFERENCE_SW_PARAMETERS).toString(), "11M6I48M");
    }

    @Test
    public void testKmerGraphSimpleReferenceRecovery() {
        // Construct the assembly graph
        final JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(5);
        final MultiDeBruijnVertex refSource = new MultiDeBruijnVertex( "AAATT".getBytes() );
        final MultiDeBruijnVertex k1 = new MultiDeBruijnVertex( "AATTT".getBytes() );
        final MultiDeBruijnVertex k2 = new MultiDeBruijnVertex( "ATTTG".getBytes() );
        final MultiDeBruijnVertex k3 = new MultiDeBruijnVertex( "TTTGG".getBytes() );
        final MultiDeBruijnVertex k4 = new MultiDeBruijnVertex( "TTGGG".getBytes() );
        final MultiDeBruijnVertex k5 = new MultiDeBruijnVertex( "TGGGC".getBytes() );
        final MultiDeBruijnVertex k6 = new MultiDeBruijnVertex( "GGGCC".getBytes() );
        final MultiDeBruijnVertex k7 = new MultiDeBruijnVertex( "GGCCC".getBytes() );
        final MultiDeBruijnVertex k8 = new MultiDeBruijnVertex( "GCCCT".getBytes() );
        final MultiDeBruijnVertex refEnd = new MultiDeBruijnVertex( "CCCTT".getBytes() );
        graph.addVertices(refSource, k1, k2, k3, k4, k5, k6, k7, k8, refEnd);
        graph.addEdges(() -> new MultiSampleEdge(true, 1, 1), refSource, k1, k2, k3, k4, k5, k6, k7, k8, refEnd);

        @SuppressWarnings("all")
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths = new JunctionTreeKBestHaplotypeFinder<>(graph, refSource, refEnd).findBestHaplotypes();

        Assert.assertEquals(paths.size(), 1);

        final Path<MultiDeBruijnVertex, MultiSampleEdge> refPath = paths.get(0);

        final String refString = "AAATTTGGGCCCTT";

        Assert.assertEquals(refPath.getBases(), refString.getBytes());
    }

    @Test
    public void testKmerGraphSimpleReferenceRecoveryWithSNP() {
        // Construct the assembly graph
        final JunctionTreeLinkedDeBruijnGraph graph = new JunctionTreeLinkedDeBruijnGraph(5);
        final MultiDeBruijnVertex refSource = new MultiDeBruijnVertex( "AAATT".getBytes() );
        final MultiDeBruijnVertex k1 = new MultiDeBruijnVertex( "AATTT".getBytes() );
        final MultiDeBruijnVertex k2 = new MultiDeBruijnVertex( "ATTTG".getBytes() );
        final MultiDeBruijnVertex k3 = new MultiDeBruijnVertex( "TTTGG".getBytes() );
        final MultiDeBruijnVertex k4 = new MultiDeBruijnVertex( "TTGGG".getBytes() );
        final MultiDeBruijnVertex k5 = new MultiDeBruijnVertex( "TGGGC".getBytes() );
        final MultiDeBruijnVertex k6 = new MultiDeBruijnVertex( "GGGCC".getBytes() );
        final MultiDeBruijnVertex k7 = new MultiDeBruijnVertex( "GGCCC".getBytes() );
        final MultiDeBruijnVertex k8 = new MultiDeBruijnVertex( "GCCCT".getBytes() );
        final MultiDeBruijnVertex refEnd = new MultiDeBruijnVertex( "CCCTT".getBytes() );
        final MultiDeBruijnVertex v3 = new MultiDeBruijnVertex( "TTTGC".getBytes() );
        final MultiDeBruijnVertex v4 = new MultiDeBruijnVertex( "TTGCG".getBytes() );
        final MultiDeBruijnVertex v5 = new MultiDeBruijnVertex( "TGCGC".getBytes() );
        final MultiDeBruijnVertex v6 = new MultiDeBruijnVertex( "GCGCC".getBytes() );
        final MultiDeBruijnVertex v7 = new MultiDeBruijnVertex( "CGCCC".getBytes() );
        graph.addVertices(refSource, k1, k2, k3, k4, k5, k6, k7, k8, refEnd, v3, v4, v5, v6, v7);
        graph.addEdges(() -> new MultiSampleEdge(true, 1, 4), refSource, k1, k2, k3, k4, k5, k6, k7, k8, refEnd);
        graph.addEdges(() -> new MultiSampleEdge(false, 1, 3), k2, v3, v4, v5, v6, v7, k8);

        @SuppressWarnings("all")
        final List<KBestHaplotype<MultiDeBruijnVertex, MultiSampleEdge>> paths = new JunctionTreeKBestHaplotypeFinder<>(graph, refSource, refEnd).findBestHaplotypes();

        Assert.assertEquals(paths.size(), 1);

        final Path<MultiDeBruijnVertex, MultiSampleEdge> altPath = paths.get(0);

        final String altString = "AAATTTGCGCCCTT";

        Assert.assertEquals(altPath.getBases(), altString.getBytes());
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


    @Test
    public void testCreateMapOfPivotalEdgesInTopoglogicalOrderBasicExample() {
        int readlength = 20;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref = "TAAACAAG"+"G"+"TTGGGTTCG"+"A"+"GCGGGGTTC"+"T"+"CTCGAAGT"+"T"+"CTTGGTAATAT"+"A"+"GGGGGCCCC"; // Reference with 5 sites all separated by at least kmer size
        String alt1 = "TAAACAAG"+"T"+"TTGGGTTCG"+"G"+"GCGGGGTTC"+"A"+"CTCGAAGT"+"C"+"CTTGGTAATAT"+"G"+"GGGGGCCCC"; // Alt with different values for all sites

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(alt1.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1);

        LinkedHashSet<MultiSampleEdge> pivotalEdges = bestPaths.createMapOfPivotalEdgesInTopologicalOrder();
        Iterator<MultiSampleEdge> edgesInOrder = pivotalEdges.iterator();
        ListIterator<String> expectedEdgeKmerTargetsFOrPivotalEdges = Arrays.asList(new String[]{"AACAAGT","GGTTCGG","GGGTTCA","CGAAGTC","TAATATG"}).listIterator();

        // Asserting that the order of edges is as expected
        Assert.assertEquals(pivotalEdges.size(), 5);
        while (edgesInOrder.hasNext()) {
            MultiSampleEdge nextEdge = edgesInOrder.next();
            String expectedKmerEdge = expectedEdgeKmerTargetsFOrPivotalEdges.next();
            String nextTarget = assembler.getEdgeTarget(nextEdge).getSequenceString();
            Assert.assertEquals(nextTarget, expectedKmerEdge);
        }

    }

    @Test
    // this test asserts that included in pivotal edges are edges that are only reachable by following uncovered reference path
    public void testCreateMapOfPivotalEdgesInTopoglogicalOrderHiddenReferenceSequence() {
        int readlength = 20;
        final JunctionTreeLinkedDeBruijnGraph assembler = new JunctionTreeLinkedDeBruijnGraph(7);
        String ref =  "TAAACAAG"+"G"+"TTGGGTTCG"+"AGCGGT"+"CTCGAAGT"+"T"+"CTTGGTAATAT"+"A"+"GGGGGCCCC"; // Reference with 5 sites all separated by at least kmer size
        String alt1 = "TAAACAAG"+"T"+"TTGGGTTCG"+"GGCGGA"+"CTCGAAGT"+"C"+"CTTGGTAATAT"+"G"+"GGGGGCCCC"; // Alt with different values for all sites
        String alt2 =                        "AGCGGTCT"+"G"+"CGAAGT"+"T"+"CTTGGTAATAT"+"A"+"GGGGGCCCC"; // A snippet of reference sequence starting after an uncovered reference edge

        // Generate some reads that do not span the entire active region
        assembler.addSequence("anonymous", getBytes(ref), true);
        for (int i = 0; i + readlength < ref.length(); i++) {
            assembler.addSequence("anonymous", getBytes(alt1.substring(i, i + readlength)), false);
        }

        for (int i = 0; i + readlength < alt2.length(); i++) {
            assembler.addSequence("anonymous", getBytes(alt2.substring(i, i + readlength)), false);
        }

        assembler.buildGraphIfNecessary();
        assembler.generateJunctionTrees();

        final JunctionTreeKBestHaplotypeFinder<MultiDeBruijnVertex, MultiSampleEdge> bestPaths = new JunctionTreeKBestHaplotypeFinder<>(assembler).setJunctionTreeEvidenceWeightThreshold(1);

        LinkedHashSet<MultiSampleEdge> pivotalEdges = bestPaths.createMapOfPivotalEdgesInTopologicalOrder();
        Iterator<MultiSampleEdge> edgesInOrder = pivotalEdges.iterator();
        ListIterator<String> expectedEdgeKmerTargetsFOrPivotalEdges = Arrays.asList(new String[]{"AACAAGT","GGTTCGG", // edges corresoinding to the start of alt1 path
                "CGGTCTG", // Edge corresponding to the G insertion (This would be dropped if we didn't include reference edges with no support in this calculation)
                "CGAAGTC",
                "TAATATA","TAATATG" // Edges correspoindg to the ref and alt respectively, the order is arbitrary because the shortest path is what is taken to determine distance
        }).listIterator();

        // Asserting that the order of edges is as expected
        Assert.assertEquals(pivotalEdges.size(), 6);
        while (edgesInOrder.hasNext()) {
            MultiSampleEdge nextEdge = edgesInOrder.next();
            String expectedKmerEdge = expectedEdgeKmerTargetsFOrPivotalEdges.next();
            String nextTarget = assembler.getEdgeTarget(nextEdge).getSequenceString();
            Assert.assertEquals(nextTarget, expectedKmerEdge);
        }
    }
}