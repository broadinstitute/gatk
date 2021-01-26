package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.Kmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;

public final class ReadThreadingGraphUnitTest extends GATKBaseTest {
    private static final boolean DEBUG = false;

    public static byte[] getBytes(final String alignment) {
        return alignment.replace("-","").getBytes();
    }

    private void assertNonUniques(final ReadThreadingGraph assembler, String... nonUniques) {
        final Set<String> actual = new HashSet<>();
        assembler.buildGraphIfNecessary();
        for ( final Kmer kmer : assembler.getNonUniqueKmers() ) actual.add(new String(kmer.bases()));
        final Set<String> expected = new HashSet<>(Arrays.asList(nonUniques));
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testSimpleHaplotypeRethreading() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(11);
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

    @Test(enabled = ! DEBUG)
    public void testNonUniqueMiddle() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "GACACACAGTCA";
        final String read1 = "GACAC---GTCA";
        final String read2 =   "CAC---GTCA";
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
        assertNonUniques(assembler, "ACA", "CAC");
    }

    @Test(enabled = ! DEBUG)
    public void testReadsCreateNonUnique() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        final String ref   = "GCAC--GTCA"; // CAC is unique
        final String read1 = "GCACACGTCA"; // makes CAC non unique because it has a duplication
        final String read2 =    "CACGTCA"; // shouldn't be allowed to match CAC as start
        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(read1), false);
        assembler.addSequence(getBytes(read2), false);
//        assembler.convertToSequenceGraph().printGraph(new File("test.dot"), 0);

        assertNonUniques(assembler, "CAC");
        //assertSingleBubble(assembler, ref, "CAAAATCGGG");
    }

    @Test(enabled = ! DEBUG)
         public void testCountingOfStartEdges() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
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
                Assert.assertEquals(edge.getMultiplicity(), 1, "Bases in the unique reference header should have multiplicity of 1");
            } else {
                Assert.assertEquals(edge.getMultiplicity(), 2, "Should have multiplicity of 2 for any edge outside the ref header but got " + edge + " " + source + " -> " + target);
            }
        }
    }

    @Test(enabled = !DEBUG)
    public void testCountingOfStartEdgesWithMultiplePrefixes() {
        final ReadThreadingGraph assembler = new ReadThreadingGraph(3);
        assembler.setIncreaseCountsThroughBranches(true);
        final String ref   = "NNNGTCAXX"; // ref has some bases before start
        final String alt1  = "NNNCTCAXX"; // alt1 has SNP right after N
        final String read  =     "TCAXX"; // starts right after SNP, but merges right before branch

        assembler.addSequence(getBytes(ref), true);
        assembler.addSequence(getBytes(alt1), false);
        assembler.addSequence(getBytes(read), false);
        assembler.buildGraphIfNecessary();
        assembler.printGraph(createTempFile("test",".dot"), 0);

        final List<String> oneCountVertices = Arrays.asList("NNN", "NNG", "NNC", "NGT", "NCT");
        final List<String> threeCountVertices = Arrays.asList("CAX", "AXX");

        for ( final MultiSampleEdge edge : assembler.edgeSet() ) {
            final MultiDeBruijnVertex source = assembler.getEdgeSource(edge);
            final MultiDeBruijnVertex target = assembler.getEdgeTarget(edge);
            final int expected = oneCountVertices.contains(target.getSequenceString()) ? 1 : (threeCountVertices.contains(target.getSequenceString()) ? 3 : 2);
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
        final ReadThreadingGraph rtgraph25 = new ReadThreadingGraph(25);
        rtgraph25.addSequence("ref", ref.getBytes(), true);
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        for ( final GATKRead read : reads ) {
            rtgraph25.addRead(read, header);
        }
        rtgraph25.buildGraphIfNecessary();
        Assert.assertTrue(rtgraph25.hasCycles());

        // test that there are no cycles detected for large kmer
        final ReadThreadingGraph rtgraph75 = new ReadThreadingGraph(75);
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
        final ReadThreadingGraph rtgraph25 = new ReadThreadingGraph(25);
        rtgraph25.addSequence("ref", ref.getBytes(), true);
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        for ( final GATKRead read : reads ) {
            rtgraph25.addRead(read, header);
        }
        rtgraph25.buildGraphIfNecessary();
    }

    @Test(enabled = !DEBUG)
    public void testNsInReadsAreNotUsedForGraph() {

        final int length = 100;
        final byte[] ref = Utils.dupBytes((byte) 'A', length);

        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(25);
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

        final SeqGraph graph = rtgraph.toSequenceGraph();
        final List<KBestHaplotype<SeqVertex, BaseEdge>> paths = new GraphBasedKBestHaplotypeFinder<>(graph, graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex()).findBestHaplotypes();
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
//        final ReadThreadingGraph assembler = new ReadThreadingGraph(min, max);
//        for ( String seq : seqs ) assembler.addSequence(seq.getBytes(), false);
//        assembler.buildGraphIfNecessary();
//        Assert.assertEquals(assembler.getKmerSize(), expectKmer);
//        assertNonUniques(assembler, expectedNonUniques.toArray(new String[]{}));
//    }

    @DataProvider(name = "DanglingTails")
    public Object[][] makeDanglingTailsData() {
        List<Object[]> tests = new ArrayList<>();

        // add 1M to the expected CIGAR because it includes the previous (common) base too
        tests.add(new Object[]{"AAAAAAAAAA", "CAAA", "5M", true, 3, -1});                  // incomplete haplotype
        tests.add(new Object[]{"AAAAAAAAAA", "CAAAAAAAAAA", "1M1I10M", true, 10, -1});     // insertion
        tests.add(new Object[]{"CCAAAAAAAAAA", "AAAAAAAAAA", "1M2D10M", true, 10, -1});    // deletion
        tests.add(new Object[]{"AAAAAAAA", "CAAAAAAA", "9M", true, 7, -1});                // 1 snp
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2, -1});                // several snps
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2, 0});                // several snps
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2, 1});                // several snps
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, 2, 2});                // several snps
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, -1, 3});                // several snps (not enough bases to match)
        tests.add(new Object[]{"AAAAAAAA", "CAAGATAA", "9M", true, -1, 4});                // several snps (not enough bases to match)
        tests.add(new Object[]{"AAAAA", "C", "1M4D1M", false, -1, -1});                    // funky SW alignment
        tests.add(new Object[]{"AAAAA", "CA", "1M3D2M", false, 1, -1});                    // very little data
        tests.add(new Object[]{"AAAAAAA", "CAAAAAC", "8M", true, -1, -1});                 // ends in mismatch
        tests.add(new Object[]{"AAAAAA", "CGAAAACGAA", "1M2I4M2I2M", false, 0, -1});       // alignment is too complex
        tests.add(new Object[]{"AAAAA", "XXXXX", "1M5I", false, -1, -1});                  // insertion

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingTails")
    public void testDanglingTails(final String refEnd,
                                  final String altEnd,
                                  final String cigar,
                                  final boolean cigarIsGood,
                                  final int mergePointDistanceFromSink,
                                  final int numLeadingMatchesAllowed) {

        final int kmerSize = 15;

        // construct the haplotypes
        final String commonPrefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        final String ref = commonPrefix + refEnd;
        final String alt = commonPrefix + altEnd;

        // create the graph and populate it
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        rtgraph.addRead(read, header);
        rtgraph.setMinMatchingBasesToDangingEndRecovery(numLeadingMatchesAllowed);
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
        final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstDownwardsReferencePath(altSink, 0, 4, false, SmithWatermanJavaAligner
                .getInstance());

        if ( result == null ) {
            Assert.assertFalse(cigarIsGood);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the goodness of the cigar agrees with our expectations
        Assert.assertEquals(ReadThreadingGraph.cigarIsOkayToMerge(result.cigar, false, true), cigarIsGood);

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

    @Test
    public void testForkedDanglingEndsWithSuffixCode() {

        final int kmerSize = 15;

        // construct the haplotypes
        final String commonPrefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        final String refEnd = "GCTAGCTAATCGTTAAGCTTTAAC";
        final String altEnd1 = "GCTAGCTAA"+ "GGCG";// two mismatches compared to the reference
        final String altEnd2 = "GCTAGCTAA"+ "GCCGATGGCT";
        final String ref = commonPrefix + refEnd;
        final String alt1 = commonPrefix + altEnd1;
        final String alt2 = commonPrefix + altEnd2;

        // create the graph and populate it
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(alt1.getBytes(), Utils.dupBytes((byte) 30, alt1.length()), alt1.length() + "M");
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(alt2.getBytes(), Utils.dupBytes((byte) 30, alt2.length()), alt2.length() + "M");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        rtgraph.setMinMatchingBasesToDangingEndRecovery(1);
        rtgraph.addRead(read2, header);
        rtgraph.addRead(read1, header);
        rtgraph.buildGraphIfNecessary();

        Assert.assertEquals(rtgraph.getSinks().size(), 3);

        // Testing a degenerate case where the wrong "reference" path is selected and assuring that when
        for (final MultiDeBruijnVertex altSink : rtgraph.getSinks()) {
            if (rtgraph.isReferenceNode(altSink) || !altSink.getSequenceString().equals("GCTAAGCCGATGGCT")) {
                continue;
            }

            // confirm that the SW alignment agrees with our expectations
            final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstDownwardsReferencePath(altSink,
                    0, 2, false, SmithWatermanJavaAligner.getInstance());
            Assert.assertNotNull(result);
            Assert.assertTrue(ReadThreadingGraph.cigarIsOkayToMerge(result.cigar, false, true));

            Assert.assertEquals(Math.min(AbstractReadThreadingGraph.longestSuffixMatch(result.referencePathString, result.danglingPathString, result.cigar.getReferenceLength() - 1), result.cigar.getCigarElement(0).getLength()), 0);

            // confirm that the tail merging works as expected
            final int mergeResult = rtgraph.mergeDanglingTail(result);
            Assert.assertTrue(mergeResult == 0);
        }
    }

    @Test
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
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
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
            final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstDownwardsReferencePath(altSink,
                    0, 4, true, SmithWatermanJavaAligner.getInstance());
            Assert.assertNotNull(result);
            Assert.assertTrue(ReadThreadingGraph.cigarIsOkayToMerge(result.cigar, false, true));

            // confirm that the tail merging works as expected
            final int mergeResult = rtgraph.mergeDanglingTail(result);
            Assert.assertTrue(mergeResult == 1);
        }

        final SeqGraph seqGraph = rtgraph.toSequenceGraph();
        seqGraph.simplifyGraph();

        final List<String> paths = new GraphBasedKBestHaplotypeFinder<>(seqGraph).findBestHaplotypes().stream()
                .map(kBestHaplotype -> new String(kBestHaplotype.getBases()))
                .distinct()
                .sorted()
                .collect(Collectors.toList());

        Assert.assertEquals(paths.size(), 3);
        Assert.assertEquals(paths, Arrays.asList(ref, alt1, alt2).stream().sorted().collect(Collectors.toList()));
    }

    @Test
    public void testWholeTailIsInsertion() {
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(10);
        final ReadThreadingGraph.DanglingChainMergeHelper result = new ReadThreadingGraph.DanglingChainMergeHelper(null, null, "AXXXXX".getBytes(), "AAAAAA".getBytes(), TextCigarCodec.decode("5I1M"));
        final int mergeResult = rtgraph.mergeDanglingTail(result);
        Assert.assertEquals(mergeResult, 0);
    }

    @Test
    public void testGetBasesForPath() {

        final int kmerSize = 4;
        final String testString = "AATGGGGCAATACTA";

        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);
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
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, -1});        // 1 snp
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, 0});        // 1 snp
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, 1});        // 1 snp
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", true, 2});        // 1 snp
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AAYCGGTTACGT", "8M", false, 3});        // 1 snp

        // One SNP failing in legacy behavior (failing due to being too close to the reference start
        tests.add(new Object[]{"XXXAACCGGTTACGT", "XAAACCGGTTACGT", "7M", false, -1});         // 1 snp
        tests.add(new Object[]{"XXXAACCGGTTACGT", "XAAACCGGTTACGT", "7M", false, 0});         // 1 snp

        // Testing that indels now work
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "XAACGGTTACGT", "4M1D4M", false, -1});   // deletion
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "XAACGGTTACGT", "4M1D4M", true, 1});   // deletion
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "XAACYCGGTTACGT", "5M1I4M", false, -1});   // insertion
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "XAACYCGGTTACGT", "5M1I4M", true, 1});   // insertion

        // Testing old/new behavior with multiple SNPs
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AYYCGGTTACGT", "8M", false, -1});        // 2 snps (only allow 1 SNP by default)
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AYYCGGTTACGT", "8M", true, 1});        // 2 snps
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGTAA", "AYCYGGTTACGTAA", "9M", false, -1});    // 2 snps (only allow 1 SNP by default)
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGTAA", "AYCYGGTTACGTAA", "9M", true, 1});    // 2 snps

        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "AYCGGTTACGT", "7M", true, -1});         // very little data
        tests.add(new Object[]{"XXXXXXXAACCGGTTACGT", "YCCGGTTACGT", "6M", true, -1});         // begins in mismatch

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DanglingHeads")
    public void testDanglingHeads(final String ref,
                                  final String alt,
                                  final String cigar,
                                  final boolean shouldBeMerged,
                                  final int numLeadingMatchesAllowed) {

        final int kmerSize = 5;

        // create the graph and populate it
        final ReadThreadingGraph rtgraph = new ReadThreadingGraph(kmerSize);
        rtgraph.addSequence("ref", ref.getBytes(), true);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(alt.getBytes(), Utils.dupBytes((byte) 30, alt.length()), alt.length() + "M");
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        rtgraph.addRead(read, header);
        rtgraph.setMinMatchingBasesToDangingEndRecovery(numLeadingMatchesAllowed);
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
        final ReadThreadingGraph.DanglingChainMergeHelper result = rtgraph.generateCigarAgainstUpwardsReferencePath(altSource, 0, 1, false, SmithWatermanJavaAligner
                .getInstance());

        if ( result == null ) {
            Assert.assertFalse(shouldBeMerged);
            return;
        }

        Assert.assertTrue(cigar.equals(result.cigar.toString()), "SW generated cigar = " + result.cigar.toString());

        // confirm that the tail merging works as expected
        final int mergeResult = numLeadingMatchesAllowed >= 0 ? rtgraph.mergeDanglingHead(result) : rtgraph.mergeDanglingHeadLegacy(result);
        Assert.assertTrue(mergeResult > 0 || !shouldBeMerged);

        // confirm that we created the appropriate bubble in the graph only if expected
        rtgraph.cleanNonRefPaths();
        final SeqGraph seqGraph = rtgraph.toSequenceGraph();
        final List<KBestHaplotype<SeqVertex, BaseEdge>> paths = new GraphBasedKBestHaplotypeFinder<>(seqGraph, seqGraph.getReferenceSourceVertex(), seqGraph.getReferenceSinkVertex()).findBestHaplotypes();
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
            final int matchSize = ReadThreadingGraph.longestSuffixMatch(seq, kmer, i);
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
