package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotype;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.KBestHaplotypeFinder;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.SeqVertex;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public final class ReadThreadingAssemblerUnitTest extends BaseTest {

    private static final boolean DEBUG = false;

    private IndexedFastaSequenceFile seq;
    private SAMFileHeader header;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        seq = new CachingIndexedFastaSequenceFile(new File(hg19_chr1_1M_Reference));
        header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
    }

    @DataProvider(name = "AssembleIntervalsData")
    public Object[][] makeAssembleIntervalsData() {
        List<Object[]> tests = new ArrayList<>();

        final String contig = "1";
        final int start = 100000;
        final int end   = 200001;
        final int windowSize = 100;
        final int stepSize = 200;
        final int nReadsToUse = 5;

        for ( int startI = start; startI < end; startI += stepSize) {
            final int endI = startI + windowSize;
            final SimpleInterval refLoc = new SimpleInterval(contig, startI, endI);
            tests.add(new Object[]{new ReadThreadingAssembler(), refLoc, nReadsToUse});
        }

        return tests.toArray(new Object[][]{});
    }

    @DataProvider(name = "AssembleIntervalsWithVariantData")
    public Object[][] makeAssembleIntervalsWithVariantData() {
        List<Object[]> tests = new ArrayList<>();

        final String contig = "1";
        final int start = 100000;
        final int end   = 101001;
        final int windowSize = 100;
        final int stepSize = 200;
        final int variantStepSize = 1;
        final int nReadsToUse = 5;

        for ( int startI = start; startI < end; startI += stepSize) {
            final int endI = startI + windowSize;
            final SimpleInterval refLoc = new SimpleInterval(contig, startI, endI);
            for ( int variantStart = windowSize / 2 - 10; variantStart < windowSize / 2 + 10; variantStart += variantStepSize ) {
                tests.add(new Object[]{new ReadThreadingAssembler(), refLoc, nReadsToUse, variantStart});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AssembleIntervalsData")
    public void testAssembleRef(final ReadThreadingAssembler assembler, final SimpleInterval loc, final int nReadsToUse) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getEnd()).getBases();

        final List<GATKRead> reads = new LinkedList<>();
        for ( int i = 0; i < nReadsToUse; i++ ) {
            final byte[] bases = refBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, refBases.length);
            final String cigar = refBases.length + "M";
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, loc.getContig(), loc.getContig(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        // TODO -- generalize to all assemblers
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertEquals(haplotypes, Collections.singletonList(refHaplotype));
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndSNP(final ReadThreadingAssembler assembler, final SimpleInterval loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getEnd()).getBases();
        final Allele refBase = Allele.create(refBases[variantSite], true);
        final Allele altBase = Allele.create((byte) (refBase.getBases()[0] == 'A' ? 'C' : 'A'), false);
        final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite, Arrays.asList(refBase, altBase));
        testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndDeletion(final ReadThreadingAssembler assembler, final SimpleInterval loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getEnd()).getBases();
        for ( int deletionLength = 1; deletionLength < 10; deletionLength++ ) {
            final Allele refBase = Allele.create(new String(refBases).substring(variantSite, variantSite + deletionLength + 1), true);
            final Allele altBase = Allele.create(refBase.getBases()[0], false);
            final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite + deletionLength, Arrays.asList(refBase, altBase));
            testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
        }
    }

    @Test(dataProvider = "AssembleIntervalsWithVariantData")
    public void testAssembleRefAndInsertion(final ReadThreadingAssembler assembler, final SimpleInterval loc, final int nReadsToUse, final int variantSite) {
        final byte[] refBases = seq.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getEnd()).getBases();
        for ( int insertionLength = 1; insertionLength < 10; insertionLength++ ) {
            final Allele refBase = Allele.create(refBases[variantSite], false);
            final Allele altBase = Allele.create(new String(refBases).substring(variantSite, variantSite + insertionLength + 1), true);
            final VariantContextBuilder vcb = new VariantContextBuilder("x", loc.getContig(), variantSite, variantSite + insertionLength, Arrays.asList(refBase, altBase));
            testAssemblyWithVariant(assembler, refBases, loc, nReadsToUse, vcb.make());
        }
    }

    private void testAssemblyWithVariant(final ReadThreadingAssembler assembler, final byte[] refBases, final SimpleInterval loc, final int nReadsToUse, final VariantContext site) {
        final String preRef = new String(refBases).substring(0, site.getStart());
        final String postRef = new String(refBases).substring(site.getEnd() + 1, refBases.length);
        final byte[] altBases = (preRef + site.getAlternateAllele(0).getBaseString() + postRef).getBytes();

//        logger.warn("ref " + new String(refBases));
//        logger.warn("alt " + new String(altBases));

        final List<GATKRead> reads = new LinkedList<>();
        for ( int i = 0; i < nReadsToUse; i++ ) {
            final byte[] bases = altBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, altBases.length);
            final String cigar = altBases.length + "M";
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, loc.getContig(), loc.getContig(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Haplotype altHaplotype = new Haplotype(altBases, false);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertEquals(haplotypes, Arrays.asList(refHaplotype, altHaplotype));
    }


    private List<Haplotype> assemble(final ReadThreadingAssembler assembler, final byte[] refBases, final SimpleInterval loc, final List<GATKRead> reads) {
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);

        final AssemblyRegion activeRegion = new AssemblyRegion(loc, null, true, 0, header);
        activeRegion.addAll(reads);
//        logger.warn("Assembling " + activeRegion + " with " + engine);
        final AssemblyResultSet assemblyResultSet =  assembler.runLocalAssembly(activeRegion, refHaplotype, refBases, loc, Collections.<VariantContext>emptyList(), null, header);
        return assemblyResultSet.getHaplotypeList();
    }

    @DataProvider(name = "SimpleAssemblyTestData")
    public Object[][] makeSimpleAssemblyTestData() {
        List<Object[]> tests = new ArrayList<>();

        final String contig  = "1";
        final int start      = 100000;
        final int windowSize = 200;
        final int end        = start + windowSize;

        final int excludeVariantsWithinXbp = 25; // TODO -- decrease to zero when the edge calling problem is fixed

        final String ref = new String(seq.getSubsequenceAt(contig, start, end).getBases());
        final SimpleInterval refLoc = new SimpleInterval(contig, start, end);

        for ( int snpPos = 0; snpPos < windowSize; snpPos++) {
            if ( snpPos > excludeVariantsWithinXbp && (windowSize - snpPos) >= excludeVariantsWithinXbp ) {
                final byte[] altBases = ref.getBytes();
                altBases[snpPos] = altBases[snpPos] == 'A' ? (byte)'C' : (byte)'A';
                final String alt = new String(altBases);
                tests.add(new Object[]{"SNP at " + snpPos, new ReadThreadingAssembler(), refLoc, ref, alt});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SimpleAssemblyTestData")
    public void testSimpleAssembly(final String name, final ReadThreadingAssembler assembler, final SimpleInterval loc, final String ref, final String alt) {
        final byte[] refBases = ref.getBytes();
        final byte[] altBases = alt.getBytes();

        final List<GATKRead> reads = new LinkedList<>();
        for ( int i = 0; i < 20; i++ ) {
            final byte[] bases = altBases.clone();
            final byte[] quals = Utils.dupBytes((byte) 30, altBases.length);
            final String cigar = altBases.length + "M";
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, loc.getContig(), loc.getContig(), loc.getStart(), bases, quals, cigar);
            reads.add(read);
        }

        final Haplotype refHaplotype = new Haplotype(refBases, true);
        final Haplotype altHaplotype = new Haplotype(altBases, false);
        final List<Haplotype> haplotypes = assemble(assembler, refBases, loc, reads);
        Assert.assertTrue(haplotypes.size() > 0, "Failed to find ref haplotype");
        Assert.assertEquals(haplotypes.get(0), refHaplotype);

        Assert.assertEquals(haplotypes.size(), 2, "Failed to find single alt haplotype");
        Assert.assertEquals(haplotypes.get(1), altHaplotype);
    }

    private static class TestAssembler {
        final ReadThreadingAssembler assembler;
        private final SAMFileHeader header;

        Haplotype refHaplotype;
        final List<GATKRead> reads = new LinkedList<>();

        private TestAssembler(final int kmerSize) {
            this.assembler = new ReadThreadingAssembler(100000, Arrays.asList(kmerSize));
            assembler.setJustReturnRawGraph(true);
            assembler.setPruneFactor(0);
            this.header= ArtificialReadUtils.createArtificialSamHeader();
        }

        public void addSequence(final byte[] bases, final boolean isRef) {
            if ( isRef ) {
                refHaplotype = new Haplotype(bases, true);
            } else {
                final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, Utils.dupBytes((byte) 30, bases.length), bases.length + "M");
                reads.add(read);
            }
        }

        public SeqGraph assemble() {
            assembler.setRemovePathsNotConnectedToRef(false); // needed to pass some of the tests
            assembler.setRecoverDanglingBranches(false); // needed to pass some of the tests
            assembler.setDebugGraphTransformations(true);
            assembler.setDebugGraphOutputPath(createTempDir("debugGraphs"));
            final SeqGraph graph = assembler.assemble(reads, refHaplotype, Collections.<Haplotype>emptyList(), header).get(0).getGraph();
            if ( DEBUG ) graph.printGraph(new File("test.dot"), 0);
            return graph;
        }
    }

    private void assertLinearGraph(final TestAssembler assembler, final String seq) {
        final SeqGraph graph = assembler.assemble();
        graph.simplifyGraph();
        Assert.assertEquals(graph.vertexSet().size(), 1);
        Assert.assertEquals(graph.vertexSet().iterator().next().getSequenceString(), seq);
    }

    private void assertSingleBubble(final TestAssembler assembler, final String one, final String two) {
        final SeqGraph graph = assembler.assemble();
        graph.simplifyGraph();
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph);
        Assert.assertEquals(paths.size(), 2);
        final Set<String> expected = new HashSet<>(Arrays.asList(one, two));
        for ( final KBestHaplotype path : paths ) {
            final String seq = new String(path.bases());
            Assert.assertTrue(expected.contains(seq));
            expected.remove(seq);
        }
    }

    @Test
    public void testRefCreation() {
        final String ref = "ACGTAACCGGTT";
        final TestAssembler assembler = new TestAssembler(3);
        assembler.addSequence(ref.getBytes(), true);
        assertLinearGraph(assembler, ref);
    }

    @Test
    public void testRefNonUniqueCreation() {
        final String ref = "GAAAAT";
        final TestAssembler assembler = new TestAssembler(3);
        assembler.addSequence(ref.getBytes(), true);
        assertLinearGraph(assembler, ref);
    }

    @Test
    public void testRefAltCreation() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref = "ACAACTGA";
        final String alt = "ACAGCTGA";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(alt.getBytes(), false);
        assertSingleBubble(assembler, ref, alt);
    }

    @Test
    public void testPartialReadsCreation() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  = "ACAACTGA";
        final String alt1 = "ACAGCT";
        final String alt2 =    "GCTGA";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(alt1.getBytes(), false);
        assembler.addSequence(alt2.getBytes(), false);
        assertSingleBubble(assembler, ref, "ACAGCTGA");
    }

    @Test
    public void testMismatchInFirstKmer() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref = "ACAACTGA";
        final String alt =   "AGCTGA";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(alt.getBytes(), false);

        final SeqGraph graph = assembler.assemble();
        graph.simplifyGraph();
        graph.removeSingletonOrphanVertices();
        final Set<SeqVertex> sources = graph.getSources();
        final Set<SeqVertex> sinks = graph.getSinks();

        Assert.assertEquals(sources.size(), 1);
        Assert.assertEquals(sinks.size(), 1);
        Assert.assertNotNull(graph.getReferenceSourceVertex());
        Assert.assertNotNull(graph.getReferenceSinkVertex());

        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph);
        Assert.assertEquals(paths.size(), 1);
    }

    @Test
    public void testStartInMiddle() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  = "CAAAATG";
        final String read =   "AAATG";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(read.getBytes(), false);
        assertLinearGraph(assembler, ref);
    }

    @Test
    public void testStartInMiddleWithBubble() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  = "CAAAATGGGG";
        final String read =   "AAATCGGG";
        assembler.addSequence(ref.getBytes(), true);
        assembler.addSequence(read.getBytes(), false);
        assertSingleBubble(assembler, ref, "CAAAATCGGG");
    }

    @Test
    public void testCreateWithBasesBeforeRefSource() {
        final TestAssembler assembler = new TestAssembler(3);
        final String ref  =  "ACTG";
        final String read =   "CTGGGACT";
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(ref), true);
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(read), false);
        assertLinearGraph(assembler, "ACTGGGACT");
    }

    @Test
    public void testSingleIndelAsDoubleIndel3Reads() {
        final TestAssembler assembler = new TestAssembler(25);
        // The single indel spans two repetitive structures
        final String ref   = "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCTCTCTGTGTGTGTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        final String read1 = "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCT----------GTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        final String read2 = "GTTTTTCCTAGGCAAATGGTTTCTATAAAATTATGTGTGTGTGTCTCT----------GTGTGTGTGTGTGTGTGTATACCTAATCTCACACTCTTTTTTCTGG";
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(ref), true);
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(read1), false);
        assembler.addSequence(ReadThreadingGraphUnitTest.getBytes(read2), false);

        final SeqGraph graph = assembler.assemble();
        final List<KBestHaplotype> paths = new KBestHaplotypeFinder(graph);
        Assert.assertEquals(paths.size(), 2);
        final byte[] refPath = paths.get(0).bases().length == ref.length() ? paths.get(0).bases() : paths.get(1).bases();
        final byte[] altPath = paths.get(0).bases().length == ref.length() ? paths.get(1).bases() : paths.get(0).bases();
        Assert.assertEquals(refPath, ReadThreadingGraphUnitTest.getBytes(ref));
        Assert.assertEquals(altPath, ReadThreadingGraphUnitTest.getBytes(read1));
    }
}
