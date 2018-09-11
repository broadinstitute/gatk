package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.*;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Stream;
import java.util.stream.Collectors;

public final class AlignmentUtilsUnitTest {
    private final static boolean DEBUG = false;
    private SAMFileHeader header;

    /** Basic aligned and mapped read. */
    private GATKRead readMapped;

    /** This read has a start position, but is flagged that it's not mapped. */
    private GATKRead readUnmappedFlag;

    /** This read says it's aligned, but to a contig not in the header. */
    private GATKRead readUnknownContig;

    @BeforeClass
    public void init() {
        header = ArtificialReadUtils.createArtificialSamHeader(3, 1, ArtificialReadUtils.DEFAULT_READ_LENGTH * 2);

        readMapped = createMappedRead("mapped", 1);

        readUnmappedFlag = createMappedRead("unmappedFlagged", 2);
        readUnmappedFlag.setIsUnmapped();

        readUnknownContig = createMappedRead("unknownContig", 3);
        readUnknownContig.setPosition("unknownContig", 3);
    }

    private GATKRead createMappedRead(String name, int start) {
        return ArtificialReadUtils.createArtificialRead(
                header,
                name,
                0,
                start,
                ArtificialReadUtils.DEFAULT_READ_LENGTH);
    }

    private final List<List<CigarElement>> makeCigarElementCombinations() {
        // this functionality can be adapted to provide input data for whatever you might want in your data
        final List<CigarElement> cigarElements = new LinkedList<>();
        for ( final int size : Arrays.asList(0, 10) ) {
            for ( final CigarOperator op : CigarOperator.values() ) {
                cigarElements.add(new CigarElement(size, op));
            }
        }

        final List<List<CigarElement>> combinations = new LinkedList<>();
        for ( final int nElements : Arrays.asList(1, 2, 3) ) {
            combinations.addAll(Utils.makePermutations(cigarElements, nElements, true));
        }

        return combinations;
    }


    /******************************************************
     * Tests for AlignmentUtils.createReadAlignedToRef()
     ******************************************************/

    private Haplotype makeHaplotypeForAlignedToRefTest(final String bases, final String cigar) {
        final Haplotype hap = new Haplotype(bases.getBytes());
        hap.setCigar(TextCigarCodec.decode(cigar));
        return hap;
    }

    private GATKRead makeReadForAlignedToRefTest(final String baseString) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, 10);
        final byte[] bases = baseString.getBytes();
        read.setBases(bases.clone());
        read.setBaseQualities(Utils.dupBytes((byte)30, read.getLength()));
        return read;
    }

    @DataProvider(name = "ReadAlignedToRefData")
    public Object[][] makeReadAlignedToRefData() {
        List<Object[]> tests = new ArrayList<>();

        final String hapBases = "ACTGAAGGTTCC";
        final Haplotype allM = makeHaplotypeForAlignedToRefTest(hapBases, hapBases.length() + "M");

        // make sure we get back a cigar of the right length
        for ( int i = -1; i < hapBases.length(); i++ ) {
            final GATKRead read = makeReadForAlignedToRefTest(hapBases);
            final byte[] bases = read.getBases();
            if ( i != -1 ) {
                bases[i] = (byte)'A';
            }
            read.setBases(bases);
            tests.add(new Object[]{read, allM, allM, 10, 10, allM.getCigar().toString()});
        }

        // make sure insertions at the front are correctly handled
        for ( int padFront = 1; padFront < 10; padFront++ ) {
            final GATKRead read = makeReadForAlignedToRefTest(Utils.dupString("N", padFront) + hapBases);
            tests.add(new Object[]{read, allM, allM, 10, 10, padFront + "I" + allM.getCigar().toString()});
        }

        // make sure insertions at the back are correctly handled
        for ( int padBack = 1; padBack < 10; padBack++ ) {
            final GATKRead read = makeReadForAlignedToRefTest(hapBases + Utils.dupString("N", padBack));
            tests.add(new Object[]{read, allM, allM, 10, 10, allM.getCigar().toString() + padBack + "I"});
        }

        // make sure refStart and hapStart are respected
        for ( int refStart = 1; refStart < 10; refStart++ ) {
            for ( int hapStart = refStart; hapStart < 10 + refStart; hapStart++ ) {
                final Haplotype hap = new Haplotype(allM.getBases());
                hap.setCigar(allM.getCigar());
                hap.setAlignmentStartHapwrtRef(hapStart);

                final GATKRead read = makeReadForAlignedToRefTest(new String(hap.getBases()));
                tests.add(new Object[]{read, hap, allM, refStart, refStart + hapStart, allM.getCigar().toString()});
            }
        }

        // example case of bad alignment because SW doesn't necessarily left-align indels
        {
            final String hap = "ACTGTGGGTTCCTCTTATTTTATTTCTACATCAATGTTCATATTTAACTTATTATTTTATCTTATTTTTAAATTTCTTTTATGTTGAGCCTTGATGAAAGCCATAGGTTCTCTCATATAATTGTATGTGTATGTATGTATATGTACATAATATATACATATATGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATATATACATATATG";
            final String hapCigar = hap.length() + "M";
            final String readBases = "ATGTACATAATATATACATATATGTATATGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTACATAATATATACGTATATGTATGTGTATGTGTATTACATAATATATACATATATGTATATATTATGTATATGTACATAATAT";
            final GATKRead read = makeReadForAlignedToRefTest(readBases);
            final int refStart = 10130100;
            final int hapStart = 500;
            final String badCigar = "31M6D211M";
            final String goodCigar = "28M6D214M";
            final Haplotype badHap = new Haplotype(hap.getBytes());
            badHap.setCigar(TextCigarCodec.decode(hapCigar));
            badHap.setAlignmentStartHapwrtRef(hapStart);

            final int expectedPos = 10130740;
            tests.add(new Object[]{read, badHap, badHap, refStart, expectedPos, goodCigar});
        }

        // example where left-align generates a leading deletion
        {
            final String ref = "CTGAACGTAACCAAAATCAATATGGATACTGAGAAATACTATTTAATAAAGACATAAATTAGACTGCTAAAAAAAATTAAAGAAATTTCAAAAGAGAATCCACCTCTTTTCCTTGCCAGTGCTCAAAAGTGAGTGTGAATCTGGTGGCTGTGGGGCTGTTTTTGGTGTGGCTCTTTGGACCAGCCTGCCTGGTAATTCAAGCCTGCCTCTCATTTCTG";
            // ref-to-hap:      ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            final String hap = "CTGAACGTAACCAAAATCAATATGGATACTGAGAAATACTATTTAATAAAGACATAAATTAGACTGCTAAAAAAAATTAAAGAAATTTCAAAAGAGAATCCACCTCTTTTCCTTGCCAGTGCTCAAAAGTGAGTGTGAATCTGGTGGCTGCGGGGCTGTTTTTGGTGTGGCTCTTTGGACCAGCCTGCCTGGTAATTCAAGCCTGCCTCTCATTTCTG";
            // hap-to-read:                                                                                                                                                               ||||.|||||||||||||||||
            // aligned read:                                                                                                                                                              GCTGCTTTTGGTGTGGCTCTTT
            final String readBases = "GCTGCTTTTGGTGTGGCTCTTT";
            final String hapCigar = hap.length() + "M";
            final GATKRead read = makeReadForAlignedToRefTest(readBases);
            final int refStart = 215239171;
            final int hapStart = 575;
            final int alignmentOffset = 154;
            final String goodCigar = "22M";
            final Haplotype badHap = new Haplotype(hap.getBytes());
            badHap.setCigar(TextCigarCodec.decode(hapCigar));
            badHap.setAlignmentStartHapwrtRef(hapStart);
            final Haplotype refHap = makeHaplotypeForAlignedToRefTest(ref, ref.length() + "M");

            final int expectedPos = refStart + hapStart + alignmentOffset;
            tests.add(new Object[]{read, badHap, refHap, refStart, expectedPos, goodCigar});
         }

        // example where the haplotype has an indel relative to reference
        {
            final String ref = "GGGATCCTGCTACAAAGGTGAAACCCAGGAGAGTGTGGAGTCCAGAGTGTTGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACAGGGGAATCCCGAAGAAATGGTGGGTCCTGGCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
            // ref-to-hap:      |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||^^||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            final String hap = "GGGATCCTGCTACAAAGGTGAAACCCAGGAGAGTGTGGAGTCCAGAGTGTTGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACGGGAATCCCGAAGAAATGGTGGGTCCTGGCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
            // hap-to-read:                                                                                                                              .|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            // aligned read:                                                                                                                             CCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC
            final String readBases = "CCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC";
            final String hapCigar = "93M2D92M";
            final GATKRead read = makeReadForAlignedToRefTest(readBases);
            final int refStart = 13011;
            final int hapStart = 553;
            final int alignmentOffset = 123;
            final String goodCigar = "64M";
            final Haplotype badHap = new Haplotype(hap.getBytes());
            badHap.setCigar(TextCigarCodec.decode(hapCigar));
            badHap.setAlignmentStartHapwrtRef(hapStart);
            final Haplotype refHap = makeHaplotypeForAlignedToRefTest(ref, ref.length() + "M");

            final int expectedPos = refStart + hapStart + alignmentOffset;
            tests.add(new Object[]{read, badHap, refHap, refStart, expectedPos, goodCigar});
         }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadAlignedToRefData")
    public void testReadAlignedToRef(final GATKRead read, final Haplotype haplotype, final Haplotype refHaplotype, final int refStart, final int expectedReadStart, final String expectedReadCigar) throws Exception {
        final GATKRead originalReadCopy = read.copy();

        if ( expectedReadCigar == null ) {
            Assert.assertNull(AlignmentUtils.createReadAlignedToRef(read, haplotype, refHaplotype, refStart, true, SmithWatermanJavaAligner
                    .getInstance()));
        } else {
            final Cigar expectedCigar = TextCigarCodec.decode(expectedReadCigar);
            final GATKRead alignedRead = AlignmentUtils.createReadAlignedToRef(read, haplotype, refHaplotype, refStart, true, SmithWatermanJavaAligner
                    .getInstance());

            Assert.assertEquals(alignedRead.getName(), originalReadCopy.getName());
            Assert.assertEquals(alignedRead.getStart(), expectedReadStart);
            Assert.assertEquals(alignedRead.getBases(), originalReadCopy.getBases());
            Assert.assertEquals(alignedRead.getBaseQualities(), originalReadCopy.getBaseQualities());
            Assert.assertEquals(alignedRead.getStart(), expectedReadStart);
            Assert.assertEquals(alignedRead.getCigar(), expectedCigar);
            Assert.assertNotNull(alignedRead.getAttributeAsInteger("HC"));
        }

        Assert.assertEquals(read, originalReadCopy, "createReadAlignedToRef seems be modifying the original read!");
    }

    private static class Mutation implements Comparable<Mutation> {
        int pos, len;
        CigarOperator operator;

        private Mutation(int pos, int len, CigarOperator operator) {
            this.pos = pos;
            this.len = len;
            this.operator = operator;
        }
        public int getNMismatches() { return len; }

        @Override
        public int compareTo(Mutation o) {
            return Integer.valueOf(pos).compareTo(o.pos);
        }

        private String apply(final String seq) {
            switch ( operator ) {
                case M:
                    final byte[] bases = seq.getBytes();
                    if ( pos < seq.length() )
                        bases[pos] = (byte)(bases[pos] == 'A' ? 'C' : 'A');
                    return new String(bases);
                case I: {
                    final String prefix = seq.substring(0, pos);
                    final String postfix = seq.substring(pos, seq.length());
                    return prefix + "GTCAGTTA".substring(0, len) + postfix;
                } case D: {
                    final String prefix = seq.substring(0, pos);
                    final String postfix = seq.substring(pos + len, seq.length());
                    return prefix + postfix;
                }default:
                    throw new IllegalStateException("Unexpected operator " + operator);
            }
        }
    }

    private static class MutatedSequence {
        int numMismatches;
        String seq;

        private MutatedSequence(int numMismatches, String seq) {
            this.numMismatches = numMismatches;
            this.seq = seq;
        }
    }

    private MutatedSequence mutateSequence(final String hapIn, final List<Mutation> mutations) {
        Collections.sort(mutations);
        int mismatches = 0;
        String hap = hapIn;
        for ( final Mutation mut : mutations ) {
            hap = mut.apply(hap);
            mismatches += mut.getNMismatches();
        }
        return new MutatedSequence(mismatches, hap);
    }

    @DataProvider(name = "ComplexReadAlignedToRef")
    public Object[][] makeComplexReadAlignedToRef() {
        List<Object[]> tests = new ArrayList<>();

        final List<Mutation> allMutations = Arrays.asList(
                new Mutation(1, 1, CigarOperator.M),
                new Mutation(2, 1, CigarOperator.M),
                new Mutation(3, 1, CigarOperator.I),
                new Mutation(7, 1, CigarOperator.D)
        );

        int i = 0;
        final String referenceBases  = "ACTGACTGACTG";
        final String paddedReference = "NNNN" + referenceBases + "NNNN";
        for ( final List<Mutation> mutations : Utils.makePermutations(allMutations, 3, false) ) {
            final MutatedSequence hap = mutateSequence(referenceBases, mutations);
            final Haplotype haplotype = new Haplotype(hap.seq.getBytes());
            final SmithWatermanAlignment align = SmithWatermanJavaAligner.getInstance()
                    .align(paddedReference.getBytes(), hap.seq.getBytes(), SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.SOFTCLIP);
            haplotype.setAlignmentStartHapwrtRef(align.getAlignmentOffset());
            haplotype.setCigar(align.getCigar());

            for ( final List<Mutation> readMutations : Utils.makePermutations(allMutations, 3, false) ) {
                final MutatedSequence readBases = mutateSequence(hap.seq, readMutations);
                final GATKRead read = makeReadForAlignedToRefTest(readBases.seq);
                tests.add(new Object[]{i++, read, paddedReference, haplotype, hap.numMismatches + readBases.numMismatches});
            }
        }

        // for convenient testing of a single failing case
        //tests.add(new Object[]{makeRead("ACCGGGACTGACTG"), reference, makeHaplotype("AAAGGACTGACTG", "1M1I11M"), 2});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ComplexReadAlignedToRef")
    public void testReadAlignedToRefComplexAlignment(final int testIndex, final GATKRead read, final String reference, final Haplotype haplotype, final int expectedMaxMismatches) throws Exception {
        final GATKRead alignedRead = AlignmentUtils.createReadAlignedToRef(read, haplotype, new Haplotype(reference.getBytes(),true), 1, true, SmithWatermanJavaAligner
                .getInstance());
        if ( alignedRead != null ) {
            final int mismatches = AlignmentUtils.getMismatchCount(alignedRead, reference.getBytes(), alignedRead.getStart() - 1).numMismatches;
            Assert.assertTrue(mismatches <= expectedMaxMismatches,
                    "Alignment of read to ref looks broken.  Expected at most " + expectedMaxMismatches + " but saw " + mismatches
                            + " for readBases " + new String(read.getBases()) + " with cigar " + read.getCigar() + " reference " + reference + " haplotype "
                            + haplotype + " with cigar " + haplotype.getCigar() + " aligned read cigar " + alignedRead.getCigar().toString() + " @ " + alignedRead.getStart());
        }
    }

    /**********************************************************
     * End of Tests for AlignmentUtils.createReadAlignedToRef()
     **********************************************************/


    @DataProvider(name = "CalcNumDifferentBasesData")
    public Object[][] makeCalcNumDifferentBasesData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{"5M", "ACGTA", "ACGTA", 0});
        tests.add(new Object[]{"5M", "ACGTA", "ACGTT", 1});
        tests.add(new Object[]{"5M", "ACGTA", "TCGTT", 2});
        tests.add(new Object[]{"5M", "ACGTA", "TTGTT", 3});
        tests.add(new Object[]{"5M", "ACGTA", "TTTTT", 4});
        tests.add(new Object[]{"5M", "ACGTA", "TTTCT", 5});
        tests.add(new Object[]{"2M3I3M", "ACGTA", "ACNNNGTA", 3});
        tests.add(new Object[]{"2M3I3M", "ACGTA", "ACNNNGTT", 4});
        tests.add(new Object[]{"2M3I3M", "ACGTA", "TCNNNGTT", 5});
        tests.add(new Object[]{"2M2D1M", "ACGTA", "ACA", 2});
        tests.add(new Object[]{"2M2D1M", "ACGTA", "ACT", 3});
        tests.add(new Object[]{"2M2D1M", "ACGTA", "TCT", 4});
        tests.add(new Object[]{"2M2D1M", "ACGTA", "TGT", 5});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CalcNumDifferentBasesData")
    public void testCalcNumDifferentBases(final String cigarString, final String ref, final String read, final int expectedDifferences) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        Assert.assertEquals(AlignmentUtils.calcNumDifferentBases(cigar, ref.getBytes(), read.getBytes()), expectedDifferences);
    }

    @DataProvider(name = "NumAlignedBasesCountingSoftClips")
    public Object[][] makeNumAlignedBasesCountingSoftClips() {
        List<Object[]> tests = new ArrayList<>();

        final EnumSet<CigarOperator> alignedToGenome = EnumSet.of(CigarOperator.M, CigarOperator.EQ, CigarOperator.X, CigarOperator.S);
        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            int n = 0;
            for ( final CigarElement elt : elements ) n += alignedToGenome.contains(elt.getOperator()) ? elt.getLength() : 0;
            tests.add(new Object[]{new Cigar(elements), n});
        }

        tests.add(new Object[]{null, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "NumAlignedBasesCountingSoftClips")
    public void testNumAlignedBasesCountingSoftClips(final Cigar cigar, final int expected) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, cigar == null ? 10 : cigar.getReadLength());
        read.setCigar(cigar);
        Assert.assertEquals(AlignmentUtils.getNumAlignedBasesCountingSoftClips(read), expected, "Cigar " + cigar + " failed NumAlignedBasesCountingSoftClips");
    }

    @DataProvider(name = "CigarHasZeroElement")
    public Object[][] makeCigarHasZeroElement() {
        List<Object[]> tests = new ArrayList<>();

        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            boolean hasZero = false;
            for ( final CigarElement elt : elements ) hasZero = hasZero || elt.getLength() == 0;
            tests.add(new Object[]{new Cigar(elements), hasZero});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "CigarHasZeroElement")
    public void testCigarHasZeroSize(final Cigar cigar, final boolean hasZero) {
        Assert.assertEquals(AlignmentUtils.cigarHasZeroSizeElement(cigar), hasZero, "Cigar " + cigar.toString() + " failed cigarHasZeroSizeElement");
    }

    @DataProvider(name = "NumHardClipped")
    public Object[][] makeNumHardClipped() {
        List<Object[]> tests = new ArrayList<>();

        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            int n = 0;
            for ( final CigarElement elt : elements ) n += elt.getOperator() == CigarOperator.H ? elt.getLength() : 0;
            tests.add(new Object[]{new Cigar(elements), n});
        }

        tests.add(new Object[]{null, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "NumHardClipped")
    public void testNumHardClipped(final Cigar cigar, final int expected) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, cigar == null ? 10 : cigar.getReadLength());
        read.setCigar(cigar);
        Assert.assertEquals(AlignmentUtils.getNumHardClippedBases(read), expected, "Cigar " + cigar + " failed num hard clips");
    }

    @DataProvider(name = "NumAlignedBlocks")
    public Object[][] makeNumAlignedBlocks() {
        List<Object[]> tests = new ArrayList<>();

        for ( final List<CigarElement> elements : makeCigarElementCombinations() ) {
            int n = 0;
            for ( final CigarElement elt : elements ) {
                switch ( elt.getOperator() ) {
                    case M:case X:case EQ: n++; break;
                    default: break;
                }
            }
            tests.add(new Object[]{new Cigar(elements), n});
        }

        tests.add(new Object[]{null, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "NumAlignedBlocks")
    public void testNumAlignedBlocks(final Cigar cigar, final int expected) {
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, cigar == null ? 10 : cigar.getReadLength());
        read.setCigar(cigar);
        Assert.assertEquals(AlignmentUtils.getNumAlignmentBlocks(read), expected, "Cigar " + cigar + " failed NumAlignedBlocks");
    }

    @DataProvider(name = "ConsolidateCigarData")
    public Object[][] makeConsolidateCigarData() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{"1M1M", "2M"});
        tests.add(new Object[]{"2M", "2M"});
        tests.add(new Object[]{"2M0M", "2M"});
        tests.add(new Object[]{"0M2M", "2M"});
        tests.add(new Object[]{"0M2M0M0I0M1M", "3M"});
        tests.add(new Object[]{"2M0M1M", "3M"});
        tests.add(new Object[]{"1M1M1M1D2M1M", "3M1D3M"});
        tests.add(new Object[]{"6M6M6M", "18M"});

        final List<CigarElement> elements = new LinkedList<>();
        int i = 1;
        for ( final CigarOperator op : CigarOperator.values() ) {
            elements.add(new CigarElement(i++, op));
        }
        for ( final List<CigarElement> ops : Utils.makePermutations(elements,  3, false) ) {
            final String expected = new Cigar(ops).toString();
            final List<CigarElement> cutElements = new LinkedList<>();
            for ( final CigarElement elt : ops ) {
                for ( int j = 0; j < elt.getLength(); j++ ) {
                    cutElements.add(new CigarElement(1, elt.getOperator()));
                }
            }

            final String actual = new Cigar(cutElements).toString();
            tests.add(new Object[]{actual, expected});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "ConsolidateCigarData")
    public void testConsolidateCigarWithData(final String testCigarString, final String expectedCigarString) {
        final Cigar testCigar = TextCigarCodec.decode(testCigarString);
        final Cigar expectedCigar = TextCigarCodec.decode(expectedCigarString);
        final Cigar actualCigar = AlignmentUtils.consolidateCigar(testCigar);
        Assert.assertEquals(actualCigar, expectedCigar);
    }

    @DataProvider(name = "SoftClipsDataProvider")
    public Object[][] makeSoftClipsDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for ( final int lengthOfLeftClip : Arrays.asList(0, 1, 10) ) {
            for ( final int lengthOfRightClip : Arrays.asList(0, 1, 10) ) {
                for ( final int qualThres : Arrays.asList(10, 20, 30) ) {
                    for ( final String middleOp : Arrays.asList("M", "D") ) {
                        for ( final int matchSize : Arrays.asList(0, 1, 10) ) {
                            final byte[] left = makeQualArray(lengthOfLeftClip, qualThres);
                            final byte[] right = makeQualArray(lengthOfRightClip, qualThres);
                            int n = 0;
                            for ( int i = 0; i < left.length; i++ ) n += left[i] > qualThres ? 1 : 0;
                            for ( int i = 0; i < right.length; i++ ) n += right[i] > qualThres ? 1 : 0;
                            tests.add(new Object[]{left, matchSize, middleOp, right, qualThres, n});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private byte[] makeQualArray(final int length, final int qualThreshold) {
        final byte[] array = new byte[length];
        for ( int i = 0; i < array.length; i++ )
            array[i] = (byte)(qualThreshold + ( i % 2 == 0 ? 1 : - 1 ));
        return array;
    }

    @Test(enabled = !DEBUG, dataProvider = "SoftClipsDataProvider")
    public void testSoftClipsData(final byte[] qualsOfSoftClipsOnLeft, final int middleSize, final String middleOp, final byte[] qualOfSoftClipsOnRight, final int qualThreshold, final int numExpected) {
        final int readLength = (middleOp.equals("D") ? 0 : middleSize) + qualOfSoftClipsOnRight.length + qualsOfSoftClipsOnLeft.length;
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
        final byte[] bases = Utils.dupBytes((byte) 'A', readLength);
        final byte[] matchBytes = middleOp.equals("D") ? new byte[]{} : Utils.dupBytes((byte)30, middleSize);
        final byte[] quals = ArrayUtils.addAll(ArrayUtils.addAll(qualsOfSoftClipsOnLeft, matchBytes), qualOfSoftClipsOnRight);

        // set the read's bases and quals
        read.setBases(bases);
        read.setBaseQualities(quals);

        final StringBuilder cigar = new StringBuilder();
        if (qualsOfSoftClipsOnLeft.length > 0 ) cigar.append(qualsOfSoftClipsOnLeft.length + "S");
        if (middleSize > 0 ) cigar.append(middleSize + middleOp);
        if (qualOfSoftClipsOnRight.length > 0 ) cigar.append(qualOfSoftClipsOnRight.length + "S");

        read.setCigar(cigar.toString());

        final int actual = AlignmentUtils.calcNumHighQualitySoftClips(read, (byte) qualThreshold);
        Assert.assertEquals(actual, numExpected, "Wrong number of soft clips detected for read " + read.toString());
    }

    ////////////////////////////////////////////
    // Test AlignmentUtils.getMismatchCount() //
    ////////////////////////////////////////////

    @DataProvider(name = "MismatchCountDataProvider")
    public Object[][] makeMismatchCountDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final int readLength = 20;
        final int lengthOfIndel = 2;
        final int locationOnReference = 10;
        final byte[] reference = Utils.dupBytes((byte)'A', readLength);
        final byte[] quals = Utils.dupBytes((byte)'A', readLength);


        for ( int startOnRead = 0; startOnRead <= readLength; startOnRead++ ) {
            for ( int basesToRead = 0; basesToRead <= readLength; basesToRead++ ) {
                for ( final int lengthOfSoftClip : Arrays.asList(0, 1, 10) ) {
                    for ( final int lengthOfFirstM : Arrays.asList(0, 3) ) {
                        for ( final char middleOp : Arrays.asList('M', 'D', 'I') ) {
                            for ( final int mismatchLocation : Arrays.asList(-1, 0, 5, 10, 15, 19) ) {

                                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, locationOnReference, readLength);

                                // set the read's bases and quals
                                final byte[] readBases = Arrays.copyOf(reference, reference.length);
                                // create the mismatch if requested
                                if ( mismatchLocation != -1 )
                                    readBases[mismatchLocation] = (byte)'C';
                                read.setBases(readBases);
                                read.setBaseQualities(quals);

                                // create the CIGAR string
                                read.setCigar(buildTestCigarString(middleOp, lengthOfSoftClip, lengthOfFirstM, lengthOfIndel, readLength));

                                // now, determine whether or not there's a mismatch
                                final boolean isMismatch;
                                if ( mismatchLocation < startOnRead || mismatchLocation >= startOnRead + basesToRead || mismatchLocation < lengthOfSoftClip ) {
                                    isMismatch = false;
                                } else if ( middleOp == 'M' || middleOp == 'D' || mismatchLocation < lengthOfSoftClip + lengthOfFirstM || mismatchLocation >= lengthOfSoftClip + lengthOfFirstM + lengthOfIndel ) {
                                    isMismatch = true;
                                } else {
                                    isMismatch = false;
                                }

                                tests.add(new Object[]{read, locationOnReference, startOnRead, basesToRead, isMismatch});
                            }
                        }
                    }
                }
            }
        }

        // Adding test to make sure soft-clipped reads go through the exceptions thrown at the beginning of the getMismatchCount method
        // todo: incorporate cigars with right-tail soft-clips in the systematic tests above.
        GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 10, 20);
        read.setBases(reference);
        read.setBaseQualities(quals);
        read.setCigar("10S5M5S");
        tests.add(new Object[]{read, 10, read.getStart(), read.getLength(), false});

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "MismatchCountDataProvider")
    public void testMismatchCountData(final GATKRead read, final int refIndex, final int startOnRead, final int basesToRead, final boolean isMismatch) {
        final byte[] reference = Utils.dupBytes((byte)'A', 100);
        final int actual = AlignmentUtils.getMismatchCount(read, reference, refIndex, startOnRead, basesToRead).numMismatches;
        Assert.assertEquals(actual, isMismatch ? 1 : 0, "Wrong number of mismatches detected for read " + read.toString());
    }

    private static String buildTestCigarString(final char middleOp, final int lengthOfSoftClip, final int lengthOfFirstM, final int lengthOfIndel, final int readLength) {
        final StringBuilder cigar = new StringBuilder();
        int remainingLength = readLength;

        // add soft clips to the beginning of the read
        if (lengthOfSoftClip > 0 ) {
            cigar.append(lengthOfSoftClip).append("S");
            remainingLength -= lengthOfSoftClip;
        }

        if ( middleOp == 'M' ) {
            cigar.append(remainingLength).append("M");
        } else {
            if ( lengthOfFirstM > 0 ) {
                cigar.append(lengthOfFirstM).append("M");
                remainingLength -= lengthOfFirstM;
            }

            if ( middleOp == 'D' ) {
                cigar.append(lengthOfIndel).append("D");
            } else {
                cigar.append(lengthOfIndel).append("I");
                remainingLength -= lengthOfIndel;
            }
            cigar.append(remainingLength).append("M");
        }

        return cigar.toString();
    }

    ////////////////////////////////////////////////////////
    // Test AlignmentUtils.calcAlignmentByteArrayOffset() //
    ////////////////////////////////////////////////////////

    @DataProvider(name = "AlignmentByteArrayOffsetDataProvider")
    public Object[][] makeAlignmentByteArrayOffsetDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final int readLength = 20;
        final int lengthOfIndel = 2;
        final int locationOnReference = 20;

        for ( int offset = 0; offset < readLength; offset++ ) {
            for ( final int lengthOfSoftClip : Arrays.asList(0, 1, 10) ) {
                for ( final int lengthOfFirstM : Arrays.asList(0, 3) ) {
                    for ( final char middleOp : Arrays.asList('M', 'D', 'I') ) {

                        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, locationOnReference, readLength);
                        // create the CIGAR string
                        read.setCigar(buildTestCigarString(middleOp, lengthOfSoftClip, lengthOfFirstM, lengthOfIndel, readLength));

                        // now, determine the expected alignment offset
                        final int expected;
                        boolean isDeletion = false;
                        if ( offset < lengthOfSoftClip ) {
                            expected = 0;
                        } else if ( middleOp == 'M' || offset < lengthOfSoftClip + lengthOfFirstM ) {
                            expected = offset - lengthOfSoftClip;
                        } else if ( offset < lengthOfSoftClip + lengthOfFirstM + lengthOfIndel ) {
                            if ( middleOp == 'D' ) {
                                isDeletion = true;
                                expected = offset - lengthOfSoftClip;
                            } else {
                                expected = lengthOfFirstM;
                            }
                        } else {
                            expected = offset - lengthOfSoftClip - (middleOp == 'I' ? lengthOfIndel : -lengthOfIndel);
                        }

                        tests.add(new Object[]{read.getCigar(), offset, expected, isDeletion, lengthOfSoftClip});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "AlignmentByteArrayOffsetDataProvider")
    public void testAlignmentByteArrayOffsetData(final Cigar cigar, final int offset, final int expectedResult, final boolean isDeletion, final int lengthOfSoftClip) {
        final int actual = AlignmentUtils.calcAlignmentByteArrayOffset(cigar, isDeletion ? -1 : offset, isDeletion, 20, 20 + offset - lengthOfSoftClip);
        Assert.assertEquals(actual, expectedResult, "Wrong alignment offset detected for cigar " + cigar.toString());
    }

    @Test
    public void testIsInsideDeletion() {
        final List<CigarElement> cigarElements = Arrays.asList(new CigarElement(5, CigarOperator.S),
                new CigarElement(5, CigarOperator.M),
                new CigarElement(5, CigarOperator.EQ),
                new CigarElement(6, CigarOperator.N),
                new CigarElement(5, CigarOperator.X),
                new CigarElement(6, CigarOperator.D),
                new CigarElement(1, CigarOperator.P),
                new CigarElement(1, CigarOperator.H));
        final Cigar cigar = new Cigar(cigarElements);
        for ( int i=-1; i <= 20; i++ ) {
            Assert.assertFalse(AlignmentUtils.isInsideDeletion(cigar, i));
        }
        for ( int i=21; i <= 26; i++ ){
            Assert.assertTrue(AlignmentUtils.isInsideDeletion(cigar, i));
        }
        for ( int i=27; i <= 28; i++ ) {
            Assert.assertFalse(AlignmentUtils.isInsideDeletion(cigar, i));
        }
    }

    ////////////////////////////////////////////////////
    // Test AlignmentUtils.readToAlignmentByteArray() //
    ////////////////////////////////////////////////////

    @DataProvider(name = "ReadToAlignmentByteArrayDataProvider")
    public Object[][] makeReadToAlignmentByteArrayDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final int readLength = 20;
        final int lengthOfIndel = 2;
        final int locationOnReference = 20;

        for ( final int lengthOfSoftClip : Arrays.asList(0, 1, 10) ) {
            for ( final int lengthOfFirstM : Arrays.asList(0, 3) ) {
                for ( final char middleOp : Arrays.asList('M', 'D', 'I') ) {

                    final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, locationOnReference, readLength);
                    // create the CIGAR string
                    read.setCigar(buildTestCigarString(middleOp, lengthOfSoftClip, lengthOfFirstM, lengthOfIndel, readLength));

                    // now, determine the byte array size
                    final int expected = readLength - lengthOfSoftClip - (middleOp == 'I' ? lengthOfIndel : (middleOp == 'D' ? -lengthOfIndel : 0));
                    final int indelBasesStart = middleOp != 'M' ? lengthOfFirstM : -1;

                    tests.add(new Object[]{read.getCigar(), expected, middleOp, indelBasesStart, lengthOfIndel});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "ReadToAlignmentByteArrayDataProvider")
    public void testReadToAlignmentByteArrayData(final Cigar cigar, final int expectedLength, final char middleOp, final int startOfIndelBases, final int lengthOfDeletion) {
        final byte[] read = Utils.dupBytes((byte)'A', cigar.getReadLength());
        final byte[] alignment = AlignmentUtils.readToAlignmentByteArray(cigar, read);

        Assert.assertEquals(alignment.length, expectedLength, "Wrong alignment length detected for cigar " + cigar.toString());

        for ( int i = 0; i < alignment.length; i++ ) {
            final byte expectedBase;
            if ( middleOp == 'D' && i >= startOfIndelBases && i < startOfIndelBases + lengthOfDeletion )
                expectedBase = PileupElement.DELETION_BASE;
            else if ( middleOp == 'I' && i == startOfIndelBases - 1 )
                expectedBase = PileupElement.A_FOLLOWED_BY_INSERTION_BASE;
            else
                expectedBase = (byte)'A';
            Assert.assertEquals(alignment[i], expectedBase, "Wrong base detected at position " + i);
        }
    }

    //////////////////////////////////////////
    // Test AlignmentUtils.leftAlignIndel() //
    //////////////////////////////////////////



    @DataProvider(name = "LeftAlignIndelDataProvider")
    public Object[][] makeLeftAlignIndelDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final byte[] repeat1Reference = "ABCDEFGHIJKLMNOPXXXXXXXXXXABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat2Reference = "ABCDEFGHIJKLMNOPXYXYXYXYXYABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat3Reference = "ABCDEFGHIJKLMNOPXYZXYZXYZXYZABCDEFGHIJKLMN".getBytes();
        final int referenceLength = repeat1Reference.length;

        for ( int indelStart = 0; indelStart < repeat1Reference.length; indelStart++ ) {
            for ( final int indelSize : Arrays.asList(0, 1, 2, 3, 4) ) {
                for ( final char indelOp : Arrays.asList('D', 'I') ) {

                    if ( indelOp == 'D' && indelStart + indelSize >= repeat1Reference.length )
                        continue;

                    final int readLength = referenceLength - (indelOp == 'D' ? indelSize : -indelSize);

                    // create the original CIGAR string
                    final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
                    read.setCigar(buildTestCigarString(indelSize == 0 ? 'M' : indelOp, 0, indelStart, indelSize, readLength));
                    final Cigar originalCigar = read.getCigar();

                    final Cigar expectedCigar1 = makeExpectedCigar1(originalCigar, indelOp, indelStart, indelSize, readLength);
                    final byte[] readString1 = makeReadString(repeat1Reference, indelOp, indelStart, indelSize, readLength, 1);
                    tests.add(new Object[]{originalCigar, expectedCigar1, repeat1Reference, readString1, 1});

                    final Cigar expectedCigar2 = makeExpectedCigar2(originalCigar, indelOp, indelStart, indelSize, readLength);
                    final byte[] readString2 = makeReadString(repeat2Reference, indelOp, indelStart, indelSize, readLength, 2);
                    tests.add(new Object[]{originalCigar, expectedCigar2, repeat2Reference, readString2, 2});

                    final Cigar expectedCigar3 = makeExpectedCigar3(originalCigar, indelOp, indelStart, indelSize, readLength);
                    final byte[] readString3 = makeReadString(repeat3Reference, indelOp, indelStart, indelSize, readLength, 3);
                    tests.add(new Object[]{originalCigar, expectedCigar3, repeat3Reference, readString3, 3});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private Cigar makeExpectedCigar1(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength) {
        if ( indelSize == 0 || indelStart < 17 || indelStart > (26 - (indelOp == 'D' ? indelSize : 0)) )
            return originalCigar;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
        read.setCigar(buildTestCigarString(indelOp, 0, 16, indelSize, readLength));
        return read.getCigar();
    }

    private Cigar makeExpectedCigar2(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength) {
        if ( indelStart < 17 || indelStart > (26 - (indelOp == 'D' ? indelSize : 0)) )
            return originalCigar;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);

        if ( indelOp == 'I' && (indelSize == 1 || indelSize == 3) && indelStart % 2 == 1 )
            read.setCigar(buildTestCigarString(indelOp, 0, Math.max(indelStart - indelSize, 16), indelSize, readLength));
        else if ( (indelSize == 2 || indelSize == 4) && (indelOp == 'D' || indelStart % 2 == 0) )
            read.setCigar(buildTestCigarString(indelOp, 0, 16, indelSize, readLength));
        else
            return originalCigar;

        return read.getCigar();
    }

    private Cigar makeExpectedCigar3(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength) {
        if ( indelStart < 17 || indelStart > (28 - (indelOp == 'D' ? indelSize : 0)) )
            return originalCigar;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);

        if ( indelSize == 3 && (indelOp == 'D' || indelStart % 3 == 1) )
            read.setCigar(buildTestCigarString(indelOp, 0, 16, indelSize, readLength));
        else if ( (indelOp == 'I' && indelSize == 4 && indelStart % 3 == 2) ||
                (indelOp == 'I' && indelSize == 2 && indelStart % 3 == 0) ||
                (indelOp == 'I' && indelSize == 1 && indelStart < 28 && indelStart % 3 == 2) )
            read.setCigar(buildTestCigarString(indelOp, 0, Math.max(indelStart - indelSize, 16), indelSize, readLength));
        else
            return originalCigar;

        return read.getCigar();
    }

    private static byte[] makeReadString(final byte[] reference, final char indelOp, final int indelStart, final int indelSize, final int readLength, final int repeatLength) {
        final byte[] readString = new byte[readLength];

        if ( indelOp == 'D' && indelSize > 0 ) {
            System.arraycopy(reference, 0, readString, 0, indelStart);
            System.arraycopy(reference, indelStart + indelSize, readString, indelStart, readLength - indelStart);
        } else if ( indelOp == 'I' && indelSize > 0 ) {
            System.arraycopy(reference, 0, readString, 0, indelStart);
            for ( int i = 0; i < indelSize; i++ ) {
                if ( i % repeatLength == 0 )
                    readString[indelStart + i] = 'X';
                else if ( i % repeatLength == 1 )
                    readString[indelStart + i] = 'Y';
                else
                    readString[indelStart + i] = 'Z';
            }
            System.arraycopy(reference, indelStart, readString, indelStart + indelSize, readLength - indelStart - indelSize);
        } else {
            System.arraycopy(reference, 0, readString, 0, readLength);
        }

        return readString;
    }

    @Test(enabled = !DEBUG, dataProvider = "LeftAlignIndelDataProvider")
    public void testLeftAlignIndelData(final Cigar originalCigar, final Cigar expectedCigar, final byte[] reference, final byte[] read, final int repeatLength) {
        final Cigar actualCigar = AlignmentUtils.leftAlignIndel(originalCigar, reference, read, 0, 0, true);
        Assert.assertTrue(expectedCigar.equals(actualCigar), "Wrong left alignment detected for cigar " + originalCigar.toString() + " to " + actualCigar.toString() + " but expected " + expectedCigar.toString() + " with repeat length " + repeatLength);
    }

    //additionally test that will not align past leftmostAllowedAlignment.  Dataset is same as provided by LeftAlignIndelDataProvider
    // with additional parameter of leftmostAllowedAlignment, and adjusted expected results
    @DataProvider(name = "LeftAlignIndelWithLimitDataProvider")
    public Object[][] makeLeftAlignIndelWithLimitDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final byte[] repeat1Reference = "ABCDEFGHIJKLMNOPXXXXXXXXXXABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat2Reference = "ABCDEFGHIJKLMNOPXYXYXYXYXYABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat3Reference = "ABCDEFGHIJKLMNOPXYZXYZXYZXYZABCDEFGHIJKLMN".getBytes();
        final int referenceLength = repeat1Reference.length;

        for (int indelStart = 0; indelStart < repeat1Reference.length; indelStart++) {
            for (final int indelSize : Arrays.asList(0, 1, 2, 3, 4)) {
                for (final char indelOp : Arrays.asList('D', 'I')) {
                    for (int leftmostAllowedAlignment = 0; leftmostAllowedAlignment < indelStart; leftmostAllowedAlignment++) {
                        if (indelOp == 'D' && indelStart + indelSize >= repeat1Reference.length)
                            continue;

                        final int readLength = referenceLength - (indelOp == 'D' ? indelSize : -indelSize);

                        // create the original CIGAR string
                        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
                        read.setCigar(buildTestCigarString(indelSize == 0 ? 'M' : indelOp, 0, indelStart, indelSize, readLength));
                        final Cigar originalCigar = read.getCigar();

                        final Cigar expectedCigar1 = makeExpectedCigar1WithLimit(originalCigar, indelOp, indelStart, indelSize, readLength, leftmostAllowedAlignment);
                        final byte[] readString1 = makeReadString(repeat1Reference, indelOp, indelStart, indelSize, readLength, 1);
                        tests.add(new Object[]{originalCigar, expectedCigar1, repeat1Reference, readString1, 1, leftmostAllowedAlignment});

                        final Cigar expectedCigar2 = makeExpectedCigar2WithLimit(originalCigar, indelOp, indelStart, indelSize, readLength, leftmostAllowedAlignment);
                        final byte[] readString2 = makeReadString(repeat2Reference, indelOp, indelStart, indelSize, readLength, 2);
                        tests.add(new Object[]{originalCigar, expectedCigar2, repeat2Reference, readString2, 2, leftmostAllowedAlignment});

                        final Cigar expectedCigar3 = makeExpectedCigar3WithLimit(originalCigar, indelOp, indelStart, indelSize, readLength, leftmostAllowedAlignment);
                        final byte[] readString3 = makeReadString(repeat3Reference, indelOp, indelStart, indelSize, readLength, 3);
                        tests.add(new Object[]{originalCigar, expectedCigar3, repeat3Reference, readString3, 3, leftmostAllowedAlignment});
                    }
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    private Cigar makeExpectedCigar1WithLimit(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength, final int leftmostAllowedAlignment) {
        if (indelSize == 0 || indelStart < 17 || indelStart > (26 - (indelOp == 'D' ? indelSize : 0)) || indelStart <= leftmostAllowedAlignment)
            return originalCigar;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
        read.setCigar(buildTestCigarString(indelOp, 0, Math.max(16, leftmostAllowedAlignment), indelSize, readLength));
        return read.getCigar();
    }

    private Cigar makeExpectedCigar2WithLimit(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength, final int leftmostAllowedAlignment) {
        if (indelStart < 17 || indelStart > (26 - (indelOp == 'D' ? indelSize : 0)) || indelStart <= leftmostAllowedAlignment)
            return originalCigar;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);

        if (indelOp == 'I' && (indelSize == 1 || indelSize == 3) && indelStart % 2 == 1)
            read.setCigar(buildTestCigarString(indelOp, 0, Math.max(Math.max(indelStart - indelSize, 16), leftmostAllowedAlignment), indelSize, readLength));
        else if ((indelSize == 2 || indelSize == 4) && (indelOp == 'D' || indelStart % 2 == 0))
            read.setCigar(buildTestCigarString(indelOp, 0, Math.max(16, leftmostAllowedAlignment), indelSize, readLength));
        else
            return originalCigar;

        return read.getCigar();
    }

    private Cigar makeExpectedCigar3WithLimit(final Cigar originalCigar, final char indelOp, final int indelStart, final int indelSize, final int readLength, final int leftmostAllowedAlignment) {
        if (indelStart < 17 || indelStart > (28 - (indelOp == 'D' ? indelSize : 0)) || indelStart <= leftmostAllowedAlignment)
            return originalCigar;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);

        if (indelSize == 3 && (indelOp == 'D' || indelStart % 3 == 1))
            read.setCigar(buildTestCigarString(indelOp, 0, Math.max(16, leftmostAllowedAlignment), indelSize, readLength));
        else if ((indelOp == 'I' && indelSize == 4 && indelStart % 3 == 2) ||
                (indelOp == 'I' && indelSize == 2 && indelStart % 3 == 0) ||
                (indelOp == 'I' && indelSize == 1 && indelStart < 28 && indelStart % 3 == 2))
            read.setCigar(buildTestCigarString(indelOp, 0, Math.max(Math.max(indelStart - indelSize, 16), leftmostAllowedAlignment), indelSize, readLength));
        else
            return originalCigar;

        return read.getCigar();
    }

    @Test(dataProvider = "LeftAlignIndelWithLimitDataProvider")
    public void testLeftAlignIndelWithLimit(final Cigar originalCigar, final Cigar expectedCigar, final byte[] reference, final byte[] read, final int repeatLength, final int leftmostAllowedAlignment) {
        final Cigar actualCigar = AlignmentUtils.leftAlignIndel(originalCigar, reference, read, 0, 0, leftmostAllowedAlignment, true);
        Assert.assertTrue(expectedCigar.equals(actualCigar), "Wrong left alignment detected for cigar " + originalCigar.toString() + " to " + actualCigar.toString() + " but expected " + expectedCigar.toString() + " with repeat length " + repeatLength + " and leftmostAllowedAlignment " + leftmostAllowedAlignment);
    }

    @DataProvider(name = "LeftAlignIndelStartOfRead")
    public Object[][] makeLeftAlignIndelStartOfReadDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        final byte[] repeat1Reference = "XXXXXXXXXXXXABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat2Reference = "XYXYXYXYXYXYABCDEFGHIJKLMNOP".getBytes();
        final byte[] repeat3Reference = "XYZXYZXYZXYZABCDEFGHIJKLMNOP".getBytes();
        final int referenceLength = repeat1Reference.length;
        final int repeatLength=12;

        for (int indelStart = 1; indelStart < repeat1Reference.length; indelStart++) {
            for (final int indelSize : Arrays.asList(0, 1, 2, 3, 4)) {
                for (final char indelOp : Arrays.asList('D', 'I')) {
                    if (indelOp == 'D' && indelStart + indelSize >= repeatLength)
                        continue;
                    final int readLength = referenceLength - (indelOp == 'D' ? indelSize : -indelSize);

                    // create the original CIGAR string
                    final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, readLength);
                    read.setCigar(buildTestCigarString(indelSize == 0 ? 'M' : indelOp, 0, indelStart, indelSize, readLength));
                    final Cigar originalCigar = read.getCigar();

                    Cigar expectedCigar = originalCigar;

                    //create expected CIGAR string if the indel will be realigned
                    if (indelStart < 13 && indelSize > 0) {
                        if (indelOp == 'D') {
                            final List<CigarElement> elements = Stream.of(new CigarElement(readLength, CigarOperator.M)).collect(Collectors.toList());

                            expectedCigar = new Cigar(elements);
                        } else {
                            final List<CigarElement> elements = Stream.of(new CigarElement(indelSize, CigarOperator.I), new CigarElement(referenceLength, CigarOperator.M)).collect(Collectors.toList());
                            expectedCigar = new Cigar(elements);
                        }
                    }

                    //create readString
                    final byte[] readString1 = makeReadStringStartOfRead(repeat1Reference, indelOp, indelStart, indelSize, readLength, 1);
                    tests.add(new Object[]{originalCigar, expectedCigar, repeat1Reference, readString1});
                    if( indelSize % 2 == 0) {
                        final byte[] readString2 = makeReadStringStartOfRead(repeat2Reference, indelOp, indelStart, indelSize, readLength, 2);
                        tests.add(new Object[]{originalCigar, expectedCigar, repeat2Reference, readString2});
                    }
                    if (indelSize % 3 ==0) {
                        final byte[] readString3 = makeReadStringStartOfRead(repeat3Reference, indelOp, indelStart, indelSize, readLength, 3);
                        tests.add(new Object[]{originalCigar, expectedCigar, repeat3Reference, readString3});
                    }
                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    private static byte[] makeReadStringStartOfRead(final byte[] reference, final char indelOp, final int indelStart, final int indelSize, final int readLength, final int repeatLength) {
        final byte[] readString = new byte[readLength];

        if ( indelOp == 'D' && indelSize > 0 ) {
            System.arraycopy(reference, 0, readString, 0, indelStart);
            System.arraycopy(reference, indelStart + indelSize, readString, indelStart, readLength - indelStart);
        } else if ( indelOp == 'I' && indelSize > 0 ) {
            System.arraycopy(reference, 0, readString, 0, indelStart);
            for ( int i = indelStart; i < indelSize+indelStart; i++ ) {
                if ( i % repeatLength == 0 )
                    readString[i] = 'X';
                else if ( i % repeatLength == 1 )
                    readString[i] = 'Y';
                else
                    readString[i] = 'Z';
            }
            System.arraycopy(reference, indelStart, readString, indelStart + indelSize, readLength - indelStart - indelSize);
        } else {
            System.arraycopy(reference, 0, readString, 0, readLength);
        }

        return readString;
    }


    @Test(dataProvider = "LeftAlignIndelStartOfRead")
    public void testLeftAlignToStartOfRead(final Cigar originalCigar, final Cigar expectedCigar, final byte[] reference, final byte[] read) {
        final Cigar actualCigar = AlignmentUtils.leftAlignIndel(originalCigar, reference, read, 0, 0, true);
        Assert.assertTrue(expectedCigar.equals(actualCigar));
    }


    //////////////////////////////////////////
    // Test AlignmentUtils.IsIndelAlignedTooFarLeft //
    //////////////////////////////////////////

    @DataProvider(name = "IsIndelAlignedTooFarLeftDataProvider")
    public Object[][] makeIsIndelALignedTooFarLeftDataProvider() {
        List<Object[]> tests = new ArrayList<>();
        for (int nM = 0; nM < 3; nM++) {
            for (int nN = 0; nN < 3; nN++) {
                for (int nEq = 0; nEq < 3; nEq++) {
                    for (int nX = 0; nX < 3; nX++) {
                        final int totalBefore = nM + nN + nEq + nX;
                        if (totalBefore == 0) {
                            continue;
                        }
                        for (int leftmostLimit = totalBefore - 2; leftmostLimit < totalBefore + 2; leftmostLimit++) {
                            if (leftmostLimit < 0) {
                                continue;
                            }
                            boolean expected = leftmostLimit > totalBefore;
                            ArrayList<CigarElement> elements = new ArrayList<>();
                            if (nM > 0) {
                                elements.add(new CigarElement(nM, CigarOperator.M));
                            }
                            if (nN > 0) {
                                elements.add(new CigarElement(nN, CigarOperator.N));
                            }
                            if (nEq > 0) {
                                elements.add(new CigarElement(nEq, CigarOperator.EQ));
                            }
                            if (nX > 0) {
                                elements.add(new CigarElement(nX, CigarOperator.X));
                            }
                            tests.add(new Object[]{new Cigar(elements), leftmostLimit, false});
                            ArrayList<CigarElement> elementsInsertion = elements;
                            elementsInsertion.add(new CigarElement(10, CigarOperator.I));
                            tests.add(new Object[]{new Cigar(elementsInsertion), leftmostLimit, expected});
                            ArrayList<CigarElement> elementsDeletion = elements;
                            elementsDeletion.add(new CigarElement(10, CigarOperator.I));
                            tests.add(new Object[]{new Cigar(elementsDeletion), leftmostLimit, expected});

                            //test clippings
                            ArrayList<CigarElement> elementsSoftClippedInsertion = elementsInsertion;
                            elementsSoftClippedInsertion.add(0, new CigarElement(2, CigarOperator.S));
                            tests.add(new Object[]{new Cigar(elementsSoftClippedInsertion), leftmostLimit, expected});
                            ArrayList<CigarElement> elementsHardClippedInsertion = elementsInsertion;
                            elementsHardClippedInsertion.add(0, new CigarElement(2, CigarOperator.H));
                            tests.add(new Object[]{new Cigar(elementsHardClippedInsertion), leftmostLimit, expected});
                            ArrayList<CigarElement> elementsPaddedInsertion = elementsInsertion;
                            elementsPaddedInsertion.add(0, new CigarElement(2, CigarOperator.P));
                            tests.add(new Object[]{new Cigar(elementsPaddedInsertion), leftmostLimit, expected});

                            ArrayList<CigarElement> elementsSoftClippedDeletion = elementsDeletion;
                            elementsSoftClippedDeletion.add(0, new CigarElement(2, CigarOperator.S));
                            tests.add(new Object[]{new Cigar(elementsSoftClippedDeletion), leftmostLimit, expected});
                            ArrayList<CigarElement> elementsHardClippedDeletion = elementsDeletion;
                            elementsHardClippedDeletion.add(0, new CigarElement(2, CigarOperator.H));
                            tests.add(new Object[]{new Cigar(elementsHardClippedDeletion), leftmostLimit, expected});
                            ArrayList<CigarElement> elementsPaddedDeletion = elementsDeletion;
                            elementsPaddedDeletion.add(0, new CigarElement(2, CigarOperator.P));
                            tests.add(new Object[]{new Cigar(elementsPaddedDeletion), leftmostLimit, expected});

                        }
                    }

                }
            }
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "IsIndelAlignedTooFarLeftDataProvider")
    public void testIsIndelAlignedTooFarLeft(final Cigar cigar, final int leftmostAllowedAlignment, final boolean expected) {
        Assert.assertEquals(AlignmentUtils.isIndelAlignedTooFarLeft(cigar, leftmostAllowedAlignment), expected);
    }

    //////////////////////////////////////////
    // Test AlignmentUtils.trimCigarByReference() //
    //////////////////////////////////////////

    @DataProvider(name = "TrimCigarData")
    public Object[][] makeTrimCigarData() {
        List<Object[]> tests = new ArrayList<>();

        for ( final CigarOperator op : Arrays.asList(CigarOperator.D, CigarOperator.EQ, CigarOperator.X, CigarOperator.M) ) {
            for ( int myLength = 1; myLength < 6; myLength++ ) {
                for ( int start = 0; start < myLength - 1; start++ ) {
                    for ( int end = start; end < myLength; end++ ) {
                        final int length = end - start + 1;

                        final List<CigarOperator> padOps = Arrays.asList(CigarOperator.D, CigarOperator.M);
                        for ( final CigarOperator padOp: padOps) {
                            for ( int leftPad = 0; leftPad < 2; leftPad++ ) {
                                for ( int rightPad = 0; rightPad < 2; rightPad++ ) {
                                    tests.add(new Object[]{
                                            (leftPad > 0 ? leftPad + padOp.toString() : "") + myLength + op.toString() + (rightPad > 0 ? rightPad + padOp.toString() : ""),
                                            start + leftPad,
                                            end + leftPad,
                                            length + op.toString()});
                                }
                            }
                        }
                    }
                }
            }
        }

        for ( final int leftPad : Arrays.asList(0, 1, 2, 5) ) {
            for ( final int rightPad : Arrays.asList(0, 1, 2, 5) ) {
                final int length = leftPad + rightPad;
                if ( length > 0 ) {
                    for ( final int insSize : Arrays.asList(1, 10) ) {
                        for ( int start = 0; start <= leftPad; start++ ) {
                            for ( int stop = leftPad; stop < length; stop++ ) {
                                final int leftPadRemaining = leftPad - start;
                                final int rightPadRemaining = stop - leftPad + 1;
                                final String insC = insSize + "I";
                                tests.add(new Object[]{
                                        leftPad + "M" + insC + rightPad + "M",
                                        start,
                                        stop,
                                        (leftPadRemaining > 0 ? leftPadRemaining + "M" : "") + insC + (rightPadRemaining > 0 ? rightPadRemaining + "M" : "")
                                });
                            }
                        }
                    }
                }
            }
        }

        tests.add(new Object[]{"3M2D4M", 0, 8, "3M2D4M"});
        tests.add(new Object[]{"3M2D4M", 2, 8, "1M2D4M"});
        tests.add(new Object[]{"3M2D4M", 2, 6, "1M2D2M"});
        tests.add(new Object[]{"3M2D4M", 3, 6, "2D2M"});
        tests.add(new Object[]{"3M2D4M", 4, 6, "1D2M"});
        tests.add(new Object[]{"3M2D4M", 5, 6, "2M"});
        tests.add(new Object[]{"3M2D4M", 6, 6, "1M"});

        tests.add(new Object[]{"2M3I4M", 0, 5, "2M3I4M"});
        tests.add(new Object[]{"2M3I4M", 1, 5, "1M3I4M"});
        tests.add(new Object[]{"2M3I4M", 1, 4, "1M3I3M"});
        tests.add(new Object[]{"2M3I4M", 2, 4, "3I3M"});
        tests.add(new Object[]{"2M3I4M", 2, 3, "3I2M"});
        tests.add(new Object[]{"2M3I4M", 2, 2, "3I1M"});
        tests.add(new Object[]{"2M3I4M", 3, 4, "2M"});
        tests.add(new Object[]{"2M3I4M", 3, 3, "1M"});
        tests.add(new Object[]{"2M3I4M", 4, 4, "1M"});

        // this doesn't work -- but I'm not sure it should
        //        tests.add(new Object[]{"2M3I4M", 2, 1, "3I"});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TrimCigarData", enabled = ! DEBUG)
    public void testTrimCigar(final String cigarString, final int start, final int length, final String expectedCigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final Cigar expectedCigar = TextCigarCodec.decode(expectedCigarString);
        final Cigar actualCigar = AlignmentUtils.trimCigarByReference(cigar, start, length);
        Assert.assertEquals(actualCigar, expectedCigar);
    }

    @DataProvider(name = "TrimCigarByBasesData")
    public Object[][] makeTrimCigarByBasesData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{"2M3I4M", 0, 8, "2M3I4M"});
        tests.add(new Object[]{"2M3I4M", 1, 8, "1M3I4M"});
        tests.add(new Object[]{"2M3I4M", 2, 8, "3I4M"});
        tests.add(new Object[]{"2M3I4M", 3, 8, "2I4M"});
        tests.add(new Object[]{"2M3I4M", 4, 8, "1I4M"});
        tests.add(new Object[]{"2M3I4M", 4, 7, "1I3M"});
        tests.add(new Object[]{"2M3I4M", 4, 6, "1I2M"});
        tests.add(new Object[]{"2M3I4M", 4, 5, "1I1M"});
        tests.add(new Object[]{"2M3I4M", 4, 4, "1I"});
        tests.add(new Object[]{"2M3I4M", 5, 5, "1M"});

        tests.add(new Object[]{"2M2D2I", 0, 3, "2M2D2I"});
        tests.add(new Object[]{"2M2D2I", 1, 3, "1M2D2I"});
        tests.add(new Object[]{"2M2D2I", 2, 3, "2D2I"});
        tests.add(new Object[]{"2M2D2I", 3, 3, "1I"});
        tests.add(new Object[]{"2M2D2I", 2, 2, "2D1I"});
        tests.add(new Object[]{"2M2D2I", 1, 2, "1M2D1I"});
        tests.add(new Object[]{"2M2D2I", 0, 1, "2M2D"});
        tests.add(new Object[]{"2M2D2I", 1, 1, "1M2D"});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TrimCigarByBasesData", enabled = !DEBUG)
    public void testTrimCigarByBase(final String cigarString, final int start, final int length, final String expectedCigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final Cigar expectedCigar = TextCigarCodec.decode(expectedCigarString);
        final Cigar actualCigar = AlignmentUtils.trimCigarByBases(cigar, start, length);
        Assert.assertEquals(actualCigar, expectedCigar);
    }

    //////////////////////////////////////////
    // Test AlignmentUtils.applyCigarToCigar() //
    //////////////////////////////////////////

    @DataProvider(name = "ApplyCigarToCigarData")
    public Object[][] makeApplyCigarToCigarData() {
        List<Object[]> tests = new ArrayList<>();

        for ( int i = 1; i < 5; i++ )
            tests.add(new Object[]{i + "M", i + "M", i + "M"});

//        * ref   : ACGTAC
//        * hap   : AC---C  - 2M3D1M
//        * read  : AC---C  - 3M
//        * result: AG---C => 2M3D
        tests.add(new Object[]{"3M", "2M3D1M", "2M3D1M"});

//        * ref   : ACxG-TA
//        * hap   : AC-G-TA  - 2M1D3M
//        * read  : AC-GxTA  - 3M1I2M
//        * result: AC-GxTA => 2M1D1M1I2M
        tests.add(new Object[]{"3M1I2M", "2M1D3M", "2M1D1M1I2M"});

//        * ref   : ACGTA
//        * hap   : A-GTA  - 1M1D3M
//        * read  : A--TA  - 1M1D2M
//        * result: A--TA => 1M2D2M
        tests.add(new Object[]{"1M1D2M", "1M1D3M", "1M2D2M"});

//        * ref   : ACG-TA
//        * hap   : A-GxTA  - 1M1D1M1I2M
//        * read  : A---TA  - 1M2D2M
//        * result: A---TA => 1M2D2M
        tests.add(new Object[]{"1M2D2M", "1M1D1M1I2M", "1M2D2M"});

//        * ref   : A-CGTA
//        * hap   : A-CGTA  - 5M
//        * read  : AxCGTA  - 1M1I4M
//        * result: AxCGTA => 1M1I4M
        tests.add(new Object[]{"1M1I4M", "5M", "1M1I4M"});

//        * ref   : ACGTA
//        * hap   : ACGTA  - 5M
//        * read  : A--TA  - 1M2D2M
//        * result: A--TA => 1M2D2M
        tests.add(new Object[]{"1M2D2M", "5M", "1M2D2M"});

//        * ref   : AC-GTA
//        * hap   : ACxGTA  - 2M1I3M
//        * read  : A--GTA  - 1M2D3M
//        * result: A--GTA => 1M1D3M
        tests.add(new Object[]{"108M14D24M2M18I29M92M1000M", "2M1I3M", "2M1I3M"});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ApplyCigarToCigarData", enabled = !DEBUG)
    public void testApplyCigarToCigar(final String firstToSecondString, final String secondToThirdString, final String expectedCigarString) {
        final Cigar firstToSecond = TextCigarCodec.decode(firstToSecondString);
        final Cigar secondToThird = TextCigarCodec.decode(secondToThirdString);
        final Cigar expectedCigar = TextCigarCodec.decode(expectedCigarString);
        final Cigar actualCigar = AlignmentUtils.applyCigarToCigar(firstToSecond, secondToThird);
        Assert.assertEquals(actualCigar, expectedCigar);
    }

    //////////////////////////////////////////
    // Test AlignmentUtils.applyCigarToCigar() //
    //////////////////////////////////////////

    @DataProvider(name = "ReadOffsetFromCigarData")
    public Object[][] makeReadOffsetFromCigarData() {
        List<Object[]> tests = new ArrayList<>();

        final int SIZE = 10;
        for ( int i = 0; i < SIZE; i++ ) {
            tests.add(new Object[]{SIZE + "M", i, i});
        }

        //          0123ii45
        // ref    : ACGT--AC
        // hap    : AC--xxAC (2M2D2I2M)
        // ref.pos: 01    45
        tests.add(new Object[]{"2M2D2I2M", 0, 0});
        tests.add(new Object[]{"2M2D2I2M", 1, 1});
        tests.add(new Object[]{"2M2D2I2M", 2, 4});
        tests.add(new Object[]{"2M2D2I2M", 3, 4});
        tests.add(new Object[]{"2M2D2I2M", 4, 4});
        tests.add(new Object[]{"2M2D2I2M", 5, 5});

        // 10132723 - 10132075 - 500 = 148
        // what's the offset of the first match after the I?
        // 108M + 14D + 24M + 2M = 148
        // What's the offset of the first base that is after the I?
        // 108M + 24M + 2M + 18I = 134M + 18I = 152 - 1 = 151
        tests.add(new Object[]{"108M14D24M2M18I29M92M", 0, 0});
        tests.add(new Object[]{"108M14D24M2M18I29M92M", 107, 107});
        tests.add(new Object[]{"108M14D24M2M18I29M92M", 108, 108 + 14}); // first base after the deletion

        tests.add(new Object[]{"108M14D24M2M18I29M92M", 132, 132+14}); // 2 before insertion
        tests.add(new Object[]{"108M14D24M2M18I29M92M", 133, 133+14}); // last base before insertion

        // entering into the insertion
        for ( int i = 0; i < 18; i++ ) {
            tests.add(new Object[]{"108M14D24M2M18I29M92M", 134+i, 148}); // inside insertion
        }
        tests.add(new Object[]{"108M14D24M2M18I29M92M", 134+18, 148}); // first base after insertion matches at same as insertion
        tests.add(new Object[]{"108M14D24M2M18I29M92M", 134+18+1, 149});
        tests.add(new Object[]{"108M14D24M2M18I29M92M", 134+18+2, 150});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadOffsetFromCigarData", enabled = !DEBUG)
    public void testReadOffsetFromCigar(final String cigarString, final int startOnCigar, final int expectedOffset) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final int actualOffset = AlignmentUtils.calcFirstBaseMatchingReferenceInCigar(cigar, startOnCigar);
        Assert.assertEquals(actualOffset, expectedOffset);
    }

    //////////////////////////////////////////
    // Test AlignmentUtils.addCigarElements() //
    //////////////////////////////////////////

    @DataProvider(name = "AddCigarElementsData")
    public Object[][] makeAddCigarElementsData() {
        List<Object[]> tests = new ArrayList<>();

        final int SIZE = 10;
        for ( final CigarOperator op : Arrays.asList(CigarOperator.I, CigarOperator.M, CigarOperator.S, CigarOperator.EQ, CigarOperator.X)) {
            for ( int start = 0; start < SIZE; start++ ) {
                for ( int end = start; end < SIZE * 2; end ++ ) {
                    for ( int pos = 0; pos < SIZE * 3; pos++ ) {
                        int length = 0;
                        for ( int i = 0; i < SIZE; i++ ) length += (i+pos) >= start && (i+pos) <= end ? 1 : 0;
                        tests.add(new Object[]{SIZE + op.toString(), pos, start, end, length > 0 ? length + op.toString() : "*"});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "AddCigarElementsData", enabled = !DEBUG)
    public void testAddCigarElements(final String cigarString, final int pos, final int start, final int end, final String expectedCigarString) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final CigarElement elt = cigar.getCigarElement(0);
        final Cigar expectedCigar = TextCigarCodec.decode(expectedCigarString);

        final List<CigarElement> elts = new LinkedList<>();
        final int actualEndPos = AlignmentUtils.addCigarElements(elts, pos, start, end, elt);

        Assert.assertEquals(actualEndPos, pos + elt.getLength());
        Assert.assertEquals(AlignmentUtils.consolidateCigar(new Cigar(elts)), expectedCigar);
    }

    @DataProvider(name = "GetBasesCoveringRefIntervalData")
    public Object[][] makeGetBasesCoveringRefIntervalData() {
        List<Object[]> tests = new ArrayList<>();

        // matches
        // 0123
        // ACGT
        tests.add(new Object[]{"ACGT", 0, 3, "4M", "ACGT"});
        tests.add(new Object[]{"ACGT", 1, 3, "4M", "CGT"});
        tests.add(new Object[]{"ACGT", 1, 2, "4M", "CG"});
        tests.add(new Object[]{"ACGT", 1, 1, "4M", "C"});

        // deletions
        // 012345
        // AC--GT
        tests.add(new Object[]{"ACGT", 0, 5, "2M2D2M", "ACGT"});
        tests.add(new Object[]{"ACGT", 1, 5, "2M2D2M", "CGT"});
        tests.add(new Object[]{"ACGT", 2, 5, "2M2D2M", null});
        tests.add(new Object[]{"ACGT", 3, 5, "2M2D2M", null});
        tests.add(new Object[]{"ACGT", 4, 5, "2M2D2M", "GT"});
        tests.add(new Object[]{"ACGT", 5, 5, "2M2D2M", "T"});
        tests.add(new Object[]{"ACGT", 0, 4, "2M2D2M", "ACG"});
        tests.add(new Object[]{"ACGT", 0, 3, "2M2D2M", null});
        tests.add(new Object[]{"ACGT", 0, 2, "2M2D2M", null});
        tests.add(new Object[]{"ACGT", 0, 1, "2M2D2M", "AC"});
        tests.add(new Object[]{"ACGT", 0, 0, "2M2D2M", "A"});

        // insertions
        // 01--23
        // ACTTGT
        tests.add(new Object[]{"ACTTGT", 0, 3, "2M2I2M", "ACTTGT"});
        tests.add(new Object[]{"ACTTGT", 1, 3, "2M2I2M", "CTTGT"});
        tests.add(new Object[]{"ACTTGT", 2, 3, "2M2I2M", "GT"});
        tests.add(new Object[]{"ACTTGT", 3, 3, "2M2I2M", "T"});
        tests.add(new Object[]{"ACTTGT", 0, 2, "2M2I2M", "ACTTG"});
        tests.add(new Object[]{"ACTTGT", 0, 1, "2M2I2M", "AC"});
        tests.add(new Object[]{"ACTTGT", 1, 2, "2M2I2M", "CTTG"});
        tests.add(new Object[]{"ACTTGT", 2, 2, "2M2I2M", "G"});
        tests.add(new Object[]{"ACTTGT", 1, 1, "2M2I2M", "C"});

        tests.add(new Object[]{"ACGT", 0, 1, "2M2I", "AC"});
        tests.add(new Object[]{"ACGT", 1, 1, "2M2I", "C"});
        tests.add(new Object[]{"ACGT", 0, 0, "2M2I", "A"});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "GetBasesCoveringRefIntervalData")
    public void testGetBasesCoveringRefInterval(final String basesString, final int refStart, final int refEnd, final String cigarString, final String expected) {
        final byte[] actualBytes = AlignmentUtils.getBasesCoveringRefInterval(refStart, refEnd, basesString.getBytes(), 0, TextCigarCodec.decode(cigarString));
        if ( expected == null )
            Assert.assertNull(actualBytes);
        else
            Assert.assertEquals(new String(actualBytes), expected);
    }

    @DataProvider(name = "StartsOrEndsWithInsertionOrDeletionData")
    public Object[][] makeStartsOrEndsWithInsertionOrDeletionData() {
        List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{"2M", false});
        tests.add(new Object[]{"1D2M", true});
        tests.add(new Object[]{"2M1D", true});
        tests.add(new Object[]{"2M1I", true});
        tests.add(new Object[]{"1I2M", true});
        tests.add(new Object[]{"1M1I2M", false});
        tests.add(new Object[]{"1M1D2M", false});
        tests.add(new Object[]{"1M1I2M1I", true});
        tests.add(new Object[]{"1M1I2M1D", true});
        tests.add(new Object[]{"1D1M1I2M", true});
        tests.add(new Object[]{"1I1M1I2M", true});
        tests.add(new Object[]{"1M1I2M1I1M", false});
        tests.add(new Object[]{"1M1I2M1D1M", false});
        tests.add(new Object[]{"1M1D2M1D1M", false});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "StartsOrEndsWithInsertionOrDeletionData")
    public void testStartsOrEndsWithInsertionOrDeletion(final String cigar, final boolean expected) {
        Assert.assertEquals(AlignmentUtils.startsOrEndsWithInsertionOrDeletion(TextCigarCodec.decode(cigar)), expected);
    }

    @Test(dataProvider = "StartsOrEndsWithInsertionOrDeletionData")
    public void testRemoveTrailingDeletions(final String cigar, final boolean expected) {

        final Cigar originalCigar = TextCigarCodec.decode(cigar);
        final Cigar newCigar = AlignmentUtils.removeTrailingDeletions(originalCigar);

        Assert.assertEquals(originalCigar.equals(newCigar), !cigar.endsWith("D"));
    }
}
