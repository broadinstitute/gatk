package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Abstract test class to test implementations of {@link SmithWatermanAligner}
 *
 * Implementations of a {@link SmithWatermanAligner} can be tested by subclassing this test and overriding
 * {@link #getAligner()} to return an instance of the class to be tested
 */
@Test
public abstract class SmithWatermanAlignerAbstractUnitTest extends GATKBaseTest {

    /**
     * @return an aligner to be tested
     */
    protected abstract SmithWatermanAligner getAligner();

    protected static void printAlignment(final byte[] ref, final byte[] read, final SmithWatermanAlignment alignment, final SWOverhangStrategy overhangStrategy) {
        final StringBuilder bread = new StringBuilder();
        final StringBuilder bref = new StringBuilder();
        final StringBuilder match = new StringBuilder();

        int i = 0;
        int j = 0;

        final int offset = alignment.getAlignmentOffset();

        Cigar cigar = alignment.getCigar();

        if ( overhangStrategy != SWOverhangStrategy.SOFTCLIP ) {

            // we need to go through all the hassle below only if we do not do soft-clipping;
            // otherwise offset is never negative
            if ( offset < 0 ) {
                for (  ; j < (-offset) ; j++ ) {
                    bread.append((char) read[j]);
                    bref.append(' ');
                    match.append(' ');
                }
                // at negative offsets, our cigar's first element carries overhanging bases
                // that we have just printed above. Tweak the first element to
                // exclude those bases. Here we create a new list of cigar elements, so the original
                // list/original cigar are unchanged (they are unmodifiable anyway!)

                final List<CigarElement> tweaked = new ArrayList<>();
                tweaked.addAll(cigar.getCigarElements());
                tweaked.set(0,new CigarElement(cigar.getCigarElement(0).getLength()+offset,
                        cigar.getCigarElement(0).getOperator()));
                cigar = new Cigar(tweaked);
            }
        }

        if ( offset > 0 ) { // note: the way this implementation works, cigar will ever start from S *only* if read starts before the ref, i.e. offset = 0
            for (; i < alignment.getAlignmentOffset() ; i++ ) {
                bref.append((char) ref[i]);
                bread.append(' ');
                match.append(' ');
            }
        }

        for ( final CigarElement e : cigar.getCigarElements() ) {
            switch (e.getOperator()) {
                case M :
                    for ( int z = 0 ; z < e.getLength() ; z++, i++, j++  ) {
                        bref.append((i< ref.length)?(char) ref[i]:' ');
                        bread.append((j < read.length)?(char) read[j]:' ');
                        match.append( ( i< ref.length && j < read.length ) ? (ref[i] == read[j] ? '.':'*' ) : ' ' );
                    }
                    break;
                case I :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append('-');
                        bread.append((char) read[j]);
                        match.append('I');
                    }
                    break;
                case S :
                    for ( int z = 0 ; z < e.getLength(); z++, j++ ) {
                        bref.append(' ');
                        bread.append((char) read[j]);
                        match.append('S');
                    }
                    break;
                case D:
                    for ( int z = 0 ; z < e.getLength(); z++ , i++ ) {
                        bref.append((char) ref[i]);
                        bread.append('-');
                        match.append('D');
                    }
                    break;
                default:
                    throw new GATKException("Unexpected Cigar element:" + e.getOperator());
            }
        }
        for (; i < ref.length; i++ ) bref.append((char) ref[i]);
        for (; j < read.length; j++ ) bread.append((char) read[j]);

        int pos = 0 ;
        final int maxlength = Math.max(match.length(), Math.max(bread.length(), bref.length()));
        while ( pos < maxlength ) {
            printCautiously(match, pos, 100);
            printCautiously(bread, pos, 100);
            printCautiously(bref, pos, 100);
            System.out.println();
            pos += 100;
        }
    }

    /** String builder's substring is extremely stupid: instead of trimming and/or returning an empty
     * string when one end/both ends of the interval are out of range, it crashes with an
     * exception. This utility function simply prints the substring if the interval is within the index range
     * or trims accordingly if it is not.
     */
    private static void printCautiously(final StringBuilder s, final int start, final int width) {
        if ( start >= s.length() ) {
            System.out.println();
            return;
        }
        final int end = Math.min(start + width, s.length());
        System.out.println(s.substring(start,end));
    }

    @DataProvider(name = "ComplexReadAlignedToRef")
    public Object[][] makeComplexReadAlignedToRef() {
        return new Object[][] {
                {"AAAGGACTGACTG", "ACTGACTGACTG", 1, "12M"}
        };
    }

    @Test(dataProvider = "ComplexReadAlignedToRef")
    public void testReadAlignedToRefComplexAlignment(final String reference, final String read, final int expectedStart, final String expectedCigar) {
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.SOFTCLIP);
    }

    @DataProvider(name = "OddNoAlignment")
    public Object[][] makeOddNoAlignment() {
        final String ref1     = "AAAGACTACTG";
        final String read1    = "AACGGACACTG";
        return new Object[][] {
                {ref1, read1, new SWParameters(50, -100, -220, -12), 1,  "2M2I3M1D4M"},
                {ref1, read1, new SWParameters(200, -50, -300, -22), 0, "11M"}
        };
    }

    @Test(dataProvider = "OddNoAlignment")
    public void testOddNoAlignment(final String reference, final String read, final SWParameters weights,
                                   final int expectedStart, final String expectedCigar) {
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, weights,
                                       SWOverhangStrategy.SOFTCLIP);
    }

    @Test
    public void testIndelsAtStartAndEnd() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match + "GGG";
        final int expectedStart = 3;
        final String expectedCigar = "5M3S";
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.SOFTCLIP);
    }

    @Test
    public void testDegenerateAlignmentWithIndelsAtBothEnds() {
        logger.warn("testDegenerateAlignmentWithIndelsAtBothEnds");
        final String ref = "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final String alt =               "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final int expectedStart = 14;
        final String expectedCigar = "31M20S";
        assertAlignmentMatchesExpected(ref, alt, expectedStart, expectedCigar, SmithWatermanAligner.STANDARD_NGS, SWOverhangStrategy.SOFTCLIP);
    }

    @Test
    public void  testForIdenticalAlignmentsWithDifferingFlankLengths() {
        //This test is designed to ensure that the indels are correctly placed
        //if the region flanking these indels is extended by a varying amount.
        //It checks for problems caused by floating point rounding leading to different
        //paths being selected.

        //Create two versions of the same sequence with different flanking regions.
        final byte[] paddedRef = "GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT"
                .getBytes();
        final byte[] paddedHap = "GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA--GGGCC---------------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT"
                .replace("-", "")
                .getBytes();
        final byte[] notPaddedRef = "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA"
                .getBytes();
        final byte[] notPaddedHap = "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA---------GGGCC--------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA"
                .replace("-", "")
                .getBytes();
        //a simplified version of the getCigar routine in the haplotype caller to align these
        final String SW_PAD = "NNNNNNNNNN";
        final String paddedsRef = SW_PAD + new String(paddedRef) + SW_PAD;
        final String paddedsHap = SW_PAD + new String(paddedHap) + SW_PAD;
        final String notPaddedsRef = SW_PAD + new String(notPaddedRef) + SW_PAD;
        final String notpaddedsHap = SW_PAD + new String(notPaddedHap) + SW_PAD;
        final SmithWatermanAlignment paddedAlignment;
        final SmithWatermanAlignment notPaddedAlignment;
        try (final SmithWatermanAligner aligner = getAligner()) {
            paddedAlignment = aligner.align(paddedsRef.getBytes(), paddedsHap.getBytes(), CigarUtils.NEW_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
            notPaddedAlignment = aligner.align(notPaddedsRef.getBytes(), notpaddedsHap.getBytes(), CigarUtils.NEW_SW_PARAMETERS, SWOverhangStrategy.SOFTCLIP);
        }
        //Now verify that the two sequences have the same alignment and not match positions.
        final Cigar rawPadded = paddedAlignment.getCigar();
        final Cigar notPadded = notPaddedAlignment.getCigar();
        final List<CigarElement> paddedC = rawPadded.getCigarElements();
        final List<CigarElement> notPaddedC = notPadded.getCigarElements();
        Assert.assertEquals(paddedC.size(), notPaddedC.size());
        for (int i = 0; i < notPaddedC.size(); i++) {
            final CigarElement pc = paddedC.get(i);
            final CigarElement npc = notPaddedC.get(i);
            if (pc.getOperator() == CigarOperator.M && npc.getOperator() == CigarOperator.M) {
                continue;
            }
            final int l1 = pc.getLength();
            final int l2 = npc.getLength();
            Assert.assertEquals(l1, l2);
            Assert.assertEquals(pc.getOperator(), npc.getOperator());
        }
    }

    @DataProvider
    public Object[][] getSubstringMatchTests(){
        return new Object[][]{
                {3, "5M", SWOverhangStrategy.SOFTCLIP},
                {0, "3D5M", SWOverhangStrategy.INDEL},
                {0, "3D5M", SWOverhangStrategy.LEADING_INDEL},
                {3, "5M", SWOverhangStrategy.IGNORE}
        };
    }

    @Test(dataProvider = "getSubstringMatchTests")
    public void testSubstringMatch(int expectedStart, String expectedCigar, SWOverhangStrategy strategy) {
        final String matchingSection = "CCCCC";
        final String reference = "AAA" + matchingSection;
        final String read = matchingSection;
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SmithWatermanAligner.ORIGINAL_DEFAULT,
                                       strategy);
    }

    protected void assertAlignmentMatchesExpected(String reference, String read, int expectedStart, String expectedCigar, SWParameters weights, SWOverhangStrategy strategy) {
        try(final SmithWatermanAligner sw = getAligner()) {
            final SmithWatermanAlignment alignment = sw.align(reference.getBytes(), read.getBytes(), weights, strategy);
            printAlignment(reference.getBytes(), read.getBytes(), alignment, strategy);
            Assert.assertEquals(alignment.getAlignmentOffset(), expectedStart);
            Assert.assertEquals(alignment.getCigar().toString(), expectedCigar);
        }
    }

    @DataProvider
    public Object[][] getSubstringMatchLong(){
        return new Object[][]{
                {359, "7M", SWOverhangStrategy.SOFTCLIP},
                {0, "1M358D6M29D", SWOverhangStrategy.INDEL},
                {0, "1M1D6M", SWOverhangStrategy.LEADING_INDEL},
                {359, "7M", SWOverhangStrategy.IGNORE}
        };
    }

    @Test(dataProvider = "getSubstringMatchLong")
    public void testSubstringMatchLong(int expectedStart, String expectedCigar, SWOverhangStrategy strategy) {
        final String reference = "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC";
        final String read      = "AAAAAAA";
        assertAlignmentMatchesExpected(reference, read, expectedStart, expectedCigar, SmithWatermanAligner.ORIGINAL_DEFAULT, strategy);
    }
}
