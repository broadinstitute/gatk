package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public final class SWPairwiseAlignmentUnitTest extends BaseTest {
    @DataProvider(name = "ComplexReadAlignedToRef")
    public Object[][] makeComplexReadAlignedToRef() {
        List<Object[]> tests = new ArrayList<>();

        final String ref1     = "ACTGACTGACTG";
        tests.add(new Object[]{"AAAGGACTGACTG", ref1, 1, "12M"});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ComplexReadAlignedToRef", enabled = true)
    public void testReadAlignedToRefComplexAlignment(final String reference, final String read, final int expectedStart, final String expectedCigar) {
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @DataProvider(name = "OddNoAlignment")
    public Object[][] makeOddNoAlignment() {
        List<Object[]> tests = new ArrayList<>();

        final String ref1     = "AAAGACTACTG";
        final String read1    = "AACGGACACTG";
        tests.add(new Object[]{ref1, read1, 50, -100, -220, -12, 1,  "2M2I3M1D4M"});
        tests.add(new Object[]{ref1, read1, 200, -50, -300, -22, 0, "11M"});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "OddNoAlignment", enabled = true)
    public void testOddNoAlignment(final String reference, final String read, final int match, final int mismatch, final int gap, final int gap_extend,
                                   final int expectedStart, final String expectedCigar) {
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(), new SWPairwiseAlignment.Parameters(match, mismatch, gap, gap_extend));
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testIndelsAtStartAndEnd() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match + "GGG";
        final int expectedStart = 3;
        final String expectedCigar = "5M3S";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes());
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testDegenerateAlignmentWithIndelsAtBothEnds() {
        logger.warn("testDegenerateAlignmentWithIndelsAtBothEnds");
        final String ref = "TGTGTGTGTGTGTGACAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final String alt =               "ACAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA";
        final int expectedStart = 14;
        final String expectedCigar = "31M20S";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(ref.getBytes(), alt.getBytes(), SWPairwiseAlignment.STANDARD_NGS);
        sw.printAlignment(ref.getBytes(), alt.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void  testForIdenticalAlignmentsWithDifferingFlankLengths() {
        //This test is designed to ensure that the indels are correctly placed
        //if the region flanking these indels is extended by a varying amount.
        //It checks for problems caused by floating point rounding leading to different
        //paths being selected.

        //Create two versions of the same sequence with different flanking regions.
        byte[] paddedRef="GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT".getBytes();
        byte[] paddedHap="GCGTCGCAGTCTTAAGGCCCCGCCTTTTCAGACAGCTTCCGCTGGGCCTGGGCCGCTGCGGGGCGGTCACGGCCCCTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA--GGGCC---------------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGACCGGGCCGAGCCGGGGGAAGGGCTCCGGTGACT".replace("-","").getBytes();
        byte[] notPaddedRef=                                                                           "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGAGGGGGCCCGGGGCCGCGTCCCTGGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA".getBytes();
        byte[] notPaddedHap=                                                                           "CTTTAAGCCTGAGCCCCGCCCCCTGGCTCCCCGCCCCCTCTTCTCCCCTCCCCCAAGCCAGCACCTGGTGCCCCGGCGGGTCGTGCGGCGCGGCGCTCCGCGGTGAGCGCCTGACCCCGA---------GGGCC--------GGGCCCTCCCCACCCTTGCGGTGGCCTCGCGGGTCCCAGGGGCGGGGCTGGAGCGGCAGCAGGGCCGGGGAGATGGGCGGTGGGGAGCGCGGGAGGGA".replace("-","").getBytes();
        //a simplified version of the getCigar routine in the haplotype caller to align these
        final String SW_PAD = "NNNNNNNNNN";
        final String paddedsRef = SW_PAD + new String(paddedRef) + SW_PAD;
        final String paddedsHap = SW_PAD + new String(paddedHap) + SW_PAD;
        final String notPaddedsRef = SW_PAD + new String(notPaddedRef) + SW_PAD;
        final String notpaddedsHap = SW_PAD + new String(notPaddedHap) + SW_PAD;
        final SWPairwiseAlignment paddedAlignment = new SWPairwiseAlignment( paddedsRef.getBytes(), paddedsHap.getBytes(), CigarUtils.NEW_SW_PARAMETERS );
        final SWPairwiseAlignment notPaddedAlignment = new SWPairwiseAlignment( notPaddedsRef.getBytes(), notpaddedsHap.getBytes(), CigarUtils.NEW_SW_PARAMETERS );
        //Now verify that the two sequences have the same alignment and not match positions.
        Cigar rawPadded = paddedAlignment.getCigar();
        Cigar notPadded= notPaddedAlignment.getCigar();
        List<CigarElement> paddedC=rawPadded.getCigarElements();
        List<CigarElement> notPaddedC=notPadded.getCigarElements();
        Assert.assertEquals(paddedC.size(), notPaddedC.size());
        for(int i=0;i<notPaddedC.size();i++)
        {
            CigarElement pc=paddedC.get(i);
            CigarElement npc=notPaddedC.get(i);
            if(pc.getOperator()== CigarOperator.M && npc.getOperator()== CigarOperator.M)
            {
                continue;
            }
            int l1=pc.getLength();
            int l2=npc.getLength();
            Assert.assertEquals(l1, l2);
            Assert.assertEquals(pc.getOperator(), npc.getOperator());
        }
    }

    @Test(enabled = true)
    public void testSubstringMatchSoftclip() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match;
        final int expectedStart = 3;
        final String expectedCigar = "5M";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.SOFTCLIP);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testSubstringMatchIndel() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match;
        final int expectedStart = 0;
        final String expectedCigar = "3D5M";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.INDEL);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testSubstringMatchLeadingIndel() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match;
        final int expectedStart = 0;
        final String expectedCigar = "3D5M";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.LEADING_INDEL);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testSubstringMatchIgnore() {
        final String match     = "CCCCC";
        final String reference = "AAA" + match;
        final String read      = match;
        final int expectedStart = 3;
        final String expectedCigar = "5M";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.IGNORE);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testSubstringMatchSoftclipLong() {
        final String reference = "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC";
        final String read      = "AAAAAAA";
        final int expectedStart = 359;
        final String expectedCigar = "7M";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.SOFTCLIP);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testSubstringMatchIndelLong() {
        final String reference = "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC";
        final String read      = "AAAAAAA";
        final int expectedStart = 0;
        final String expectedCigar = "1M358D6M29D";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.INDEL);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testSubstringMatchLeadingIndelLong() {
        final String reference = "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC";
        final String read      = "AAAAAAA";
        final int expectedStart = 0;
        final String expectedCigar = "1M1D6M";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.LEADING_INDEL);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }

    @Test(enabled = true)
    public void testSubstringMatchLeadingIgnoreLong() {
        final String reference = "ATAGAAAATAGTTTTTGGAAATATGGGTGAAGAGACATCTCCTCTTATGGAAAAAGGGATTCTAGAATTTAACAATAAATATTCCCAACTTTCCCCAAGGCTTTAAAATCTACCTTGAAGGAGCAGCTGATGTATTTCTAGAACAGACTTAGGTGTCTTGGTGTGGCCTGTAAAGAGATACTGTCTTTCTCTTTTGAGTGTAAGAGAGAAAGGACAGTCTACTCAATAAAGAGTGCTGGGAAAACTGAATATCCACACACAGAATAATAAAACTAGATCCTATCTCTCACCATATACAAAGATCAACTCAAAACAAATTAAAGACCTAAATGTAAGACAAGAAATTATAAAACTACTAGAAAAAAACACAAGGGAAATGCTTCAGGACATTGGC";
        final String read      = "AAAAAAA";
        final int expectedStart = 359;
        final String expectedCigar = "7M";
        final SWPairwiseAlignment sw = new SWPairwiseAlignment(reference.getBytes(), read.getBytes(),
                                                               SWPairwiseAlignment.ORIGINAL_DEFAULT,
                                                               SWPairwiseAlignment.OverhangStrategy.IGNORE);
        sw.printAlignment(reference.getBytes(), read.getBytes());
        Assert.assertEquals(sw.getAlignmentStart2wrt1(), expectedStart);
        Assert.assertEquals(sw.getCigar().toString(), expectedCigar);
    }
}
