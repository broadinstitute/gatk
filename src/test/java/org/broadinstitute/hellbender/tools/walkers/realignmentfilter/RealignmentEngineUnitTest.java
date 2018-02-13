package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class RealignmentEngineUnitTest {
    private static final SAMFileHeader HEADER = ArtificialReadUtils.createArtificialSamHeader();
    private static final int REF_INDEX = 0;
    private static final String SOURCE = "source";
    private static final String CONTIG = "chr1";
    private static final int TOL = 1;

    @Test
    public void testSupportsVariant() {
        final GATKRead refRead = makeRead("refRead", 1, "ACGTACGT", "8M");

        // G -> C at position 3
        final GATKRead subst1 = makeRead("subst1", 1, "ACCTACGT", "8M");
        final VariantContext vcSubst1 = makeVariant(3, "G", "C");

        // A -> T at position 5
        final GATKRead subst2 = makeRead("subst2", 5, "TCGTGGGG", "8M");
        final VariantContext vcSubst2 = makeVariant(5, "A", "T");

        // G -> GA at position 3
        final GATKRead ins1 = makeRead("ins1", 1, "ACGATACG", "3M1I4M");
        final VariantContext vcIns1 = makeVariant(3, "G", "GA");

        // AC -> A at position 5
        final GATKRead del1 = makeRead("del1", 5, "AGTGGGGG", "1M1D7M");
        final VariantContext vcDel1 = makeVariant(5, "AC", "A");

        // read[n] supports vcs[n]
        final List<GATKRead> reads = Arrays.asList(subst1, subst2, ins1, del1, refRead);
        final List<VariantContext> vcs = Arrays.asList(vcSubst1, vcSubst2, vcIns1, vcDel1);

        for (int m = 0; m < reads.size(); m++) {
            for (int n = 0; n < vcs.size(); n++) {
                Assert.assertEquals(RealignmentEngine.supportsVariant(reads.get(m), vcs.get(n), TOL), m == n);
            }
        }
    }

    @Test
    public void testCheckRealignments() {
        final BwaMemAlignment aln90 = makeAlignment(0,90, 1, false);
        final BwaMemAlignment aln80 = makeAlignment(1,80, 1, true);
        final BwaMemAlignment aln70 = makeAlignment(0,70, 1, false);
        Assert.assertTrue(RealignmentEngine.checkAlignments(Arrays.asList(aln70), 1).isGood());
        Assert.assertTrue(RealignmentEngine.checkAlignments(Arrays.asList(aln90, aln80), 5).isGood());
        Assert.assertTrue(RealignmentEngine.checkAlignments(Arrays.asList(aln90, aln80, aln70), 5).isGood());
        Assert.assertFalse(RealignmentEngine.checkAlignments(Arrays.asList(aln90, aln80, aln70), 11).isGood());
    }

    @Test
    public void testFindPlausiblePairs() {
        final int maxFragmentLength = 1000;
        final BwaMemAlignment aln1 = makeAlignment(0, 90, 100, true);
        final BwaMemAlignment aln2 = makeAlignment(0, 90, 100, false);
        final BwaMemAlignment aln3 = makeAlignment(0, 90, 300, true);
        final BwaMemAlignment aln4 = makeAlignment(0, 90, 300, false);
        final BwaMemAlignment aln5 = makeAlignment(1, 90, 100, true);
        final BwaMemAlignment aln6 = makeAlignment(1, 90, 100, false);
        final BwaMemAlignment aln7 = makeAlignment(1, 90, 300, true);
        final BwaMemAlignment aln8 = makeAlignment(0, 90, 2000, false);
        final BwaMemAlignment aln9 = makeAlignment(0, 90, 2100, true);
        final BwaMemAlignment aln10 = makeAlignment(0, 90, 2200, true);

        // reads are 1,3,5,7,9; mates are 2,4,6,8,10 -- only pair is 9 & 10
        final List<BwaMemAlignment> readAlignments1 = Arrays.asList(aln1, aln3, aln5, aln7, aln9);
        final List<BwaMemAlignment> mateAlignments1 = Arrays.asList(aln2, aln4, aln6, aln8, aln10);
        final List<Pair<BwaMemAlignment, BwaMemAlignment>> pairs1 = RealignmentEngine.findPlausiblePairs(readAlignments1, mateAlignments1, maxFragmentLength);
        Assert.assertEquals(pairs1.size(),1);
        Assert.assertEquals(pairs1.get(0).getFirst(), aln9);
        Assert.assertEquals(pairs1.get(0).getSecond(), aln10);

        // reads are 1,2,3,4,5; mates are 6,7,8,9,10 -- only pair is 5 & 7
        final List<BwaMemAlignment> readAlignments2 = Arrays.asList(aln1, aln2, aln3, aln4, aln5);
        final List<BwaMemAlignment> mateAlignments2 = Arrays.asList(aln6, aln7, aln8, aln9, aln10);
        final List<Pair<BwaMemAlignment, BwaMemAlignment>> pairs2 = RealignmentEngine.findPlausiblePairs(readAlignments2, mateAlignments2, maxFragmentLength);
        Assert.assertEquals(pairs2.size(),1);
        Assert.assertEquals(pairs2.get(0).getFirst(), aln5);
        Assert.assertEquals(pairs2.get(0).getSecond(), aln7);

        // reads are 1,2,5,6,9; mates are 3,4,7,8,10 -- pairs are 1 & 3, 2 & 4, 5 & 7, 9 & 10
        final List<BwaMemAlignment> readAlignments3 = Arrays.asList(aln1, aln2, aln5, aln6, aln9);
        final List<BwaMemAlignment> mateAlignments3 = Arrays.asList(aln3, aln4, aln7, aln8, aln10);
        final List<Pair<BwaMemAlignment, BwaMemAlignment>> pairs3 = RealignmentEngine.findPlausiblePairs(readAlignments3, mateAlignments3, maxFragmentLength);
        Assert.assertEquals(pairs3.size(),4);
    }

    private GATKRead makeRead(final String name, final int alignmentStart, final String baseString, final String cigar ) {
        final byte[] bases = baseString.getBytes();
        final byte[] qual = new byte[bases.length];
        Arrays.fill(qual, (byte) 30);
        return ArtificialReadUtils.createArtificialRead(HEADER, name, REF_INDEX, alignmentStart, bases, qual, cigar);
    }

    private VariantContext makeVariant(final int position, final String ref, final String alt) {
        final int end = position + ref.length() - 1;
        final Collection<Allele> alleles = Arrays.asList(Allele.create(ref, true), Allele.create(alt, false));
        return new VariantContextBuilder(SOURCE, CONTIG, position, end, alleles).make();
    }

    // checkAlignments only uses the alignerScore, so every else is a placeholder here
    private BwaMemAlignment makeAlignment(final int refId, final int score, final int refStart, final boolean isReverseStrand) {
        return new BwaMemAlignment(isReverseStrand ? 16 : 4, refId, refStart, refStart + 100,
        1, 100, 50, 1, score, score - 30,
        "100M", null, null, 0, refStart + 200, refStart + 300 );
    }
}