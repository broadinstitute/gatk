package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class FragmentDepthPerAlleleBySampleUnitTest extends GATKBaseTest {

    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final String SAMPLE = "sample1";
    private static final int POSITION = 1000;



    // return a Pair of left reads and right reads
    private Pair<List<GATKRead>, List<GATKRead>> makeReadPairs(final SAMFileHeader header, final String baseName, final boolean leftOverlap, final boolean rightOverlap, final int numPairs) {
        // nthe bases don't matter because we assign likelihoods manually
        final byte[] bases = "AAAAAAAAAA".getBytes();
        final byte[] quals = new byte[bases.length];
        Arrays.fill(quals, (byte) 30);
        final int contigIndex = 0;

        final List<GATKRead> leftReads = new ArrayList<>();
        final List<GATKRead> rightReads = new ArrayList<>();
        for (int n = 0; n < numPairs; n++) {
            final String name = baseName + n;

            // the attention to start position isn't actually necessary because "overlap" is expressed via the likelihoods
            // in that lack of overlap is reflected in uninformative likelihoods that are equal between ref and alt
            final int leftStart = leftOverlap ? POSITION - bases.length + 2 : 1;
            final int rightStart = rightOverlap ? POSITION - 2 : POSITION + 2;

            final GATKRead leftRead = ArtificialReadUtils.createArtificialRead(header, name, contigIndex, leftStart, bases, quals);
            final GATKRead rightRead = ArtificialReadUtils.createArtificialRead(header, name, contigIndex, rightStart, bases, quals);
            leftRead.setIsFirstOfPair();
            rightRead.setIsSecondOfPair();
            leftReads.add(leftRead);
            rightReads.add(rightRead);
        }
        return ImmutablePair.of(leftReads, rightReads);
    }

    @Test
    public void testDescription(){
        Assert.assertEquals(new FragmentDepthPerAlleleBySample().getKeyNames(), Collections.singletonList(GATKVCFConstants.FRAGMENT_ALLELE_DEPTHS));
        Assert.assertEquals(new FragmentDepthPerAlleleBySample().getDescriptions(), Collections.singletonList(VCFStandardHeaderLines.getFormatLine(GATKVCFConstants.FRAGMENT_ALLELE_DEPTHS)));
    }

    // provide # of pairs where left/right reads are alt/neither, neither/alt, alt/alt, ref/neither, neither/ref, ref/ref
    @DataProvider(name = "ReadDepthData")
    public Object[][] readDepthData() {
        return new Object[][]{
                {10, 10, 10, 10, 10, 10},
                {0, 0, 0, 0, 0, 0},
                {1, 2, 3, 4, 5, 6},
                {1, 2, 3, 1, 2, 3},
                {1, 1, 2, 3, 5, 8},
                {2, 3, 5, 7, 11, 13}
        };
    }

    @Test(dataProvider = "ReadDepthData")
    public void testUsingReads(final int leftAlt, final int rightAlt, final int bothAlt, final int leftRef, final int rightRef, final int bothRef){
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 100000000);
        final int[] expectedFAD = {leftRef + rightRef + bothRef, leftAlt + rightAlt + bothAlt}; // regular AD would multiply bothRef and bothAlt by 2

        final List<GATKRead> altReads = new ArrayList<>(leftAlt + rightAlt + 2 * bothAlt);
        final List<GATKRead> refReads = new ArrayList<>(leftRef + rightRef + 2 * bothRef);
        final List<GATKRead> uninformativeReads = new ArrayList<>(leftAlt + rightAlt + leftRef + rightRef);

        final Pair<List<GATKRead>, List<GATKRead>> onlyLeftSupportsAlt = makeReadPairs(header, "leftAlt", true, false, leftAlt);
        altReads.addAll(onlyLeftSupportsAlt.getLeft());
        uninformativeReads.addAll(onlyLeftSupportsAlt.getRight());

        final Pair<List<GATKRead>, List<GATKRead>> onlyRightSupportsAlt = makeReadPairs(header, "rightAlt", false, true, rightAlt);
        altReads.addAll(onlyRightSupportsAlt.getRight());
        uninformativeReads.addAll(onlyRightSupportsAlt.getLeft());

        final Pair<List<GATKRead>, List<GATKRead>> bothSupportAlt = makeReadPairs(header, "bothAlt", true, true, bothAlt);
        altReads.addAll(bothSupportAlt.getLeft());
        altReads.addAll(bothSupportAlt.getRight());

        final Pair<List<GATKRead>, List<GATKRead>> onlyLeftSupportsRef = makeReadPairs(header, "leftRef", true, false, leftRef);
        refReads.addAll(onlyLeftSupportsRef.getLeft());
        uninformativeReads.addAll(onlyLeftSupportsAlt.getRight());

        final Pair<List<GATKRead>, List<GATKRead>> onlyRightSupportsRef = makeReadPairs(header, "rightRef", false, true, rightRef);
        refReads.addAll(onlyRightSupportsRef.getRight());
        uninformativeReads.addAll(onlyRightSupportsAlt.getLeft());

        final Pair<List<GATKRead>, List<GATKRead>> bothSupportRef = makeReadPairs(header, "bothRef", true, true, bothRef);
        refReads.addAll(bothSupportRef.getLeft());
        refReads.addAll(bothSupportRef.getRight());

        final AlleleLikelihoods<GATKRead, Allele> readLikelihoods = ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, uninformativeReads,
                -10.0, -10.0, -1.0, REF, ALT);

        final AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods = readLikelihoods.groupEvidence(GATKRead::getName, Fragment::createAndAvoidFailure);


        final Genotype genotype = new GenotypeBuilder(SAMPLE, ALLELES).DP(leftAlt + rightAlt + 2 * bothAlt + leftRef + rightRef + 2 * bothRef).make();
        final VariantContext vc = new VariantContextBuilder("test", header.getSequence(0).getContig(), POSITION, POSITION, ALLELES).genotypes(Arrays.asList(genotype)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(genotype);
        new FragmentDepthPerAlleleBySample().annotate(null, null, vc, genotype, gb, readLikelihoods, fragmentLikelihoods, null);
        final int[] fad = VariantContextGetters.getAttributeAsIntArray(gb.make(), GATKVCFConstants.FRAGMENT_ALLELE_DEPTHS, () -> new int[] {0,0}, 0);
        Assert.assertEquals(fad, expectedFAD);

    }

}
