package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ClippedReadsCountTest extends GATKBaseTest {
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final String SAMPLE = "sample1";

    private GATKRead makeRead(final byte readBase) {
        final int readLength = 10;
        return ArtificialAnnotationUtils.makeRead(readBase, (byte)30, 50);
    }
    private GATKRead makeSoftClippedRead(final byte readBase, final boolean leftClip){
        return ArtificialAnnotationUtils.makeSoftClippedRead(readBase,(byte)10, 50, leftClip);
    }

    @Test
    public void testDescription(){
        String[] constants = {GATKVCFConstants.SOFT_CLIP_LEFT_COUNT_KEY, GATKVCFConstants.SOFT_CLIP_RIGHT_COUNT_KEY};
        VCFFormatHeaderLine[] hlines    = {GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.SOFT_CLIP_LEFT_COUNT_KEY),
                              GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.SOFT_CLIP_RIGHT_COUNT_KEY)};
        Assert.assertEquals(new ClippedReadsCount().getKeyNames(), new ArrayList<>(Arrays.asList(constants)));
        Assert.assertEquals(new ClippedReadsCount().getDescriptions(), new ArrayList<>(Arrays.asList(hlines)));
    }

    @DataProvider(name = "SoftClippingData")
    public Object[][][] readDepthData() {
        return new Object[][][]{{
            {20,0}, {0,17}},
            {{0,0}, {0,0}}
        };
    }

    @Test(dataProvider = "SoftClippingData")
    public void testUsingReads(final Object[] refDepth, final Object[] altDepth){
        final int[] expectedLeftCount = {(int)refDepth[0], (int)altDepth[0]};
        final int[] expectedRightCount = {(int)refDepth[1], (int)altDepth[1]};


        final int dpDepth = 50; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> refReads = IntStream.range(0, dpDepth - expectedLeftCount[0] - expectedRightCount[0]).mapToObj(i -> makeRead(REF.getBases()[0])).collect(Collectors.toList());
        final List<GATKRead> refClippedReads = IntStream.range(0, expectedLeftCount[0]).mapToObj(i -> makeSoftClippedRead(REF.getBases()[0], true)).collect(Collectors.toList());
        final List<GATKRead> altClippedReads = IntStream.range(0, expectedRightCount[1]).mapToObj(i -> makeSoftClippedRead(ALT.getBases()[0], false)).collect(Collectors.toList());
        refReads.addAll(refClippedReads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altClippedReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "1", 10005, 10005, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new ClippedReadsCount().annotate(new ReferenceContext(new ReferenceFileSource(Path.of(b37_reference_20_21)),
                new SimpleInterval("20", 10005, 10005)), vc, gAC, gb, likelihoods);
        final int[] scl = (int [])gb.make().getExtendedAttribute("SCL");
        final int[] scr = (int [])gb.make().getExtendedAttribute("SCR");
        Assert.assertEquals(scl, expectedLeftCount);
        Assert.assertEquals(scr, expectedRightCount);

    }
}