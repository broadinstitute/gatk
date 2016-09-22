package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public final class BaseQualitySumPerAlleleBySampleUnitTest {
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);

    private static final String SAMPLE_1 = "sample1";

    private GATKRead makeRead(final int baseQual) {
        return AnnotationArtificialData.makeRead(baseQual, 50);
    }

    @Test
    public void testUsableRead() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(5, 1, 10000);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, 76);
        read.setMappingQuality(60);
        Assert.assertTrue(BaseQualitySumPerAlleleBySample.isUsableRead(read));

        read.setMappingQuality(0);
        Assert.assertFalse(BaseQualitySumPerAlleleBySample.isUsableRead(read));

        read.setMappingQuality(QualityUtils.MAPPING_QUALITY_UNAVAILABLE);
        Assert.assertFalse(BaseQualitySumPerAlleleBySample.isUsableRead(read));
    }

    @Test
    public void testDescriptions() {
        Assert.assertEquals(new BaseQualitySumPerAlleleBySample().getKeyNames(), Collections.singletonList(GATKVCFConstants.QUALITY_SCORE_SUM_KEY), "annots");
        Assert.assertEquals(new BaseQualitySumPerAlleleBySample().getDescriptions(), Collections.singletonList(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.QUALITY_SCORE_SUM_KEY)));
    }

    @Test
    public void testUsingReads(){
        final int refDepth = 20;
        final int altDepth = 17;

        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE_1, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;
        final byte baseQual = 23;

        final List<GATKRead> refReads = Collections.nCopies(refDepth, makeRead(baseQual));
        final List<GATKRead> altReads = Collections.nCopies(altDepth, makeRead(baseQual));
        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.makeLikelihoods(SAMPLE_1, refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new BaseQualitySumPerAlleleBySample().annotate(null, vc, gAC, gb, likelihoods);
        final Integer[] quals = (Integer[]) gb.make().getAnyAttribute(GATKVCFConstants.QUALITY_SCORE_SUM_KEY);
        final Integer[] extectedAD = {refDepth * baseQual, altDepth * baseQual};
        Assert.assertEquals(quals, extectedAD);

        //now test a no-op
        final GenotypeBuilder gb1 = new GenotypeBuilder(gAC);
        new BaseQualitySumPerAlleleBySample().annotate(null, vc, null, gb1, likelihoods);  //null genotype
        Assert.assertFalse(gb1.make().hasAD());
    }
}
