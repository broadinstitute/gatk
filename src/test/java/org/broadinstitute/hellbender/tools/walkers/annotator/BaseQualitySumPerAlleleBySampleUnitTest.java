package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
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
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final Allele A = Allele.create("A", true);
        final Allele C = Allele.create("C");

        final List<Allele> AC = Arrays.asList(A, C);
        final int readDepthRef = 20;
        final int readDepthAlt = 17;

        final String sample1 = "sample1";
        final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(sample1, AC).DP(dpDepth).make();

        final double log10PError = -5;
        final byte baseQual = 23;
        final int readLen = 10;

        for (int i = 0; i < readDepthAlt; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(readLen + "M"), "read1_" + i);
            read.setBaseQualities(Utils.dupBytes(baseQual, readLen));
            read.setMappingQuality(20);
            map.add(read, A, -10.0);
            map.add(read, C, -1.0);      //try to fool it - add another likelihood to same read
        }
        for (int i = 0; i < readDepthRef; i++) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(readLen + "M"), "read2_" + i);
            read.setBaseQualities(Utils.dupBytes(baseQual, readLen));
            read.setMappingQuality(20);
            map.add(read, A, -1.0);
            map.add(read, C, -100.0);  //try to fool it - add another likelihood to same read
        }

        //throw in one non-informative read
        final GATKRead badRead = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(readLen + "M"), "read3");
        badRead.setBaseQualities(Utils.dupBytes(baseQual, readLen));
        badRead.setMappingQuality(20);
        map.add(badRead, A, -1.0);
        map.add(badRead, C, -1.1); //maybe it's ref, maybe it's alt, too close to call -> not informative

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, AC).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new BaseQualitySumPerAlleleBySample().annotate(null, vc, gAC, gb, map);
        final Integer[] quals = (Integer[]) gb.make().getAnyAttribute(GATKVCFConstants.QUALITY_SCORE_SUM_KEY);
        final Integer[] extectedAD = {readDepthRef * baseQual, readDepthAlt * baseQual};
        Assert.assertEquals(quals, extectedAD);

        //now test a no-op
        final GenotypeBuilder gb1 = new GenotypeBuilder(gAC);
        new BaseQualitySumPerAlleleBySample().annotate(null, vc, null, gb1, map);  //null genotype
        Assert.assertFalse(gb1.make().hasAD());
    }
}
