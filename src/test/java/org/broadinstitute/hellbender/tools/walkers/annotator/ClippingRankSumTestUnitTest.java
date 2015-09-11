package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;

public final class ClippingRankSumTestUnitTest {

    private final String sample1 = "NA1";
    private final String sample2 = "NA2";

    private VariantContext makeVC( final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private GATKRead makeRead(final int hardClip, final int mq) {
        Cigar cigar = hardClip == 0 ? TextCigarCodec.decode("10M") : TextCigarCodec.decode("10M" + hardClip + "H");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mq);
        return read;
    }
    @Test
    public void testClipping(){
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final int[] hardAlts = {1, 2};
        final int[] hardRefs = {10, 0};
        final GATKRead read1 = makeRead(hardAlts[0], 30);
        final GATKRead read2 = makeRead(hardAlts[1], 30);
        final GATKRead read3 = makeRead(hardRefs[0], 30);
        final GATKRead read4 = makeRead(hardRefs[1], 30);
        map.add(read1, alleleAlt, -1.0);
        map.add(read1, alleleRef, -100.0);

        map.add(read2, alleleAlt, -1.0);
        map.add(read2, alleleRef, -100.0);

        map.add(read3, alleleAlt, -100.0);
        map.add(read3, alleleRef, -1.0);

        map.add(read4, alleleAlt, -100.0);
        map.add(read4, alleleRef, -1.0);

        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap = Collections.singletonMap(sample1, map);


        final ReferenceContext ref= null;
        final VariantContext vc= makeVC(alleleRef, alleleAlt);
        final InfoFieldAnnotation ann = new ClippingRankSumTest();

        final Map<String, Object> annotate = ann.annotate(ref, vc, stratifiedPerReadAlleleLikelihoodMap);

        final double val= MannWhitneyU.runOneSidedTest(false, Arrays.asList(hardAlts[0], hardAlts[1]),
                                                              Arrays.asList(hardRefs[0], hardRefs[1])).getLeft();
        final String valStr= String.format("%.3f", val);
        Assert.assertEquals(annotate.get(GATKVCFConstants.CLIPPING_RANK_SUM_KEY), valStr);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), GATKVCFConstants.CLIPPING_RANK_SUM_KEY);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), GATKVCFConstants.CLIPPING_RANK_SUM_KEY);


    }
}
