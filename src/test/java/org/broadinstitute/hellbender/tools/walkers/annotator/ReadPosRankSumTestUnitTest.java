package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;

public final class ReadPosRankSumTestUnitTest extends BaseTest {
    private final String sample1 = "NA1";
    private final String sample2 = "NA2";

    private VariantContext makeVC( final String contig, final long position, final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr(contig).start(position).stop(position).genotypes(testGC).make();
    }

    private GATKRead makeRead(final String contig, final int start, final int mq) {
        Cigar cigar = TextCigarCodec.decode("10M");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mq);
        read.setPosition(contig, start);
        return read;
    }
    @Test
    public void testReadPos(){
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final String contig = "1";

        final Allele alleleRef = Allele.create("T", true);
        final Allele alleleAlt = Allele.create("A", false);

        final int[] startAlts = {3, 4};
        final int[] startRefs = {1, 2};
        final GATKRead read1 = makeRead(contig, startAlts[0],  30);
        final GATKRead read2 = makeRead(contig, startAlts[1], 30);
        final GATKRead read3 = makeRead(contig, startRefs[0], 30);
        final GATKRead read4 = makeRead(contig, startRefs[1], 30);
        map.add(read1, alleleAlt, -1.0);
        map.add(read1, alleleRef, -100.0);

        map.add(read2, alleleAlt, -1.0);
        map.add(read2, alleleRef, -100.0);

        map.add(read3, alleleAlt, -100.0);
        map.add(read3, alleleRef, -1.0);

        map.add(read4, alleleAlt, -100.0);
        map.add(read4, alleleRef, -1.0);

        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap = Collections.singletonMap(sample1, map);

        final InfoFieldAnnotation ann = new ReadPosRankSumTest();
        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), GATKVCFConstants.READ_POS_RANK_SUM_KEY);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), GATKVCFConstants.READ_POS_RANK_SUM_KEY);


        final ReferenceContext ref= null;

        final long position = 5L;  //middle of the read
        final VariantContext vc= makeVC(contig, position, alleleRef, alleleAlt);

        final Map<String, Object> annotate = ann.annotate(ref, vc, stratifiedPerReadAlleleLikelihoodMap);
        final double val= MannWhitneyU.runOneSidedTest(false,
                Arrays.asList(position - startAlts[0], position - startAlts[1]),
                Arrays.asList(position - startRefs[0], position - startRefs[1])).getLeft();
        final String valStr= String.format("%.3f", val);
        Assert.assertEquals(annotate.get(GATKVCFConstants.READ_POS_RANK_SUM_KEY), valStr);


        final long positionEnd = 8L;  //past middle
        final VariantContext vcEnd= makeVC(contig, positionEnd, alleleRef, alleleAlt);

        //Note: past the middle of the read we compute the position from the end.
        final Map<String, Object> annotateEnd = ann.annotate(ref, vcEnd, stratifiedPerReadAlleleLikelihoodMap);
        final double valEnd= MannWhitneyU.runOneSidedTest(false,
                Arrays.asList(startAlts[0], startAlts[1]),
                Arrays.asList(startRefs[0], startRefs[1])).getLeft();
        final String valStrEnd= String.format("%.3f", valEnd);
        Assert.assertEquals(annotateEnd.get(GATKVCFConstants.READ_POS_RANK_SUM_KEY), valStrEnd);

        final long positionPastEnd = 20L;  //past middle
        final VariantContext vcPastEnd= makeVC(contig, positionPastEnd, alleleRef, alleleAlt);

        //Note: past the end of the read, there's nothing
        final Map<String, Object> annotatePastEnd = ann.annotate(ref, vcPastEnd, stratifiedPerReadAlleleLikelihoodMap);
        Assert.assertNull(annotatePastEnd);

    }
}
