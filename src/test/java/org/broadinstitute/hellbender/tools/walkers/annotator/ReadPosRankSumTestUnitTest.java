package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_ReadPosRankSumTest;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class ReadPosRankSumTestUnitTest extends BaseTest {
    private final String sample1 = "NA1";
    private final String sample2 = "NA2";
    private static final String CONTIG = "1";

    private static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    private VariantContext makeVC(final long position) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(REF, ALT)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(REF, ALT)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(REF, ALT)).chr(CONTIG).start(position).stop(position).genotypes(testGC).make();
    }

    private GATKRead makeRead(final int start, final int mq) {
        Cigar cigar = TextCigarCodec.decode("10M");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mq);
        read.setPosition(CONTIG, start);
        return read;
    }

    @Test
    public void testReadPos(){
        final InfoFieldAnnotation ann = new ReadPosRankSumTest();
        final String key =  GATKVCFConstants.READ_POS_RANK_SUM_KEY;
        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        final int[] startAlts = {3, 4};
        final int[] startRefs = {1, 2};
        final List<GATKRead> refReads = Arrays.asList(makeRead(startRefs[0], 30), makeRead(startRefs[1], 30));
        final List<GATKRead> altReads = Arrays.asList(makeRead(startAlts[0], 30), makeRead(startAlts[1], 30));
        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key);

        final ReferenceContext ref= null;

        final long position = 5L;  //middle of the read
        final VariantContext vc= makeVC(position);

        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);
        final double zScore = mannWhitneyU.test(new double[]{position - startAlts[0], position - startAlts[1]}, new double[]{position - startRefs[0], position - startRefs[1]}, MannWhitneyU.TestType.FIRST_DOMINATES).getZ();
        final String zScoreStr = String.format("%.3f", zScore);
        Assert.assertEquals(annotate.get(key), zScoreStr);
        
        final long positionEnd = 8L;  //past middle
        final VariantContext vcEnd= makeVC(positionEnd);

        //Note: past the middle of the read we compute the position from the end.
        final Map<String, Object> annotateEnd = ann.annotate(ref, vcEnd, likelihoods);
        final double zScoreEnd = mannWhitneyU.test(new double[]{startAlts[0], startAlts[1]}, new double[]{startRefs[0], startRefs[1]}, MannWhitneyU.TestType.FIRST_DOMINATES).getZ();
        final String zScoreEndStr = String.format("%.3f", zScoreEnd);
        Assert.assertEquals(annotateEnd.get(key), zScoreEndStr);

        final long positionPastEnd = 20L;  //past middle
        final VariantContext vcPastEnd= makeVC(positionPastEnd);

        //Note: past the end of the read, there's nothing
        final Map<String, Object> annotatePastEnd = ann.annotate(ref, vcPastEnd, likelihoods);
        Assert.assertTrue(annotatePastEnd.isEmpty());
    }

    @Test
    public void testReadPos_Raw(){
        final AS_RankSumTest ann= new AS_ReadPosRankSumTest();
        final String key1 = GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY;
        final String key2 = GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY;
        final int[] startAlts = {3, 4};
        final int[] startRefs = {1, 2};
        final int readLength = 10;

        final List<GATKRead> refReads = Arrays.asList(makeRead(startRefs[0], 30), makeRead(startRefs[1], 30));
        final List<GATKRead> altReads = Arrays.asList(makeRead(startAlts[0], 30), makeRead(startAlts[1], 30));
        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key1);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key2);

        final ReferenceContext ref= null;

        final long position = 5L;  //middle of the read
        final VariantContext vc= makeVC(position);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotateNonRaw = ann.annotate(ref, vc, likelihoods);
        final String expected = startAlts[0] + ",1," + startAlts[1] + ",1" + AS_RankSumTest.PRINT_DELIM + startRefs[0] + ",1," + startRefs[1] + ",1";
        Assert.assertEquals(annotateRaw.get(key1), expected);
        Assert.assertEquals(annotateNonRaw.get(key1), expected);
        
        final long positionEnd = 8L;  //past middle
        final VariantContext vcEnd= makeVC(positionEnd);

        //Note: past the middle of the read we compute the position from the end.
        final Map<String, Object> annotateEndRaw = ann.annotateRawData(ref, vcEnd, likelihoods);
        final Map<String, Object> annotateEndNonRaw = ann.annotate(ref, vcEnd, likelihoods);
        final String refS = (startRefs[0]+readLength-positionEnd-1)+ ",1," +(startRefs[1]+readLength-positionEnd-1) + ",1";
        final String altS = (positionEnd-startAlts[1]) + ",1," + (positionEnd-startAlts[0]) + ",1";
        Assert.assertEquals(annotateEndRaw.get(key1), refS + AS_RankSumTest.PRINT_DELIM + altS );
        Assert.assertEquals(annotateEndNonRaw.get(key1), refS + AS_RankSumTest.PRINT_DELIM + altS );

        final long positionPastEnd = 20L;  //past middle
        final VariantContext vcPastEnd= makeVC(positionPastEnd);

        //Note: past the end of the read, there's nothing
        final Map<String, Object> annotatePastEndRaw = ann.annotateRawData(ref, vcPastEnd, likelihoods);
        final Map<String, Object> annotatePastEndNonRaw = ann.annotate(ref, vcPastEnd, likelihoods);
        Assert.assertTrue(annotatePastEndRaw.isEmpty());
        Assert.assertTrue(annotatePastEndNonRaw.isEmpty());
    }

    @DataProvider(name = "dataIsUsableRead")
    private Object[][] dataIsUsableRead(){
        return new Object[][]{
                {"20M6D2M", 10, 1, 0, true},
                {"20M6D2M", 10, 1, 27, true},
                {"20M6D2M", 10, 1, 29, false},
                {"1I20M1S", 10, 1, 0, true},
                {"1I20M1S", 0, 1, 0, false},
                {"1I20M1S", QualityUtils.MAPPING_QUALITY_UNAVAILABLE, 1, 0, false},
                {"1I20M1S", 10, 1, 22, false},
                {"1I20M1S", 10, 21, 42, false},
                {"1I20M1H", 10, 1, 21, false},
        };
    }

    @Test(dataProvider = "dataIsUsableRead")
    public void testIsUsableRead(final String cigarString, final int mappingQuality, final int start, final int refLoc, final boolean isUsable ) {
        final ReadPosRankSumTest readPosRankSumTest = new ReadPosRankSumTest();
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mappingQuality);
        read.setPosition(CONTIG, start);
        Assert.assertEquals(readPosRankSumTest.isUsableRead(read, refLoc), isUsable);
    }
}
