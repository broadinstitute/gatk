package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class AS_ReadPosRankSumTestUnitTest extends ReducibleAnnotationBaseTest {
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
        final AS_ReadPosRankSumTest as_readPosRankSumTest = new AS_ReadPosRankSumTest();
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mappingQuality);
        read.setPosition("1", start);
        Assert.assertEquals(as_readPosRankSumTest.isUsableRead(read, refLoc), isUsable);
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
                ArtificialAnnotationUtils.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key2);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key2);

        final ReferenceContext ref= null;

        final long position = 5L;  //middle of the read
        final VariantContext vc= makeVC(position);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotateNonRaw = ann.annotate(ref, vc, likelihoods);

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        MannWhitneyU.Result expectedAlt = mannWhitneyU.test(new double[]{1.0, 2.0},new double[]{3.0, 4.0}, MannWhitneyU.TestType.FIRST_DOMINATES);
        String firstExpected = "|"+String.format("%.1f",Math.round(Math.floor((expectedAlt.getZ() )/0.1))*0.1)+",1";
        String firstNonRawExpected = String.format("%.3f",expectedAlt.getZ());

        Assert.assertEquals(annotateRaw.get(key1), firstExpected);
        Assert.assertEquals(annotateNonRaw.get(key2), firstNonRawExpected);

        final long positionEnd = 8L;  //past middle
        final VariantContext vcEnd= makeVC(positionEnd);

        //Note: past the middle of the read we compute the position from the end.
        final Map<String, Object> annotateEndRaw = ann.annotateRawData(ref, vcEnd, likelihoods);
        final Map<String, Object> annotateEndNonRaw = ann.annotate(ref, vcEnd, likelihoods);

        expectedAlt = mannWhitneyU.test(new double[]{(positionEnd-startAlts[1]),(positionEnd-startAlts[0]) },new double[]{(startRefs[0]+readLength-positionEnd-1), (startRefs[1]+readLength-positionEnd-1)}, MannWhitneyU.TestType.FIRST_DOMINATES);
        String secondExpected = "|"+String.format("%.1f",Math.round(Math.floor((expectedAlt.getZ() )/0.1))*0.1)+",1";
        String secondNonRawExpected = String.format("%.3f",expectedAlt.getZ());

        Assert.assertEquals(annotateEndRaw.get(key1), secondExpected );
        Assert.assertEquals(annotateEndNonRaw.get(key2), secondNonRawExpected );

        final long positionPastEnd = 20L;  //past middle
        final VariantContext vcPastEnd= makeVC(positionPastEnd);

        //Note: past the end of the read, there's nothing
        final Map<String, Object> annotatePastEndRaw = ann.annotateRawData(ref, vcPastEnd, likelihoods);
        final Map<String, Object> annotatePastEndNonRaw = ann.annotate(ref, vcPastEnd, likelihoods);
        Assert.assertEquals(annotatePastEndRaw.get(key1), "|");
        Assert.assertEquals(annotatePastEndNonRaw.get(key2), null);
    }
    
    @Override
    protected List<Annotation> getAnnotationsToUse() {
        return Collections.singletonList(new AS_ReadPosRankSumTest());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY;
    }

    @Override
    protected String getKey() {
        return GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY;
    }

}
