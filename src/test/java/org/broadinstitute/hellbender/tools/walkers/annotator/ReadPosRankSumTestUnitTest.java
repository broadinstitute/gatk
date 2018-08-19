package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class ReadPosRankSumTestUnitTest extends GATKBaseTest {
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
                ArtificialAnnotationUtils.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

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

    //Basic aligned read
    private GATKRead allMatch;

    //Read with insertion and deletion
    private GATKRead twoIndels;

    //Read with soft clips at start
    private GATKRead softClipStart;

    //Read with hard clips at start
    private GATKRead hardClipStart;

    //Read with low quality tail
    private GATKRead lowQualTail;

    //Read with low quality tail, partially soft clipped
    private GATKRead lowQualClippedTail;

    //Read with low quality bases at start
    private GATKRead lowQualStart;

    //Read with low quality bases, partially soft clipped at both ends
    private GATKRead lowQualBothEnds;

    @BeforeClass
    public void init() {
        List<CigarElement> cigarElements_allMatch = new LinkedList<>();
        cigarElements_allMatch.add(new CigarElement(151, CigarOperator.M));
        allMatch = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_allMatch));

        List<CigarElement> cigarElements_2indels = new LinkedList<>();
        cigarElements_2indels.add(new CigarElement(66, CigarOperator.M));
        cigarElements_2indels.add(new CigarElement(10, CigarOperator.I));
        cigarElements_2indels.add(new CigarElement(7, CigarOperator.M));
        cigarElements_2indels.add(new CigarElement(10, CigarOperator.D));
        cigarElements_2indels.add(new CigarElement(68, CigarOperator.M));
        twoIndels = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_2indels));

        List<CigarElement> cigarElements_softClipStart = new LinkedList<>();
        cigarElements_softClipStart.add(new CigarElement(17, CigarOperator.S));
        cigarElements_softClipStart.add(new CigarElement(134, CigarOperator.M));
        softClipStart = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_softClipStart));

        List<CigarElement> cigarElements_hardClipStart = new LinkedList<>();
        cigarElements_hardClipStart.add(new CigarElement(17, CigarOperator.H));
        cigarElements_hardClipStart.add(new CigarElement(134, CigarOperator.M));
        hardClipStart = ArtificialReadUtils.createArtificialRead(new Cigar(cigarElements_hardClipStart));


        final byte [] bases_lowQualTail = {'A', 'C', 'T', 'G', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals_lowQualTail = {30, 15, 25, 30, 2, 2, 2, 2, 2, 2};
        lowQualTail = ArtificialReadUtils.createArtificialRead(bases_lowQualTail, quals_lowQualTail, "10M");

        final byte [] bases_lowQualClippedTail = {'A', 'C', 'T', 'G', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals_lowQualClippedTail = {30, 15, 25, 30, 2, 2, 2, 2, 2, 2};
        lowQualClippedTail = ArtificialReadUtils.createArtificialRead(bases_lowQualClippedTail, quals_lowQualClippedTail, "8M2S");

        final byte [] bases_lowQualStart = {'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'T', 'G'};
        final byte [] quals_lowQualStart = {2, 2, 2, 2, 2, 2, 30, 15, 25, 30};
        lowQualStart = ArtificialReadUtils.createArtificialRead(bases_lowQualStart, quals_lowQualStart, "10M");

        final byte [] bases_lowQualBothEnds = {'A', 'A', 'A', 'A', 'A', 'A', 'G', 'C', 'T', 'G', 'A', 'A', 'A', 'A', 'A', 'A'};
        final byte [] quals_lowQualBothEnds = { 2, 2, 2, 2, 2, 2, 30, 15, 25, 30, 2, 2, 2, 2, 2, 2};
        lowQualBothEnds = ArtificialReadUtils.createArtificialRead(bases_lowQualBothEnds, quals_lowQualBothEnds, "2S12M2S");

    }

    @DataProvider(name = "makeGetFinalVariantReadPositionTestReads")
    public Object[][] makeFinalPosTestReads() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 10, 10});
        tests.add(new Object[] {allMatch, 140, 10});
        tests.add(new Object[] {twoIndels, 10, 10});
        tests.add(new Object[] {twoIndels, 140, 10});
        tests.add(new Object[] {hardClipStart, 20, 20});
        tests.add(new Object[] {hardClipStart, 110, 6});  //this is what the code produces as-is
        tests.add(new Object[] {softClipStart, 10, 10});
        tests.add(new Object[] {softClipStart, 140, 10});
        tests.add(new Object[] {lowQualTail, 3, 0});
        tests.add(new Object[] {lowQualTail, 2, 2});
        tests.add(new Object[] {lowQualClippedTail, 3, 0});
        tests.add(new Object[] {lowQualClippedTail, 2, 2});
        tests.add(new Object[] {lowQualStart, 7, -4});   //this is what the code produces as-is, but should be 1
        tests.add(new Object[] {lowQualStart, 8, -5}); //this is what the code produces as-is, but should be 1
        tests.add(new Object[] {lowQualBothEnds, 7, -4});   //this is what the code produces as-is, but should be 1
        tests.add(new Object[] {lowQualBothEnds, 8, -5});    //this is what the code produces as-is, but should be 1
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "makeGetFinalVariantReadPositionTestReads")
    public void testGetFinalVariantReadPosition(GATKRead read, int variantPosition, int expected) throws Exception {
        Assert.assertEquals(ReadPosRankSumTest.getFinalVariantReadPosition(read, variantPosition), expected);
    }

    @DataProvider(name = "getNumClippedBasesAtStartTestReads")
    public Object[][] numClipStart() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 0});
        tests.add(new Object[] {twoIndels, 0});
        tests.add(new Object[] {softClipStart, 0});
        tests.add(new Object[] {hardClipStart, 17});
        tests.add(new Object[] {lowQualTail, 0});
        tests.add(new Object[] {lowQualClippedTail, 0});
        tests.add(new Object[] {lowQualStart, 6});
        tests.add(new Object[] {lowQualBothEnds, 6});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getNumClippedBasesAtStartTestReads")
    public void testGetNumClippedBasesAtStart(GATKRead read, int expected) throws Exception {
        Assert.assertEquals(ReadPosRankSumTest.getNumClippedBasesAtStart(read),expected);
    }

    @DataProvider(name = "getNumAlignedBasesTestReads")
    public Object[][] numAligned() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 151});
        tests.add(new Object[] {twoIndels, 151});
        tests.add(new Object[] {softClipStart, 151});
        tests.add(new Object[] {hardClipStart, 117});  //This is what the code produces, but it's wrong
        tests.add(new Object[] {lowQualTail, 4});
        tests.add(new Object[] {lowQualClippedTail, 4});
        tests.add(new Object[] {lowQualStart, 4});
        tests.add(new Object[] {lowQualBothEnds, 4});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getNumAlignedBasesTestReads")
    public void testGetNumAlignedBases(GATKRead read, int expected) throws Exception {
        Assert.assertEquals(ReadPosRankSumTest.getNumAlignedBases(read),expected);
    }

    @DataProvider(name = "getNumClippedBasesAtEndTestReads")
    public Object[][] numClipEnd() {
        final List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {allMatch, 0});
        tests.add(new Object[] {twoIndels, 0});
        tests.add(new Object[] {softClipStart, 0});
        tests.add(new Object[] {hardClipStart, 0});
        tests.add(new Object[] {lowQualTail, 6});
        tests.add(new Object[] {lowQualClippedTail, 6});
        tests.add(new Object[] {lowQualStart, 0});
        tests.add(new Object[] {lowQualBothEnds, 6});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getNumClippedBasesAtEndTestReads")
    public void testGetNumClippedBasesAtEnd(GATKRead read, int expected) throws Exception {
        Assert.assertEquals(ReadPosRankSumTest.getNumClippedBasesAtEnd(read), expected);
    }
}
