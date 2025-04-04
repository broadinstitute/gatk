package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ReblockGVCF;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.BiPredicate;

import static htsjdk.variant.vcf.VCFConstants.MAX_GENOTYPE_QUAL;
import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.makeHomRef;
import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.makeHomRefAlt;
import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.makeDeletion;
import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.makeNonRef;
import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.makeVariantContext;

public class GVCFWriterUnitTest extends GATKBaseTest {

    private static final String CHR1 = "1";
    private static final String CHR2 = "2";
    private static final List<Allele> ALLELES = ImmutableList.of(Allele.REF_G, Allele.NON_REF_ALLELE);
    private static final String SAMPLE_NAME = "XXYYZZ";


    private static final List<Number> standardPartition = ImmutableList.of(1, 10, 20);
    private static final List<Number> highConfLowConf = ImmutableList.of(20,100);


    @Test
    public void testHeaderWriting() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);
        writer.writeHeader(new VCFHeader());
        Assert.assertTrue(mockWriter.headerSet);
        Assert.assertTrue(mockWriter.headerWritten);
    }

    @Test
    public void testHeaderSetting(){
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);
        writer.setHeader(new VCFHeader());
        Assert.assertTrue(mockWriter.headerSet);
        Assert.assertFalse(mockWriter.headerWritten);
    }

    @Test
    public void testClose() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);
        writer.close();
        Assert.assertTrue(mockWriter.closed);
    }

    @Test
    public void testCloseEmitsLastVariant() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        Assert.assertEquals(mockWriter.emitted.size(), 0);

        writer.close();
        Assert.assertTrue(mockWriter.closed);
        Assert.assertEquals(mockWriter.emitted.size(), 1);
    }

    @Test
    public void testCloseDoesntEmitsLastVariantWhenNonRef() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeNonRef(CHR1, 1));
        Assert.assertEquals(mockWriter.emitted.size(), 1);

        writer.close();
        Assert.assertTrue(mockWriter.closed);
        Assert.assertEquals(mockWriter.emitted.size(), 1);
    }

    @Test
    public void testCrossingContigBoundaryRef() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        writer.add(makeHomRef(2));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef(CHR2, 3, 0));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 2, false);

        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(1), CHR2, 3, 3, false);
    }

    @Test
    public void testCrossingContigBoundaryToLowerPositionsRef() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(30));
        writer.add(makeHomRef(31));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef(CHR2, 10, 0));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 30, 31, false);
        writer.add(makeNonRef(CHR2, 11));
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(1), CHR2, 10, 10, false);
        assertGoodVC(mockWriter.emitted.get(2), CHR2, 11, 11, true);
    }

    @Test
    public void testCrossingContigBoundaryFromNonRefToLowerPositionsRef() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeNonRef(CHR1, 20));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        writer.add(makeHomRef(CHR2, 10, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 20, 20, true);
        writer.add(makeNonRef(CHR2, 11));
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(1), CHR2, 10, 10, false);
        assertGoodVC(mockWriter.emitted.get(2), CHR2, 11, 11, true);
    }

    @Test
    public void testCrossingContigBoundaryNonRef() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        writer.add(makeHomRef(2));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeNonRef(CHR2, 3));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), CHR2, 3, 3, true);
    }

    @Test
    public void testCrossingContigBoundaryNonRefThenNonRef() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeNonRef(CHR1, 1));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        writer.add(makeNonRef(CHR2, 1));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 1, true);
        assertGoodVC(mockWriter.emitted.get(1), CHR2, 1, 1, true);
    }

    @SuppressWarnings("unchecked")
    private static void assertGoodVCwithPPs(final VariantContext vc, final String contig, final int start, final int stop, final boolean nonRef) {
        final Genotype g = vc.getGenotype(SAMPLE_NAME);
        Assert.assertTrue(g.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        final List<Integer> PPs = (List<Integer>)vc.getGenotype(0).getAnyAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY);
        if (!nonRef) {
            Assert.assertTrue(PPs.size() == 3);
            Assert.assertTrue(g.hasGQ());
            Assert.assertTrue(PPs.get(1) == g.getGQ());
        }
        assertGoodVC(vc, contig, start, stop, nonRef);
    }

    private static void assertGoodVC(final VariantContext vc, final String contig, final int start, final int stop, final boolean nonRef) {
        Assert.assertEquals(vc.getContig(), contig);
        Assert.assertEquals(vc.getStart(), start);
        Assert.assertEquals(vc.getEnd(), stop);
        if ( nonRef ) {
            Assert.assertNotEquals(vc.getAlternateAllele(0), Allele.NON_REF_ALLELE);
        } else {
            Assert.assertEquals(vc.getNAlleles(), 2);
            Assert.assertEquals(vc.getAlternateAllele(0), Allele.NON_REF_ALLELE);
            Assert.assertEquals(vc.getAttributeAsInt(VCFConstants.END_KEY, -1), stop);
            Assert.assertTrue(vc.hasGenotypes());
            Assert.assertTrue(vc.hasGenotype(SAMPLE_NAME));
            Assert.assertEquals(vc.getGenotypes().size(), 1);
            final Genotype g = vc.getGenotype(SAMPLE_NAME);
            Assert.assertFalse(g.hasAD());
            Assert.assertTrue(g.hasLikelihoods());
            Assert.assertTrue(g.hasPL());
            Assert.assertEquals(g.getPL().length, 3);
            Assert.assertTrue(g.hasDP());
            Assert.assertTrue(g.hasGQ());
        }
    }

    @Test
    public void testVariantForcesNonRef() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        writer.add(makeHomRef(2));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeNonRef(CHR1, 3));
        writer.add(makeHomRef(4));
        writer.add(makeHomRef(5));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), CHR1, 3, 3, true);
        writer.close();
        assertGoodVC(mockWriter.emitted.get(2), CHR1, 4, 5, false);
    }

    @Test
    public void testEmittingTwoBands() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        writer.add(makeHomRef(2));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef(3, 50));
        writer.add(makeHomRef(4, 50));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), CHR1, 3, 4, false);
    }


    @Test
    public void testBandingUsingPP() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        int[] PPs1 = {0,63,128};
        int[] PPs2 = {0,67,145};
        writer.add(makeVariantContext(new VariantContextBuilder("test", CHR1, 10000, 10000,
                ALLELES), Arrays.asList(Allele.REF_G, Allele.REF_G), 2, PPs1));
        writer.add(makeVariantContext(new VariantContextBuilder("test", CHR1, 10001, 10001, ALLELES), Arrays.asList(Allele.REF_G, Allele.REF_G), 21, PPs2));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVCwithPPs(mockWriter.emitted.get(0), CHR1, 10000, 10001, false);
    }


    @Test
    public void testNonContiguousBlocks() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        writer.add(makeHomRef(2));
        writer.add(makeHomRef(10));
        writer.add(makeHomRef(11));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), CHR1, 10, 11, false);
    }

    @Test
    public void testInputBlocks() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, highConfLowConf);

        writer.add(makeHomRef("20", 1, 16, 600));
        writer.add(makeHomRef("20", 601, 0, 620));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 620, false);
    }

    @Test
    public void testDeletion() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        writer.add(makeHomRef(2));
        writer.add(makeDeletion(3, 3));
        writer.add(makeHomRef(4));
        writer.add(makeHomRef(5));
        writer.add(makeHomRef(6));
        writer.add(makeHomRef(7));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), CHR1, 3, 5, true);
        assertGoodVC(mockWriter.emitted.get(2), CHR1, 6, 7, false);
    }

    @Test
    public void testHomRefAlt() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition);

        writer.add(makeHomRef(1));
        writer.add(makeHomRef(2));
        writer.add(makeHomRefAlt(3));
        writer.add(makeHomRef(4));
        writer.add(makeHomRef(5));
        writer.add(makeHomRef(6));
        writer.add(makeHomRef(7));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(0), CHR1, 1, 2, false);
        Assert.assertFalse(mockWriter.emitted.get(1).hasAttribute(VCFConstants.END_KEY));
        Assert.assertFalse(mockWriter.emitted.get(1).hasAttribute("BLOCK_SIZE"));
        assertGoodVC(mockWriter.emitted.get(2), CHR1, 4, 7, false);
    }

    @DataProvider(name = "GoodBandPartitionData")
    public Object[][] makeBandPartitionData() {
        return new Object[][]{
                {Collections.singletonList(1), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,MAX_GENOTYPE_QUAL+1))},
                {Arrays.asList(1, 2, 3), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,2), Range.closedOpen(2,3),Range.closedOpen(3,MAX_GENOTYPE_QUAL+1))},
                {Arrays.asList(1, 10), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,10), Range.closedOpen(10, MAX_GENOTYPE_QUAL+1))},
                {Arrays.asList(1, 10, 30), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,10), Range.closedOpen(10, 30),Range.closedOpen(30,MAX_GENOTYPE_QUAL+1))},
                {Arrays.asList(1, 10, MAX_GENOTYPE_QUAL - 1), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,10), Range.closedOpen(10, MAX_GENOTYPE_QUAL - 1),Range.closedOpen(MAX_GENOTYPE_QUAL - 1,MAX_GENOTYPE_QUAL+1))},
                {Arrays.asList(1, 10, MAX_GENOTYPE_QUAL), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,10), Range.closedOpen(10, MAX_GENOTYPE_QUAL), Range.closedOpen(MAX_GENOTYPE_QUAL, MAX_GENOTYPE_QUAL+1))},
                {Arrays.asList(1, 10, MAX_GENOTYPE_QUAL + 1), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,10), Range.closedOpen(10, MAX_GENOTYPE_QUAL+1))},
                {Collections.singletonList(VCFConstants.MAX_GENOTYPE_QUAL + 1), Arrays.asList(Range.closedOpen(0,MAX_GENOTYPE_QUAL+1))}
        };
    }

    @Test(dataProvider = "GoodBandPartitionData")
    public void testGoodPartitions(final List<Number> partitions, List<Range<Integer>> expected) {
        final GVCFBlockCombiner combiner = new GVCFBlockCombiner(partitions, false);
        Assert.assertEquals(new ArrayList<>(combiner.gqPartitions.asMapOfRanges().values()), expected);
        Assert.assertEquals(new ArrayList<>(combiner.gqPartitions.asMapOfRanges().keySet()), expected);
    }

    @DataProvider(name = "BadBandPartitionData")
    public Object[][] makeBadBandPartitionData() {
        return new Object[][]{
                {null},
                {Collections.emptyList()},
                {Arrays.asList(10, 1, 30)},
                {Arrays.asList(-1, 1)},
                {Arrays.asList(1, null, 10)},
                {Arrays.asList(1, 1, 10)},
                {Arrays.asList(1, 10, MAX_GENOTYPE_QUAL+2)}
        };
    }

    @Test(dataProvider = "BadBandPartitionData", expectedExceptions = IllegalArgumentException.class)
    public void testBadPartitionsThrowException(final List<Number> partitions){
        GVCFBlockCombiner combiner = new GVCFBlockCombiner(partitions, false); // we should explode here
    }

    @Test
    public void testCheckError(){
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter gvcfWriter = new GVCFWriter(mockWriter, standardPartition);
        mockWriter.error = false;
        Assert.assertEquals(gvcfWriter.checkError(), mockWriter.checkError());
        mockWriter.error = true;
        Assert.assertEquals(gvcfWriter.checkError(), mockWriter.checkError());
    }

    @Test
    public void testToVCFHeaderLine() {
        final Range<Integer> band = Range.closedOpen(10,20);
        Assert.assertEquals(GVCFBlockCombiner.rangeToVCFHeaderLine(band).getKey(), "GVCFBlock10-20", "Wrong key for " + band);
        Assert.assertEquals(GVCFBlockCombiner.rangeToVCFHeaderLine(band).getValue(), "minGQ=10(inclusive),maxGQ=20(exclusive)", "Wrong value for" + band);
    }

    @DataProvider(name = "toWriteToDisk")
    public Object[][] toWriteToDisk(){

        return new Object[][]{
                {Arrays.asList(makeHomRef(1, 33), makeHomRef(2, 31), makeHomRef(3, 30)), //merge multiple sites into one band
                        Collections.singletonList(new MinimalData(CHR1, 1, 3, 30))},
                {Arrays.asList(makeHomRef(CHR1, 1, 10), makeHomRef(CHR2, 2, 10)), //don't merge across chromosomes
                        Arrays.asList(new MinimalData(CHR1, 1, 1, 10), new MinimalData(CHR2, 2, 2, 10))},
                {Arrays.asList(makeHomRef(1,30), makeHomRef(2,15),makeHomRef(3,15)), //multiple bands stay distinct
                    Arrays.asList(new MinimalData(CHR1, 1, 1, 30), new MinimalData(CHR1, 2,3,15))},
                {Collections.singletonList(makeDeletion(100, 5)), Collections.singletonList(new MinimalData(CHR1, 100, 104, 50))},
                {Arrays.asList(makeHomRef(1, 15), makeHomRef(2, 16), makeNonRef(CHR1, 3), makeHomRef(4, 15)),
                        Arrays.asList(new MinimalData(CHR1, 1, 2, 15), new MinimalData(CHR1, 3, 3, 30), new MinimalData(CHR1, 4, 4, 15))}

        };
    }

    /**
     * data class to hold a SimpleInterval and a GQ together for easier testing
     */
    private static class MinimalData {
        public final SimpleInterval location;
        public final int gq;

        public MinimalData(String chr, int start, int end, int gq) {
            this.location = new SimpleInterval(chr, start, end);
            this.gq = gq;
        }
    }

    @Test(dataProvider = "toWriteToDisk")
    public void writeGVCFToDisk(List<VariantContext> variants, List<MinimalData> expected) {
        final List<Number> gqPartitions = Arrays.asList(1, 10, 30);
        final File outputFile =  createTempFile("generated", ".g.vcf");

        try (VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(outputFile.toPath(), null, false);
             GVCFWriter gvcfWriter = new GVCFWriter(writer, gqPartitions))
        {
            gvcfWriter.writeHeader(getMinimalVCFHeader());
            variants.forEach(gvcfWriter::add);
        }
        assertGVCFIsParseableAndVariantsMatch(outputFile, expected);
    }

    private static VCFHeader getMinimalVCFHeader() {
        final Set<VCFHeaderLine> headerlines = new LinkedHashSet<>();
        VCFStandardHeaderLines.addStandardFormatLines(headerlines, true,
                VCFConstants.GENOTYPE_KEY, VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY, VCFConstants.GENOTYPE_PL_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS);

        VCFStandardHeaderLines.addStandardInfoLines(headerlines, true,
                VCFConstants.DEPTH_KEY,
                VCFConstants.RMS_MAPPING_QUALITY_KEY,
                VCFConstants.MAPPING_QUALITY_ZERO_KEY );

        Arrays.asList(GATKVCFConstants.BASE_QUAL_RANK_SUM_KEY,
               GATKVCFConstants.CLIPPING_RANK_SUM_KEY,
               GATKVCFConstants.MLE_ALLELE_COUNT_KEY,
               GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY,
               GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY,
               GATKVCFConstants.READ_POS_RANK_SUM_KEY)
               .forEach( c -> headerlines.add(GATKVCFHeaderLines.getInfoLine(c)));

        headerlines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY));
        return new VCFHeader(headerlines, Collections.singleton(SAMPLE_NAME));
    }

    private static void assertGVCFIsParseableAndVariantsMatch(File variantFile, List<MinimalData> expected) {
        Assert.assertTrue(variantFile.exists());
        try ( FeatureDataSource<VariantContext> input = new FeatureDataSource<>(variantFile))
        {
            final List<VariantContext> variants = Lists.newArrayList(input.iterator());
            Assert.assertEquals(variants.size(), expected.size());
            assertForEachPair(variants, expected, (vc, md) -> new SimpleInterval(vc).equals(md.location));
            assertForEachPair(variants, expected, (vc, md) -> vc.getGenotype(0).getGQ() == md.gq);
        }
    }

    /**
     * assert that actual and expected have the same length l and
     * assert that assertion(actual[i], expected[i]) is true for every 0 <= i < l
     */
    public static <A, B> void assertForEachPair(Iterable<A> actual, Iterable<B> expected, BiPredicate<A,B> assertion){
        final Iterator<A> iteratorActual = actual.iterator();
        final Iterator<B> iteratorExpected = expected.iterator();
        while( iteratorActual.hasNext() && iteratorExpected.hasNext()){
            final A a = iteratorActual.next();
            final B b = iteratorExpected.next();
            Assert.assertTrue(assertion.test(a, b), "Assertion failed for " + a + " and " + b);
        }
        if( iteratorActual.hasNext() || iteratorExpected.hasNext()){
            Assert.fail("Actual and Expected are different lengths");
        }

    }

    @Test
    public void testAgainstExampleGVCF() throws IOException {
        final List<Number> gqPartitions = Arrays.asList(1, 10, 20, 30, 40, 50, 60);
        final File comparisonFile = new File("src/test/resources/org/broadinstitute/hellbender/utils/variant/writers/small.g.vcf");
        final File outputFile = createTempFile("generated", ".g.vcf");
        final VariantContextBuilder chr = new VariantContextBuilder().chr(CHR1);
        final Allele REF_C = Allele.create("C", true);
        final Allele REF_G = Allele.create("G", true);
        final Allele C = Allele.create("C", false);

        final GenotypeBuilder block1GenotypeBuilder = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF_C, REF_C))
                .DP(35)
                .GQ(57)
                .PL(new int[]{0, 57, 855});
        final VariantContextBuilder block1 = new VariantContextBuilder(null, "1", 14663, 14663, Arrays.asList(REF_C,
                                                                                                              Allele.NON_REF_ALLELE))
                .genotypes(block1GenotypeBuilder.make());

        final VariantContextBuilder block2 = new VariantContextBuilder(null, "1", 14667, 14667, Arrays.asList(REF_G,
                                                                                                              Allele.NON_REF_ALLELE))
                .genotypes(new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF_G,REF_G))
                        .DP(34)
                        .GQ(60)
                        .PL(new int[]{0, 60, 900})
                        .make());

        final VariantContext snp = new VariantContextBuilder(null, "1", 14673, 14673, Arrays.asList(REF_G, C,
                                                                                                    Allele.NON_REF_ALLELE))
                .log10PError(541.77 / -10 )
                .attribute(GATKVCFConstants.BASE_QUAL_RANK_SUM_KEY, "-1.814")
                .attribute(GATKVCFConstants.CLIPPING_RANK_SUM_KEY, "2.599")
                .attribute(VCFConstants.DEPTH_KEY,30)
                .attribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY, new int[]{1,0})
                .attribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, new double[]{0.500, 0.00})
                .attribute(VCFConstants.RMS_MAPPING_QUALITY_KEY, 26.87)
                .attribute(VCFConstants.MAPPING_QUALITY_ZERO_KEY, 0)
                .attribute(GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY, 0.245)
                .attribute(GATKVCFConstants.READ_POS_RANK_SUM_KEY, 0.294)
                .genotypes(new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF_G, C))
                        .AD(new int[]{7,23,0})
                        .DP(30)
                        .GQ(99)
                        .PL(new int[] {570, 0, 212, 591, 281, 872})
                        .attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, new int[]{4,3,23,0})
                        .make())
                .make();

        final GenotypeBuilder block3genotypeBuilder = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF_G, REF_G))
                .DP(29)
                .GQ(54)
                .PL(new int[]{0, 54, 810});

        final VariantContextBuilder block3 = new VariantContextBuilder(null, "1", 14674, 14674, Arrays.asList(REF_G,
                                                                                                              Allele.NON_REF_ALLELE))
                .genotypes(block3genotypeBuilder.make());

        try (VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(outputFile.toPath(), null, false);
             GVCFWriter gvcfWriter = new GVCFWriter(writer, gqPartitions))
        {
            gvcfWriter.writeHeader(getMinimalVCFHeader());

            gvcfWriter.add(block1.genotypes(block1GenotypeBuilder.DP(35).make()).make()); // add with different DP's so that the DP and MIN_DP are correct
            gvcfWriter.add(block1.start(14664).stop(14664).genotypes(block1GenotypeBuilder.DP(35).make()).make());
            gvcfWriter.add(block1.start(14665).stop(14665).genotypes(block1GenotypeBuilder.DP(33).make()).make());
            gvcfWriter.add(block1.start(14666).stop(14666).genotypes(block1GenotypeBuilder.DP(35).make()).make());

            for(int i = 14667; i < 14673; i++){
                gvcfWriter.add(block2.start(i).stop(i).make());
            }

            gvcfWriter.add(snp);

            gvcfWriter.add(block3.start(14674).stop(14674).genotypes(block3genotypeBuilder.DP(28).PL(new int[] { 1, 54, 980}).make()).make()); //vary both DP and PL so that MIN_DP and MIN_PL are both correct
            gvcfWriter.add(block3.start(14675).stop(14675).genotypes(block3genotypeBuilder.DP(29).PL(new int[] { 0, 56, 900}).make()).make());
            gvcfWriter.add(block3.start(14676).stop(14676).genotypes(block3genotypeBuilder.DP(29).PL(new int[] { 0, 59, 810}).make()).make());
        }

        IntegrationTestSpec.assertEqualTextFiles(outputFile, comparisonFile);
    }

    @Test
    public void testOverlappingDeletions() {
        final Allele ref1 = Allele.create("TACACACACATACACACACAC", true);
        final Allele alt1 = Allele.create("T", false);
        final Allele ref2 = Allele.create("TACACACACACTACTA", true);
        final Allele ref3 = Allele.create("T", true);
        final VariantContext deletion1 = new VariantContextBuilder(null, "1", 10000, 10020, Arrays.asList(ref1, alt1,
                Allele.NON_REF_ALLELE))
                .log10PError(1000 / -10 )
                .genotypes(new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(ref1, alt1))
                        .AD(new int[]{7,23,0})
                        .DP(30)
                        .GQ(99)
                        .PL(new int[] {40, 0, 212, 591, 281, 872})
                        .attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, new int[]{4,3,23,0})
                        .make())
                .make();

        final VariantContext deletion2 = new VariantContextBuilder(null, "1", 10010, 10025, Arrays.asList(ref2, alt1,
                Allele.NON_REF_ALLELE))
                .log10PError(10 / -10 )
                .genotypes(new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(ref2, alt1))
                        .AD(new int[]{7,23,0})
                        .DP(30)
                        .GQ(99)
                        .PL(new int[] {40, 0, 212, 591, 281, 872})
                        .attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, new int[]{4,3,23,0})
                        .make())
                .make();

        final VariantContext origRefBlock = makeHomRef("1", 10026, 60, 10050);

        //Let's say that these are "low quality" deletions below the RGQ threshold that get converted to homRefs with all zero PLs
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(ref3, ref3)).DP(30).GQ(0).PL(new int[3]).attribute(VCFConstants.END_KEY, 10025);
        final VariantContextBuilder vcb = new VariantContextBuilder(null, "1", 10010, 10025, Arrays.asList(ref3, Allele.NON_REF_ALLELE))
                .genotypes(gb.make()).attribute(VCFConstants.END_KEY, 10025);
        final VariantContext block2 = (new HomRefBlock(vcb.make(), 20, 100, 2).toVariantContext(SAMPLE_NAME, false));

        final MockVcfWriter mockWriter = new MockVcfWriter();
        final GVCFWriter writer = new GVCFWriter(mockWriter, Arrays.asList(20,100));
        writer.add(deletion1);
        writer.add(block2);
        writer.add(origRefBlock);
        writer.close();
        Assert.assertTrue(mockWriter.emitted.size() == 3);
        Assert.assertTrue(mockWriter.emitted.get(1).getEnd()+1 == mockWriter.emitted.get(2).getStart());
        //The first two blocks overlap, which is fine, but the important thing is that there's no "hole" between the first deletion and the final block

        //note that this is not a reblocking GVCF writer, so starts and ends won't be changed except for ref block merges
        Assert.assertTrue(mockWriter.emitted.get(1).getStart() == block2.getStart());
    }

}
