package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class GVCFWriterUnitTest extends BaseTest {
    private static class MockWriter implements VariantContextWriter {
        final List<VariantContext> emitted = new ArrayList<>();
        boolean headerWritten = false;
        boolean closed = false;
        boolean error = false;

        @Override
        public void writeHeader(VCFHeader header) {
            headerWritten = true;
        }

        @Override
        public void close() {
            closed = true;
        }

        @Override
        public boolean checkError() {
            return error;
        }

        @Override
        public void add(VariantContext vc) {
            emitted.add(vc);
        }
    }

    private MockWriter mockWriter;
    private final List<Integer> standardPartition = Arrays.asList(1, 10, 20);
    private final Allele REF = Allele.create("N", true);
    private final Allele ALT = Allele.create("A");
    private final List<Allele> ALLELES = Arrays.asList(REF, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
    private final String SAMPLE_NAME = "XXYYZZ";

    @BeforeMethod
    public void setUp() throws Exception {
        mockWriter = new MockWriter();
    }

    @Test
    public void testHeaderWriting() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);
        writer.writeHeader(new VCFHeader());
        Assert.assertTrue(mockWriter.headerWritten);
    }

    @Test
    public void testClose() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);
        writer.close();
        Assert.assertTrue(mockWriter.closed);
    }

    private VariantContext makeHomRef(final String contig, final int start, final int GQ) {
        final VariantContextBuilder vcb = new VariantContextBuilder("test", contig, start, start, ALLELES);
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, REF));
        gb.GQ(GQ);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    private VariantContext makeHomRefAlt(final String contig, final int start, final int GQ) {
        final VariantContextBuilder vcb = new VariantContextBuilder("test", contig, start, start, Arrays.asList(REF, ALT));
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, REF));
        gb.GQ(GQ);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    private VariantContext makeNonRef(final String contig, final int start, final int GQ) {
        final VariantContextBuilder vcb = new VariantContextBuilder("test", contig, start, start, Arrays.asList(REF, ALT));
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, ALT));
        gb.GQ(GQ);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    private VariantContext makeDeletion(final String contig, final int start, final int size) {
        final String del = Utils.dupChar('A', size);
        final String alt = del.substring(0, 1);
        final VariantContext vc = GATKVariantContextUtils.makeFromAlleles("test", contig, start, Arrays.asList(del, alt));
        final VariantContextBuilder vcb = new VariantContextBuilder(vc);
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(vc.getReference(), vc.getAlternateAllele(0)));
        gb.GQ(50);
        gb.DP(10);
        gb.AD(new int[]{1, 2});
        gb.PL(new int[]{0, 10, 100});
        return vcb.genotypes(gb.make()).make();
    }

    @Test
    public void testCloseEmitsLastVariant() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);

        writer.close();
        Assert.assertTrue(mockWriter.closed);
        Assert.assertEquals(mockWriter.emitted.size(), 1);
    }

    @Test
    public void testCloseDoesntEmitsLastVariantWhenNonRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeNonRef("20", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);

        writer.close();
        Assert.assertTrue(mockWriter.closed);
        Assert.assertEquals(mockWriter.emitted.size(), 1);
    }

    @Test
    public void testCrossingContigBoundaryRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        writer.add(makeHomRef("20", 2, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef("21", 3, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);

        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(1), "21", 3, 3, false);
    }

    @Test
    public void testCrossingContigBoundaryToLowerPositionsRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 30, 30));
        writer.add(makeHomRef("20", 31, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef("21", 10, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), "20", 30, 31, false);
        writer.add(makeNonRef("21", 11, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(1), "21", 10, 10, false);
        assertGoodVC(mockWriter.emitted.get(2), "21", 11, 11, true);
    }

    @Test
    public void testCrossingContigBoundaryFromNonRefToLowerPositionsRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeNonRef("20", 20, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        writer.add(makeHomRef("21", 10, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        assertGoodVC(mockWriter.emitted.get(0), "20", 20, 20, true);
        writer.add(makeNonRef("21", 11, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(1), "21", 10, 10, false);
        assertGoodVC(mockWriter.emitted.get(2), "21", 11, 11, true);
    }

    @Test
    public void testCrossingContigBoundaryNonRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        writer.add(makeHomRef("20", 2, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeNonRef("21", 3, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "21", 3, 3, true);
    }

    @Test
    public void testCrossingContigBoundaryNonRefThenNonRef() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeNonRef("20", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 1);
        writer.add(makeNonRef("21", 1, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 1, true);
        assertGoodVC(mockWriter.emitted.get(1), "21", 1, 1, true);
    }

    private void assertGoodVC(final VariantContext vc, final String contig, final int start, final int stop, final boolean nonRef) {
        Assert.assertEquals(vc.getContig(), contig);
        Assert.assertEquals(vc.getStart(), start);
        Assert.assertEquals(vc.getEnd(), stop);
        if ( nonRef ) {
            Assert.assertNotEquals(vc.getAlternateAllele(0), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        } else {
            Assert.assertEquals(vc.getNAlleles(), 2);
            Assert.assertEquals(vc.getAlternateAllele(0), GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
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
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 30));
        writer.add(makeHomRef("20", 2, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeNonRef("20", 3, 30));
        writer.add(makeHomRef("20", 4, 30));
        writer.add(makeHomRef("20", 5, 30));
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 3, 3, true);
        writer.close();
        assertGoodVC(mockWriter.emitted.get(2), "20", 4, 5, false);
    }

    @Test
    public void testEmittingTwoBands() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        Assert.assertEquals(mockWriter.emitted.size(), 0);
        writer.add(makeHomRef("20", 3, 50));
        writer.add(makeHomRef("20", 4, 50));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 3, 4, false);
    }

    @Test
    public void testNonContiguousBlocks() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        writer.add(makeHomRef("20", 10, 0));
        writer.add(makeHomRef("20", 11, 0));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 2);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 10, 11, false);
    }

    @Test
    public void testDeletion() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        writer.add(makeDeletion("20", 3, 3));
        writer.add(makeHomRef("20", 4, 0));
        writer.add(makeHomRef("20", 5, 0));
        writer.add(makeHomRef("20", 6, 0));
        writer.add(makeHomRef("20", 7, 0));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        assertGoodVC(mockWriter.emitted.get(1), "20", 3, 5, true);
        assertGoodVC(mockWriter.emitted.get(2), "20", 6, 7, false);
    }

    @Test
    public void testHomRefAlt() {
        final GVCFWriter writer = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);

        writer.add(makeHomRef("20", 1, 0));
        writer.add(makeHomRef("20", 2, 0));
        writer.add(makeHomRefAlt("20", 3, 0));
        writer.add(makeHomRef("20", 4, 0));
        writer.add(makeHomRef("20", 5, 0));
        writer.add(makeHomRef("20", 6, 0));
        writer.add(makeHomRef("20", 7, 0));
        writer.close();
        Assert.assertEquals(mockWriter.emitted.size(), 3);
        assertGoodVC(mockWriter.emitted.get(0), "20", 1, 2, false);
        Assert.assertFalse(mockWriter.emitted.get(1).hasAttribute("END"));
        Assert.assertFalse(mockWriter.emitted.get(1).hasAttribute("BLOCK_SIZE"));
        assertGoodVC(mockWriter.emitted.get(2), "20", 4, 7, false);
    }

    @DataProvider(name = "GoodBandPartitionData")
    public Object[][] makeBandPartitionData() {
        return new Object[][]{
                {Arrays.asList(1), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,Integer.MAX_VALUE))},
                {Arrays.asList(1, 2, 3), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,2), Range.closedOpen(2,3),Range.closedOpen(3,Integer.MAX_VALUE))},
                {Arrays.asList(1, 10), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,10), Range.closedOpen(10, Integer.MAX_VALUE))},
                {Arrays.asList(1, 10, 30), Arrays.asList(Range.closedOpen(0,1), Range.closedOpen(1,10), Range.closedOpen(10, 30),Range.closedOpen(30,Integer.MAX_VALUE))},
        };
    }

    @Test(dataProvider = "GoodBandPartitionData")
    public void testGoodPartitions(final List<Integer> partitions, List<Range<Integer>> expected) {
        final RangeMap<Integer, Range<Integer>> ranges = GVCFWriter.parsePartitions(partitions);
        Assert.assertEquals(new ArrayList<>(ranges.asMapOfRanges().keySet()), expected);

    }

    @DataProvider(name = "BadBandPartitionData")
    public Object[][] makeBadBandPartitionData() {
        return new Object[][]{
                {null},
                {Collections.emptyList()},
                {Arrays.asList(10, 1, 30)},
                {Arrays.asList(-1, 1)},
                {Arrays.asList(1, null, 10)}
        };
    }

    @Test(dataProvider = "BadBandPartitionData", expectedExceptions = IllegalArgumentException.class)
    public void testBadPartitionsThrowException(final List<Integer> partitions){
        GVCFWriter.parsePartitions(partitions); // we should explode here
    }

    @Test
    public void testCheckError(){
        GVCFWriter gvcfWriter = new GVCFWriter(mockWriter, standardPartition, HomoSapiensConstants.DEFAULT_PLOIDY);
        mockWriter.error = false;
        Assert.assertEquals(gvcfWriter.checkError(), mockWriter.checkError());
        mockWriter.error = true;
        Assert.assertEquals(gvcfWriter.checkError(), mockWriter.checkError());
    }

    @Test
    public void testToVCFHeaderLine() {
        final Range<Integer> band = Range.closedOpen(10,20);
        Assert.assertEquals(GVCFWriter.rangeToVCFHeaderLine(band).getKey(), "GVCFBlock10-20", "Wrong key for " + band);
        Assert.assertEquals(GVCFWriter.rangeToVCFHeaderLine(band).getValue(), "minGQ=10(inclusive),maxGQ=20(exclusive)", "Wrong value for" + band);

    }

}
