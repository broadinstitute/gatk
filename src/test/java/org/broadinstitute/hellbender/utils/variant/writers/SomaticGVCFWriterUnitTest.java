package org.broadinstitute.hellbender.utils.variant.writers;

import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.testutils.VariantContextTestUtils.makeSomaticRef;

public class SomaticGVCFWriterUnitTest {
    private static final List<Number> standardPartition = ImmutableList.of(-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5);
    private static final List<Number> precisionTwoPartition = ImmutableList.of(-.55, -.5, -.45, -.4);
    private static final List<Number> precisionThreePartition = ImmutableList.of(-.55, -.501, -.5, -.45, -.4);
    private static final Allele REF = Allele.create("G", true);
    private static final List<Allele> ALLELES = ImmutableList.of(REF, Allele.NON_REF_ALLELE);
    private static final String SAMPLE_NAME = "XXYYZZ";


    @Test
    public void testValueBinning() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        SomaticGVCFWriter writer = new SomaticGVCFWriter(mockWriter, standardPartition);
        //derives partitionPrecision 1 from standardPartition values
        Assert.assertTrue(writer.convertLODtoInt(2.3) == 23);
        Assert.assertTrue(writer.convertLODtoInt(-2.3) == -23);
        Assert.assertTrue(writer.convertLODtoInt(2.0) == 20);

        writer = new SomaticGVCFWriter(mockWriter, standardPartition, 2);
        Assert.assertTrue(writer.convertLODtoInt(2.3) == 230);
        Assert.assertTrue(writer.convertLODtoInt(-2.3) == -230);
        Assert.assertTrue(writer.convertLODtoInt(2.0) == 200);

        writer = new SomaticGVCFWriter(mockWriter, standardPartition, 3);
        Assert.assertTrue(writer.convertLODtoInt(2.3) == 2300);
        Assert.assertTrue(writer.convertLODtoInt(2.33) == 2330);
        Assert.assertTrue(writer.convertLODtoInt(2.333) == 2333);
        Assert.assertTrue(writer.convertLODtoInt(2.3333) == 2333);
        Assert.assertTrue(writer.convertLODtoInt(2.3337) == 2334);
        Assert.assertTrue(writer.convertLODtoInt(-2.3) == -2300);
        Assert.assertTrue(writer.convertLODtoInt(2.0) == 2000);
    }

    @Test
    public void testAddingAndMerging() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        final SomaticGVCFWriter writer = new SomaticGVCFWriter(mockWriter, standardPartition);
        final GenotypeBuilder gb = new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(REF, REF));
        int pos = 1;
        final VariantContextBuilder vcb = new VariantContextBuilder("source", "contig", pos, pos, ALLELES);

        Genotype g = gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, -3.4).make();
        Assert.assertFalse(writer.gvcfBlockCombiner.genotypeCanBeMergedInCurrentBlock(g));  //should be false if there's no current block
        writer.add(vcb.genotypes(gb.make()).make());

        vcb.start(++pos).stop(pos);
        g = gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, -3.2).make();
        Assert.assertTrue(writer.gvcfBlockCombiner.genotypeCanBeMergedInCurrentBlock(g));
        writer.add(vcb.genotypes(gb.make()).make());

        //test that inclusive lower bounds function properly
        vcb.start(++pos).stop(pos);
        g = gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, -3.0).make();
        Assert.assertFalse(writer.gvcfBlockCombiner.genotypeCanBeMergedInCurrentBlock(g));
        writer.add(vcb.genotypes(gb.make()).make());

        vcb.start(++pos).stop(pos);
        g = gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, -2.7).make();
        Assert.assertTrue(writer.gvcfBlockCombiner.genotypeCanBeMergedInCurrentBlock(g));
        writer.add(vcb.genotypes(gb.make()).make());

        vcb.start(++pos).stop(pos);
        g = gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 600.0).make();
        Assert.assertFalse(writer.gvcfBlockCombiner.genotypeCanBeMergedInCurrentBlock(g));
        writer.add(vcb.genotypes(gb.make()).make());

        vcb.start(++pos).stop(pos);
        g = gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, 601.1).make();
        Assert.assertTrue(writer.gvcfBlockCombiner.genotypeCanBeMergedInCurrentBlock(g));
        writer.add(vcb.genotypes(gb.make()).make());

        writer.close();
        Assert.assertTrue(mockWriter.closed);
        Assert.assertEquals(mockWriter.emitted.size(), 3);
    }

    @Test
    public void testPrecision() {
        final MockVcfWriter mockWriter = new MockVcfWriter();
        SomaticGVCFWriter writer = new SomaticGVCFWriter(mockWriter, precisionTwoPartition);
        Assert.assertTrue(((SomaticGVCFBlockCombiner)writer.gvcfBlockCombiner).partitionPrecision == 2);

        writer = new SomaticGVCFWriter(mockWriter, precisionThreePartition);
        Assert.assertTrue(((SomaticGVCFBlockCombiner)writer.gvcfBlockCombiner).partitionPrecision == 3);
        writer.add(makeSomaticRef("chr1", 1, -0.500005, 10));
        writer.close();
        VariantContext vc = mockWriter.emitted.get(0);
        //partitionPrecision does not affect the precision of the minLOD for the block
        Assert.assertEquals(vc.getGenotype(0).getExtendedAttribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY), -0.500005);
    }
}