package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class OxoGReadCountsUnitTest {
    @Test
    public void testUsableRead() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(5, 1, 10000);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, 76);

        read.setMappingQuality(60);
        Assert.assertTrue(OxoGReadCounts.isUsableRead(read));

        read.setMappingQuality(0);
        Assert.assertFalse(OxoGReadCounts.isUsableRead(read));

        read.setMappingQuality(QualityUtils.MAPPING_QUALITY_UNAVAILABLE);
        Assert.assertFalse(OxoGReadCounts.isUsableRead(read));
    }

    @Test
    public void testDescriptions() {
        Assert.assertEquals(new OxoGReadCounts().getKeyNames(), Arrays.asList(GATKVCFConstants.OXOG_ALT_F1R2_KEY,
                GATKVCFConstants.OXOG_ALT_F2R1_KEY,
                GATKVCFConstants.OXOG_REF_F1R2_KEY,
                GATKVCFConstants.OXOG_REF_F2R1_KEY,
                GATKVCFConstants.OXOG_FRACTION_KEY), "annots");
        Assert.assertEquals(new OxoGReadCounts().getDescriptions(),
                Arrays.asList(
                        GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_ALT_F1R2_KEY),
                        GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_ALT_F2R1_KEY),
                        GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_REF_F1R2_KEY),
                        GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_REF_F2R1_KEY),
                        GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_FRACTION_KEY))
        );
    }

    private GATKRead makeRead(final Allele ref, final Allele alt, final boolean isRefRead, final boolean isF1R2Read, final PerReadAlleleLikelihoodMap map, int name){
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(10 + "M"), "random_read_" + isRefRead + "_" + isF1R2Read + "_" + name);
        read.setMappingQuality(20);
        if (isF1R2Read){
            F1R2(read);
        } else {
            F2R1(read);
        }
        if (isRefRead){
            ref(map, ref, alt, read);
        } else {
            alt(map, ref, alt, read);
        }
        return read;
    }


    private static final Allele refA = Allele.create("A", true);
    private static final Allele refC = Allele.create("C", true);
    private static final Allele refG = Allele.create("G", true);
    private static final Allele refT = Allele.create("T", true);

    private static final Allele altA = Allele.create("A", false);
    private static final Allele altC = Allele.create("C", false);
    private static final Allele altG = Allele.create("G", false);
    private static final Allele altT = Allele.create("T", false);

    private static final String sample1 = "sample1";
    private static final int dpDepth = 30; //Note: using a different value on purpose so that we can check that reads are preferred over DP


    @DataProvider(name = "allPairs")
    public Object[][] allPairs() {
        final List<Object[]> tests = new ArrayList<>();

        for (final Allele ref : Arrays.asList(refA, refC, refG, refT)){
            for (final Allele alt : Arrays.asList(altA, altC, altG, altT)){
                if (!ref.getBaseString().equals(alt.getBaseString())){
                    tests.add(new Object[]{ref, alt});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "allPairs")
    public void testUsingReads(final Allele refAllele, final Allele altAllele){
        final int alt_F1R2 = 1;
        final int alt_F2R1 = 2;
        final int ref_F1R2 = 4;
        final int ref_F2R1 = 8;

        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final Genotype g = new GenotypeBuilder(sample1, alleles).DP(dpDepth).make();
        final VariantContext vc = makeReads(alt_F1R2, alt_F2R1, ref_F1R2, ref_F2R1, map, refAllele, altAllele, alleles, g);
        final GenotypeBuilder gb = new GenotypeBuilder(g);
        new OxoGReadCounts().annotate(null, vc, g, gb, map);
        final int actual_alt_F1R2 = (int) gb.make().getExtendedAttribute(GATKVCFConstants.OXOG_ALT_F1R2_KEY);
        Assert.assertEquals(actual_alt_F1R2, alt_F1R2, GATKVCFConstants.OXOG_ALT_F1R2_KEY);

        final int actual_alt_F2R1 = (int) gb.make().getExtendedAttribute(GATKVCFConstants.OXOG_ALT_F2R1_KEY);
        Assert.assertEquals(actual_alt_F2R1, alt_F2R1, GATKVCFConstants.OXOG_ALT_F2R1_KEY);

        final int actual_ref_F2R1 = (int) gb.make().getExtendedAttribute(GATKVCFConstants.OXOG_REF_F2R1_KEY);
        Assert.assertEquals(actual_ref_F2R1, ref_F2R1, GATKVCFConstants.OXOG_REF_F2R1_KEY);

        final int actual_ref_F1R2 = (int) gb.make().getExtendedAttribute(GATKVCFConstants.OXOG_REF_F1R2_KEY);
        Assert.assertEquals(actual_ref_F1R2, ref_F1R2, GATKVCFConstants.OXOG_REF_F1R2_KEY);

        final double actual_fraction = (double) gb.make().getExtendedAttribute(GATKVCFConstants.OXOG_FRACTION_KEY);
        final double num = refAllele.equals(refA) || refAllele.equals(refC) ? alt_F2R1 : alt_F1R2;
        final double expectedFraction = num / (alt_F1R2 + alt_F2R1);
        Assert.assertEquals(actual_fraction, expectedFraction, GATKVCFConstants.OXOG_FRACTION_KEY);

        //now test a no-op
        final GenotypeBuilder gb1 = new GenotypeBuilder(g);
        new OxoGReadCounts().annotate(null, vc, null, gb1, map);  //null genotype
        Assert.assertFalse(gb1.make().hasAD());
    }

    private VariantContext makeReads(int alt_F1R2, int alt_F2R1, int ref_F1R2, int ref_F2R1, PerReadAlleleLikelihoodMap map, Allele refAllele, Allele altAllele, List<Allele> alleles, Genotype g) {
        for (int i = 0; i < alt_F1R2; i++) {
            makeRead(refAllele, altAllele, false, true, map, i);
        }
        for (int i = 0; i < alt_F2R1; i++) {
            makeRead(refAllele, altAllele, false, false, map, i);
        }
        for (int i = 0; i < ref_F1R2; i++) {
            makeRead(refAllele, altAllele, true, true, map, i);
        }
        for (int i = 0; i < ref_F2R1; i++) {
            makeRead(refAllele, altAllele, true, false, map, i);
        }
        //throw in one non-informative read
        final GATKRead badRead = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(10 + "M"));
        badRead.setMappingQuality(20);
        map.add(badRead, refAllele, -1.0);
        map.add(badRead, altAllele, -1.1); //maybe it's ref, maybe it's alt, too close to call -> not informative

        return new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(Arrays.asList(g)).make();
    }

    private void ref(PerReadAlleleLikelihoodMap map, Allele a, Allele c, GATKRead read) {
        map.add(read, a, -1.0);
        map.add(read, c, -100.0);  //try to fool it - add another likelihood to same read
    }

    private void alt(PerReadAlleleLikelihoodMap map, Allele a, Allele c, GATKRead read) {
        map.add(read, a, -10.0);
        map.add(read, c, -1.0);      //try to fool it - add another likelihood to same read
    }

    private void F2R1(GATKRead read) {
        read.setIsReverseStrand(false);
        read.setIsPaired(true);
        read.setIsSecondOfPair();
    }

    private void F1R2(GATKRead read) {
        read.setIsPaired(true);
        read.setIsReverseStrand(false);
        read.setIsFirstOfPair();
    }
}
