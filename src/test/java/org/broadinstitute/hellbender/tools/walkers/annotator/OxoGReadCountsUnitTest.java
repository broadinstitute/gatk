package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
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
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

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

        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final Genotype g = new GenotypeBuilder(sample1, alleles).DP(dpDepth).make();
        final Pair<VariantContext, ReadLikelihoods<Allele>> pair = makeReads(alt_F1R2, alt_F2R1, ref_F1R2, ref_F2R1, refAllele, altAllele, alleles, g);
        final VariantContext vc = pair.getLeft();
        final ReadLikelihoods<Allele> likelihoods = pair.getRight();
        final GenotypeBuilder gb = new GenotypeBuilder(g);
        new OxoGReadCounts().annotate(null, vc, g, gb, likelihoods);
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
        new OxoGReadCounts().annotate(null, vc, null, gb1, likelihoods);  //null genotype
        Assert.assertFalse(gb1.make().hasAD());
    }

    private Pair<VariantContext, ReadLikelihoods<Allele>> makeReads(int alt_F1R2, int alt_F2R1, int ref_F1R2, int ref_F2R1, Allele refAllele, Allele altAllele, List<Allele> alleles, Genotype g) {
        final List<GATKRead> altReads = Stream.concat(IntStream.range(0, alt_F1R2).mapToObj(i -> makeRead(false, true, i)),
                IntStream.range(0, alt_F2R1).mapToObj(i -> makeRead(false, false, i))).collect(Collectors.toList());
        final List<GATKRead> refReads = Stream.concat(IntStream.range(0, ref_F1R2).mapToObj(i -> makeRead(true, true, i)),
                IntStream.range(0, ref_F2R1).mapToObj(i -> makeRead(true, false, i))).collect(Collectors.toList());
        final GATKRead badRead = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(10 + "M"));
        badRead.setMappingQuality(20);

        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.makeLikelihoods(sample1, refReads, altReads, Arrays.asList(badRead), -100.0, -10.0, -1.1, refAllele, altAllele);

        return ImmutablePair.of(new VariantContextBuilder("test", "20", 10, 10, alleles).genotypes(Arrays.asList(g)).make(),
                likelihoods);
    }

    private GATKRead makeRead(final boolean isRefRead, final boolean isF1R2Read, int name){
        final GATKRead read = ArtificialReadUtils.createArtificialRead(TextCigarCodec.decode(10 + "M"), "random_read_" + isRefRead + "_" + isF1R2Read + "_" + name);
        read.setMappingQuality(20);
        if (isF1R2Read){
            read.setIsPaired(true);
            read.setIsReverseStrand(false);
            read.setIsFirstOfPair();
        } else {
            read.setIsReverseStrand(false);
            read.setIsPaired(true);
            read.setIsSecondOfPair();
        }
        return read;
    }
}
