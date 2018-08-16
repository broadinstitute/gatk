package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter.OrientationBiasUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public final class OxoGReadCountsUnitTest {
    @Test
    public void testUsableRead() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(5, 1, 10000);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 1, 76);

        read.setMappingQuality(60);
        Assert.assertTrue(ReadUtils.readHasReasonableMQ(read));

        read.setMappingQuality(0);
        Assert.assertFalse(ReadUtils.readHasReasonableMQ(read));

        read.setMappingQuality(QualityUtils.MAPPING_QUALITY_UNAVAILABLE);
        Assert.assertFalse(ReadUtils.readHasReasonableMQ(read));
    }

    @Test
    public void testDescriptions() {
        Assert.assertEquals(new OxoGReadCounts().getKeyNames(), Arrays.asList(GATKVCFConstants.F1R2_KEY, GATKVCFConstants.F2R1_KEY));
        Assert.assertEquals(new OxoGReadCounts().getDescriptions(),
                Arrays.asList(
                        GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F1R2_KEY),
                        GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F2R1_KEY))
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
        final int altF1R2 = 1;
        final int altF2R1 = 2;
        final int refF1R2 = 4;
        final int refF2R1 = 8;

        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final Genotype g = new GenotypeBuilder(sample1, alleles).DP(dpDepth).make();
        final Pair<VariantContext, ReadLikelihoods<Allele>> pair = makeReads(altF1R2, altF2R1, refF1R2, refF2R1, refAllele, altAllele, alleles, g);
        final VariantContext vc = pair.getLeft();
        final ReadLikelihoods<Allele> likelihoods = pair.getRight();
        final GenotypeBuilder gb = new GenotypeBuilder(g);
        new OxoGReadCounts().annotate(null, vc, g, gb, likelihoods);

        Assert.assertEquals(OrientationBiasUtils.getF1R2(gb.make()), new int[] {refF1R2, altF1R2});

        Assert.assertEquals(OrientationBiasUtils.getF2R1(gb.make()), new int[] {refF2R1, altF2R1});

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
                ArtificialAnnotationUtils.makeLikelihoods(sample1, refReads, altReads, Arrays.asList(badRead), -100.0, -10.0, -1.1, refAllele, altAllele);

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
