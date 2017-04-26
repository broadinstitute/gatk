package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_MappingQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RankSumTest;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public final class MappingQualityRankSumTestUnitTest {
    private  static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    private static final String sample1 = "NA1";
    private static final String sample2 = "NA2";

    private VariantContext makeVC( final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30,0,190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // sample1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(sample1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // sample2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(sample2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    private GATKRead makeRead(final int mq) {
        Cigar cigar = TextCigarCodec.decode("10M");
        final GATKRead read = ArtificialReadUtils.createArtificialRead(cigar);
        read.setMappingQuality(mq);
        return read;
    }

    @Test
    public void testMQ(){
        final InfoFieldAnnotation ann = new MappingQualityRankSumTest();
        final String key = GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY;
        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        final int[] altMappingQualities = {10, 20};
        final int[] refMappingQualities = {100, 110};
        final List<GATKRead> refReads = Arrays.stream(refMappingQualities).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> altReads = Arrays.stream(altMappingQualities).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);
        final VariantContext vc = makeVC(REF, ALT);

        final Map<String, Object> annotate = ann.annotate(null, vc, likelihoods);

        final double zScore = mannWhitneyU.test(new double[]{altMappingQualities[0], altMappingQualities[1]}, new double[]{refMappingQualities[0], refMappingQualities[1]}, MannWhitneyU.TestType.FIRST_DOMINATES).getZ();
        final String zScoreStr = String.format("%.3f", zScore);
        Assert.assertEquals(annotate.get(key), zScoreStr);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key);
    }

    @Test
    public void testAS_MQRaw(){
        final AS_RankSumTest ann = new AS_MappingQualityRankSumTest();
        final String key1 = GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY;
        final String key2 = GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY;

        final int[] altMappingQualities = {10, 20};
        final int[] refMappingQualities = {100, 110};
        final List<GATKRead> refReads = Arrays.stream(refMappingQualities).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> altReads = Arrays.stream(altMappingQualities).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                AnnotationArtificialData.makeLikelihoods(sample1, refReads, altReads, -100.0, -100.0, REF, ALT);

        final ReferenceContext ref= null;
        final VariantContext vc= makeVC(REF, ALT);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);

        final String expected = refMappingQualities[0] + ",1," + refMappingQualities[1] + ",1" + AS_RankSumTest.PRINT_DELIM + altMappingQualities[0] + ",1," + altMappingQualities[1] + ",1";
        Assert.assertEquals(annotate.get(key1), expected);
        Assert.assertEquals(annotateRaw.get(key1), expected);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key1);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key2);
    }
}
