package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public final class BaseQualityRankSumTestUnitTest {

    private static final String SAMPLE_1 = "NA1";
    private static final String SAMPLE_2 = "NA2";

    private static final Allele REF = Allele.create("T", true);
    private static final Allele ALT = Allele.create("A", false);

    public static GATKRead makeRead(final int qual) {
        return ArtificialAnnotationUtils.makeRead(qual, 50);
    }

    public static VariantContext makeVC(final Allele refAllele, final Allele altAllele) {
        final double[] genotypeLikelihoods1 = {30, 0, 190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // SAMPLE_1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(30).make());
        // SAMPLE_2 -> A/T with GQ 40
        testGC.add(new GenotypeBuilder(SAMPLE_2).alleles(Arrays.asList(refAllele, altAllele)).PL(genotypeLikelihoods1).GQ(40).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    @Test
    public void testBaseQual() {
        final InfoFieldAnnotation ann = new BaseQualityRankSumTest();
        final String key = GATKVCFConstants.BASE_QUAL_RANK_SUM_KEY;
        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        final int[] altBaseQuals = {10, 20};
        final int[] refBaseQuals = {50, 60};
        final List<GATKRead> refReads = Arrays.stream(refBaseQuals).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> altReads = Arrays.stream(altBaseQuals).mapToObj(i -> makeRead(i)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE_1, refReads, altReads, -100.0, -100.0, REF, ALT);

        final ReferenceContext ref = null;
        final VariantContext vc = makeVC(REF, ALT);

        final Map<String, Object> annotate = ann.annotate(ref, vc, likelihoods);

        final double zScore = mannWhitneyU.test(new double[]{altBaseQuals[0], altBaseQuals[1]}, new double[]{refBaseQuals[0], refBaseQuals[1]}, MannWhitneyU.TestType.FIRST_DOMINATES).getZ();
        final String zScoreStr = String.format("%.3f", zScore);
        Assert.assertEquals(annotate.get(key), zScoreStr);

        Assert.assertEquals(ann.getDescriptions().size(), 1);
        Assert.assertEquals(ann.getDescriptions().get(0).getID(), key);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key);
    }

    @Test
    public void testEmptyIfNoGenotypes() throws Exception {
        final BaseQualityRankSumTest ann = new BaseQualityRankSumTest();

        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of(SAMPLE_1, reads);
        final SampleList sampleList = new IndexedSampleList(Arrays.asList(SAMPLE_1));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.NO_CALL));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        final Map<String, Object> annotate = ann.annotate(null, when(mock(VariantContext.class).getGenotypesOrderedByName()).thenReturn(Collections.<Genotype>emptyList()).getMock(), likelihoods);
        Assert.assertTrue(annotate.isEmpty());
    }

    @Test
    public void testNullLikelihoodsReturnsEmpty() {
        final BaseQualityRankSumTest ann = new BaseQualityRankSumTest();
        Assert.assertEquals(ann.annotate(null, mock(VariantContext.class), null), Collections.emptyMap());
    }
}