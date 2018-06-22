package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQualityRankSumTestUnitTest;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by emeryj on 8/11/17.
 */
public class AS_BaseQualityRankSumTestUnitTest extends ReducibleAnnotationBaseTest {

    private static final String SAMPLE_1 = "NA1";
    private static final String SAMPLE_2 = "NA2";

    private static final Allele REF = Allele.create("T", true);
    private static final Allele ALT1 = Allele.create("A", false);
    private static final Allele ALT2 = Allele.create("AA", false);


    private static VariantContext makeVC(final Allele refAllele, final Allele alt1Allele, final Allele alt2Allele) {
        final double[] genotypeLikelihoods1 = {30, 0, 190};
        final GenotypesContext testGC = GenotypesContext.create(2);
        // SAMPLE_1 -> A/T with GQ 30
        testGC.add(new GenotypeBuilder(SAMPLE_1).alleles(Arrays.asList(refAllele, alt1Allele, alt2Allele)).PL(genotypeLikelihoods1).GQ(30).make());

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, alt1Allele, alt2Allele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }


    @Override
    protected List<Annotation> getAnnotationsToUse() {
        return Collections.singletonList(new AS_BaseQualityRankSumTest());
    }

    @Override
    protected String getRawKey() {
        return GATKVCFConstants.AS_RAW_BASE_QUAL_RANK_SUM_KEY;
    }

    @Override
    protected String getKey() {
        return GATKVCFConstants.AS_BASE_QUAL_RANK_SUM_KEY;
    }

    @Test
    public void testBaseQualRawAnnotate() {
        final AS_RankSumTest ann =  new AS_BaseQualityRankSumTest();
        final String key1 = GATKVCFConstants.AS_RAW_BASE_QUAL_RANK_SUM_KEY;
        final String key2 = GATKVCFConstants.AS_BASE_QUAL_RANK_SUM_KEY;

        final int[] alt1BaseQuals = {10, 20, 30};
        final int[] alt2BaseQuals = {30, 40, 61};
        final int[] refBaseQuals = {50, 60};
        final List<GATKRead> refReads = Arrays.stream(refBaseQuals).mapToObj(i -> BaseQualityRankSumTestUnitTest.makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> alt1Reads = Arrays.stream(alt1BaseQuals).mapToObj(i -> BaseQualityRankSumTestUnitTest.makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> alt2Reads = Arrays.stream(alt2BaseQuals).mapToObj(i -> BaseQualityRankSumTestUnitTest.makeRead(i)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeTriAllelicLikelihoods(SAMPLE_1, refReads, alt1Reads, alt2Reads, new ArrayList<GATKRead>(), -100.0, -100.0, -100,0, REF, ALT1, ALT2);

        final ReferenceContext ref = null;
        final VariantContext vc = makeVC(REF, ALT1, ALT2);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotateNonRaw = ann.annotate(ref, vc, likelihoods);

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();

        final MannWhitneyU.Result expectedAlt1 = mannWhitneyU.test(Arrays.stream(alt1BaseQuals).asDoubleStream().toArray(), Arrays.stream(refBaseQuals).asDoubleStream().toArray(), MannWhitneyU.TestType.FIRST_DOMINATES);
        final MannWhitneyU.Result expectedAlt2 = mannWhitneyU.test(Arrays.stream(alt2BaseQuals).asDoubleStream().toArray(), Arrays.stream(refBaseQuals).asDoubleStream().toArray(), MannWhitneyU.TestType.FIRST_DOMINATES);

        // Binning the results to reflect the Histogram behavior
        String firstExpected = String.format("%.1f",Math.round(Math.floor((expectedAlt1.getZ() )/0.1))*0.1);
        String secondExpected = String.format("%.1f",Math.round(Math.floor((expectedAlt2.getZ() )/0.1))*0.1);

        // Note, when we output the raw annotated RankSum score, we output the MannWhitneyU test Z value as a histogram for each alt allele
        final String expectedAnnotation = AS_RankSumTest.PRINT_DELIM + firstExpected + ",1" + AS_RankSumTest.PRINT_DELIM + secondExpected + ",1";

        final MannWhitneyU.Result annotateResult = mannWhitneyU.test(Stream.concat(Arrays.stream(alt1BaseQuals).boxed(), Arrays.stream(alt2BaseQuals).boxed()).mapToDouble(i->(double)i).toArray(), Arrays.stream(refBaseQuals).asDoubleStream().toArray(), MannWhitneyU.TestType.FIRST_DOMINATES);
        final double annotateZScore = annotateResult.getZ();

        Assert.assertEquals(annotateRaw.get(key1),    expectedAnnotation);
        Assert.assertEquals(annotateNonRaw.get(key2), String.format("%.3f",annotateZScore));

        Assert.assertEquals(ann.getRawDescriptions().size(), 1);
        Assert.assertEquals(ann.getRawDescriptions().get(0).getID(), key1);
        Assert.assertEquals(ann.getKeyNames().size(), 1);
        Assert.assertEquals(ann.getKeyNames().get(0), key2);
    }

    // This test is expected to fail because the underlying VC upon which the rankSumTests are calculated has multiple samples, which we would like to throw an error in concordance with GATK3
    @Test (expectedExceptions = IllegalStateException.class)
    public void testMultipleSamplesL() {
        final AS_RankSumTest ann =  new AS_BaseQualityRankSumTest();

        final int[] alt1BaseQuals = {10, 20, 30};
        final int[] refBaseQuals = {50, 60};
        final List<GATKRead> refReads = Arrays.stream(refBaseQuals).mapToObj(i -> BaseQualityRankSumTestUnitTest.makeRead(i)).collect(Collectors.toList());
        final List<GATKRead> alt1Reads = Arrays.stream(alt1BaseQuals).mapToObj(i -> BaseQualityRankSumTestUnitTest.makeRead(i)).collect(Collectors.toList());
        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE_1, refReads, alt1Reads,  -100.0, -100.0, REF, ALT1);

        final ReferenceContext ref = null;
        final VariantContext vc = BaseQualityRankSumTestUnitTest.makeVC(REF, ALT1);

        final Map<String, Object> annotateRaw = ann.annotateRawData(ref, vc, likelihoods);
        final Map<String, Object> annotateNonRaw = ann.annotate(ref, vc, likelihoods);
    }

}