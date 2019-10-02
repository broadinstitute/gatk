package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class RMSMappingQualityUnitTest {

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull() throws Exception {
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new RMSMappingQuality();
        cov.annotate(referenceContext, vc, null); //vc can't be null
    }

    @Test
    public void testDescriptions() throws Exception {
        final RMSMappingQuality cov = new RMSMappingQuality();
        Assert.assertEquals(cov.getDescriptions().size(), 1);
        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.RMS_MAPPING_QUALITY_KEY);
        Assert.assertEquals(cov.getRawDescriptions().size(), 1);
        Assert.assertEquals(cov.getRawDescriptions().get(0).getID(), GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY);
        RMSMappingQuality annotationClass = new RMSMappingQuality();
        Assert.assertEquals(annotationClass.getPrimaryRawKey(), GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY);
        Assert.assertEquals(annotationClass.getKeyNames(), Sets.newHashSet(VCFConstants.RMS_MAPPING_QUALITY_KEY, GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
    }

    @Test
    public void testNullLikelihoods() throws Exception {
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final RMSMappingQuality cov = new RMSMappingQuality();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);
        Assert.assertTrue(annotate.isEmpty());

        Assert.assertEquals(cov.getDescriptions().get(0).getID(), VCFConstants.RMS_MAPPING_QUALITY_KEY);
        Assert.assertEquals(cov.getRawDescriptions().get(0).getID(), GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY);
    }

    /**
     * @return get a new variant context with only alleles and position filled in
     */
    private static VariantContext makeVC() {
        final GenotypesContext testGC = GenotypesContext.create(2);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr("1").start(15L).stop(15L).genotypes(testGC).make();
    }

    @Test
    public void testLikelihoods(){

        final Allele REF = Allele.create("A", true);
        final Allele ALT = Allele.create("C", true);

        final int[] MQs = {1,2,3,4,5,6,7,8,9,10, QualityUtils.MAPPING_QUALITY_UNAVAILABLE};
        final List<Integer> MQsList = Arrays.asList(ArrayUtils.toObject(MQs));

        //MQ 255 are excluded from the calculations, we test it here.
        final List<Integer> MQsListOK = new ArrayList<>(MQsList);
        //NOTE: if we just call remove(i), Java thinks i is an index.
        //A workaround for this overloading bogosity to to call removeAll and pass a collection
        //(casting i to (Object) would work too but it's more error prone)
        MQsListOK.removeAll(Collections.singleton(QualityUtils.MAPPING_QUALITY_UNAVAILABLE));


        final List<GATKRead> reads = Arrays.stream(MQs).mapToObj(mq -> ArtificialAnnotationUtils.makeRead(20, mq)).collect(Collectors.toList());
        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods("sample1", reads, -1.0, REF, ALT);


        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new RMSMappingQuality().annotate(referenceContext, vc, likelihoods);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(VCFConstants.RMS_MAPPING_QUALITY_KEY), "annots");
        final double rms= MathUtils.sumOfSquares(MQsListOK); //only those are MQ0
        Assert.assertNull(annotate.get(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        Assert.assertEquals(annotate.get(VCFConstants.RMS_MAPPING_QUALITY_KEY), String.format("%.2f", Math.sqrt(rms/(reads.size()-1))));
    }

    @Test
    public void testEmptyLikelihoods() throws Exception {
        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.create("A")));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(sampleList, alleleList, readsBySample);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new RMSMappingQuality().annotate(referenceContext, vc, likelihoods);
        Assert.assertTrue(annotate.isEmpty());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAllNull_AS() throws Exception {
        final VariantContext vc= null;
        final ReferenceContext referenceContext= null;
        final InfoFieldAnnotation cov = new AS_RMSMappingQuality();
        cov.annotate(referenceContext, vc, null); //vc can't be null
    }

    @Test
    public void testDescriptions_AS() throws Exception {
        final ReducibleAnnotation cov = new AS_RMSMappingQuality();
        Assert.assertEquals(cov.getRawDescriptions().size(), 1);
        Assert.assertEquals(cov.getRawDescriptions().get(0).getID(), GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY);
    }

    @Test
    public void testNullLikelihoods_AS() throws Exception {
        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final AS_RMSMappingQuality cov = new AS_RMSMappingQuality();
        final Map<String, Object> annotate = cov.annotate(referenceContext, vc, null);
        Assert.assertTrue(annotate.isEmpty());

        Assert.assertEquals(cov.getRawDescriptions().size(), 1);
        Assert.assertEquals(cov.getRawDescriptions().get(0).getID(), GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY);
    }

    @Test
    public void testLikelihoods_AS(){

        final Allele REF = Allele.create("A", true);
        final Allele ALT = Allele.create("C");

        final int[] MQs = {1,2,3,4,5,6,7,8,9,10, QualityUtils.MAPPING_QUALITY_UNAVAILABLE};
        final List<Integer> MQsList = Arrays.asList(ArrayUtils.toObject(MQs));

        //MQ 255 are excluded from the calculations, we test it here.
        final List<Integer> MQsListOK = new ArrayList<>(MQsList);
        //NOTE: if we just call remove(i), Java thinks i is an index.
        //A workaround for this overloading bogosity to to call removeAll and pass a collection
        //(casting i to (Object) would work too but it's more error prone)
        MQsListOK.removeAll(Collections.singleton(QualityUtils.MAPPING_QUALITY_UNAVAILABLE));

        final List<GATKRead> reads = Arrays.stream(MQs).mapToObj(mq -> ArtificialAnnotationUtils.makeRead(30, mq)).collect(Collectors.toList());

        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods("sample1", reads, -10.0, REF, ALT);

        final VariantContext vc = makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new AS_RMSMappingQuality().annotateRawData(referenceContext, vc, likelihoods);
        Assert.assertEquals(annotate.size(), 1, "size");
        Assert.assertEquals(annotate.keySet(), Collections.singleton(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY), "annots");
        final double rms= MathUtils.sumOfSquares(MQsListOK); //only those are MQ0
        final String[] split =((String)annotate.get(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY)).split(AnnotationUtils.ALLELE_SPECIFIC_SPLIT_REGEX);
        Assert.assertEquals(split.length, 2);
        Assert.assertEquals(split[0], String.format("%.2f", rms));
        Assert.assertEquals(split[1], String.format("%.2f", 0.0));
    }

    @Test
    public void testLikelihoodsEmpty_AS() throws Exception {
        final List<GATKRead> reads = Collections.emptyList();
        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(Allele.create("A")));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(sampleList, alleleList, readsBySample);

        final VariantContext vc= makeVC();
        final ReferenceContext referenceContext= null;
        final Map<String, Object> annotate = new AS_RMSMappingQuality().annotateRawData(referenceContext, vc, likelihoods);
        final String[] split =((String)annotate.get(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY)).split(AnnotationUtils.ALLELE_SPECIFIC_SPLIT_REGEX);
        Assert.assertEquals(split.length, 2);
        Assert.assertEquals(split[0], String.format("%.2f", 0.0));
        Assert.assertEquals(split[1], String.format("%.2f", 0.0));
    }

    @Test
    public void testCombineAndFinalize() {
        final List<Allele> vcAlleles = Arrays.asList(Allele.create("A", true), Allele.create("T", false));
        final List<ReducibleAnnotationData<?>> combinedVCdata = new ArrayList<>();
        combinedVCdata.add(new ReducibleAnnotationData<>("33640,10"));  //10 MQ58 reads
        combinedVCdata.add(new ReducibleAnnotationData<>("36000,10"));  //10 MQ60 reads

        final RMSMappingQuality annotator = RMSMappingQuality.getInstance();

        final Map<String, Object> combined = annotator.combineRawData(vcAlleles, combinedVCdata);
        final String combinedListString = (String)combined.get(annotator.getPrimaryRawKey());
        Assert.assertEquals(combinedListString, "69640,20");

        final VariantContext vc = new VariantContextBuilder(makeVC())
                .attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, combinedListString)
                .make();
        final VariantContext originalVC = null;
        final Map<String, Object> output = new RMSMappingQuality().finalizeRawData(vc, originalVC);
        Assert.assertEquals(Double.parseDouble((String)output.get("MQ")), 59.01);
    }

    @Test
    public void testCombineAndFinalizeHighMQSquared() {
        final List<Allele> vcAlleles = Arrays.asList(Allele.create("A", true), Allele.create("T", false));
        final List<ReducibleAnnotationData<?>> combinedVCdata = new ArrayList<>();
        // Test that updating the annotation works when Integer.MAX_VALUE is exceeded, both for small and large updates:
        combinedVCdata.add(new ReducibleAnnotationData<>("10125000000,5000000"));  //5,000,000 MQ45 reads
        combinedVCdata.add(new ReducibleAnnotationData<>("2601,1"));  //1 MQ51 read
        combinedVCdata.add(new ReducibleAnnotationData<>("1800000000,500000"));  //500,000 MQ60 reads

        final RMSMappingQuality annotator = RMSMappingQuality.getInstance();

        final Map<String, Object> combined = annotator.combineRawData(vcAlleles, combinedVCdata);
        final String combinedListString = (String)combined.get(annotator.getPrimaryRawKey());
        Assert.assertEquals(combinedListString, "11925002601,5500001");

        final VariantContext vc = new VariantContextBuilder(makeVC())
                .attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, combinedListString)
                .make();
        final VariantContext originalVC = null;
        final Map<String, Object> output = new RMSMappingQuality().finalizeRawData(vc, originalVC);
        Assert.assertEquals(Double.parseDouble((String)output.get("MQ")), 46.56);
    }

    @Test
    public void testFinalizeRawData(){
        // NOTE: RMSMappingQuality should ignore homRef depth of 20 and use only the 13 variants to calculate score.
        final VariantContext vc = new VariantContextBuilder(makeVC())
                .attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, "43732,13")
                .attribute(VCFConstants.DEPTH_KEY, 20)
                .make();
        final VariantContext originalVC = null;
        final Map<String, Object> output = new RMSMappingQuality().finalizeRawData(vc, originalVC);
        Assert.assertEquals(output.get("MQ"), "58.00");
    }

    @Test
    public void testFinalizeHighMQSquaredRawData(){
        // Test that RMS Mapping Quality is correctly computed when Integer.MAX_VALUE is exceeded.
        final VariantContext vc = new VariantContextBuilder(makeVC())
                .attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, "3415207168,1749038")
                .make();
        final VariantContext originalVC = null;
        final Map<String, Object> output = new RMSMappingQuality().finalizeRawData(vc, originalVC);
        Assert.assertEquals(output.get("MQ"), "44.19");
    }

    @Test
    public void testFinalizeRawMQ(){
        final VariantContext vc = new VariantContextBuilder(makeVC())
                .attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, "2000,20")
                .attribute(VCFConstants.DEPTH_KEY, 20)
                .make();
        final VariantContext output = new RMSMappingQuality().finalizeRawMQ(vc);
        Assert.assertFalse(output.hasAttribute( GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        Assert.assertTrue(output.hasAttribute(VCFConstants.RMS_MAPPING_QUALITY_KEY));
        Assert.assertEquals(output.getAttributeAsDouble(VCFConstants.RMS_MAPPING_QUALITY_KEY, -1.0), 10.0, 0.01);
    }


    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBadRawMQ(){
        final VariantContext vc = new VariantContextBuilder(makeVC())
                .attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, "fiftyeight")
                .make();
        new RMSMappingQuality().finalizeRawMQ(vc);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoDepth(){
        final VariantContext vc = new VariantContextBuilder(makeVC())
                .attribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, "20000")
                .make();
        new RMSMappingQuality().finalizeRawMQ(vc);
    }

    @Test
    public void testNoRawMQ(){
        final String OTHER_KEY = "SOME_OTHER_ANNOTATION";
        final String VALUE = "value";
        final VariantContext vc = new VariantContextBuilder(makeVC())
                    .attribute(OTHER_KEY, VALUE)
                    .make();
        final VariantContext output = new RMSMappingQuality().finalizeRawMQ(vc);
        Assert.assertFalse(output.hasAttribute(VCFConstants.RMS_MAPPING_QUALITY_KEY));
        Assert.assertFalse(output.hasAttribute(GATKVCFConstants.RAW_MAPPING_QUALITY_WITH_DEPTH_KEY));
        Assert.assertEquals(output.getAttribute(OTHER_KEY), VALUE);
    }

    @Test
    public void testGetNumReads(){
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Genotype homRefDP = new GenotypeBuilder("sample1", Arrays.asList(refAllele, refAllele)).DP(5).make();
        final Genotype homRefMinDP = new  GenotypeBuilder("sample2", Arrays.asList(refAllele, refAllele))
                .DP(11)
                .attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, 10).make();
        final Genotype hetAlt = new GenotypeBuilder("sample3", Arrays.asList(refAllele, altAllele)).DP(100).make();
        final VariantContext vc = new VariantContextBuilder().alleles(Arrays.asList(refAllele, altAllele))
                .chr("1")
                .start(15L)
                .stop(15L)
                .genotypes(homRefDP, homRefMinDP, hetAlt)
                .attribute(VCFConstants.DEPTH_KEY, 5 + 10 + 100)
                .make();

        Assert.assertEquals(RMSMappingQuality.getNumOfReads(vc, null), 115 - 5 - 10);
    }

    @Test
    public void testGetNumReadsReturnsNegativeOneWhenDataIsBad(){
        final VariantContext emptyVariantContext = makeVC();
        Assert.assertEquals(RMSMappingQuality.getNumOfReads(emptyVariantContext, null), -1);

    }

    @Test
    public void testGetNumReadsFallbackBehavior() {
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final VariantContext vc = new VariantContextBuilder().alleles(Arrays.asList(refAllele, altAllele))
                .chr("1").start(15L).stop(15L).make();

        final List<GATKRead> refReads = IntStream.range(0, 5).mapToObj(i -> ArtificialAnnotationUtils.makeRead(30, 5)).collect(Collectors.toList());
        final List<GATKRead> altReads = IntStream.range(0, 10).mapToObj(i -> ArtificialAnnotationUtils.makeRead(30, 5)).collect(Collectors.toList());
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = ArtificialAnnotationUtils.makeLikelihoods("sample1", refReads, altReads, -100, -100 ,refAllele, altAllele);

        // Testing the likelihoods map
        Assert.assertEquals(RMSMappingQuality.getNumOfReads(vc, likelihoods), 15);

        // Testing that unavailable mapping quality gets filtered out
        altReads.get(0).setMappingQuality(QualityUtils.MAPPING_QUALITY_UNAVAILABLE);
        Assert.assertEquals(RMSMappingQuality.getNumOfReads(vc, likelihoods), 14);
    }


    @Test
    public void testGetNumReadsFallbackReturnsNegativeOneWhenDataIsBad(){
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final VariantContext vc = new VariantContextBuilder().alleles(Arrays.asList(refAllele, altAllele))
                .chr("1").start(15L).stop(15L).make();
        final List<GATKRead> refReads = Collections.emptyList();
        final List<GATKRead> altReads = Collections.emptyList();
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = ArtificialAnnotationUtils.makeLikelihoods("sample1", refReads, altReads, -100, -100 ,refAllele, altAllele);
        Assert.assertEquals(RMSMappingQuality.getNumOfReads(vc,likelihoods), -1);
    }
}
