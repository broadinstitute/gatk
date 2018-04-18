package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public final class VariantAnnotatorEngineUnitTest extends GATKBaseTest {
    @Test
    public void testCombineAnnotations() throws Exception {
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(Collections.singletonList(new AS_RMSMappingQuality()), dbSNPBinding, features, false);

        final Allele REF = Allele.create("A", true);
        final Allele ALT = Allele.create("T");
        final List<Allele> alleles = Arrays.asList(new Allele[]{REF, ALT});


        final int[] MQs = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, QualityUtils.MAPPING_QUALITY_UNAVAILABLE};
        final int[] MQ2s = {1, 2, 3, 4, 5, 6, 7, 8, 9, QualityUtils.MAPPING_QUALITY_UNAVAILABLE};
        final List<Integer> MQsList = Arrays.asList(ArrayUtils.toObject(MQs));


        final List<GATKRead> altReads = Arrays.stream(MQs).mapToObj(mq -> ArtificialAnnotationUtils.makeRead(30, mq)).collect(Collectors.toList());
        final List<GATKRead> refReads = Arrays.stream(MQ2s).mapToObj(mq -> ArtificialAnnotationUtils.makeRead(30, mq)).collect(Collectors.toList());

        final ReadLikelihoods<Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods("sample1", refReads, altReads, -10.0, -9.0, REF, ALT);

        final VariantContext vc = ArtificialAnnotationUtils.makeVC();
        final ReducibleAnnotationData<Object> myData = new ReducibleAnnotationData<>("-10.0");

        new AS_RMSMappingQuality().calculateRawData(vc, likelihoods, myData);
        Map<String, List<?>> testAnnotationData = new HashMap<>();
        testAnnotationData.put(new AS_RMSMappingQuality().getRawKeyName(), Collections.singletonList(myData));

        Map<String, Object> value = vae.combineAnnotations(alleles, testAnnotationData);
        Assert.assertEquals(value.get(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY), "285.00|385.00");
    }

    @Test
    public void testFinalizeAnnotations() throws Exception {
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(Collections.singletonList(new AS_RMSMappingQuality()), dbSNPBinding, features, false);
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");
        final Genotype genotype = new  GenotypeBuilder("sample2", Arrays.asList(refAllele, altAllele))
                .AD(new int[]{8,9}).make();

        final VariantContext vc = new VariantContextBuilder(ArtificialAnnotationUtils.makeVC()).attribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, "285.00|385.00").attribute(GATKVCFConstants.ALLELE_BALANCE_KEY, "10").genotypes(genotype).make();
        final VariantContext result = vae.finalizeAnnotations(vc, vc);
        Assert.assertNull(result.getAttribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY));
        Assert.assertNotNull(result.getAttribute(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY));
        Assert.assertNotNull(result.getAttribute(GATKVCFConstants.ALLELE_BALANCE_KEY));
        Assert.assertEquals(result.getAttribute(GATKVCFConstants.AS_RMS_MAPPING_QUALITY_KEY), "6.54");
    }


    @Test
    public void testEmpty(){
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(Collections.emptyList(), dbSNPBinding, features, false);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

//    @Test
//    public void testExclude(){
//        final List<String> annotationGroupsToUse= Collections.emptyList();
//        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
//        final List<String> annotationsToExclude= annotationsToUse;
//        final FeatureInput<VariantContext> dbSNPBinding = null;
//        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
//        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(Collections.singletonList(new AS_RMSMappingQuality()), dbSNPBinding, features, false);
//        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
//        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
//    }

//    @Test
//    public void testIncludeGroupExcludeIndividual(){
//        final List<String> annotationGroupsToUse= Collections.singletonList(StandardAnnotation.class.getSimpleName());
//        final List<String> annotationsToUse = Collections.emptyList();
//        final List<String> annotationsToExclude= Arrays.asList(Coverage.class.getSimpleName());
//        final FeatureInput<VariantContext> dbSNPBinding = null;
//        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
//        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
//        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
//        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
//
//        //check that Coverage is out
//        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
//    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullFeatures(){
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<Annotation> annotationsToUse = Arrays.asList(new Coverage());
        final List<FeatureInput<VariantContext>> features = null;
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine( annotationsToUse, dbSNPBinding, features, false);
    }

//    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
//    public void testBadExcludeAnnot(){
//        final List<String> annotationsToExclude= Arrays.asList("fred");
//        final FeatureInput<VariantContext> dbSNPBinding = null;
//        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
//        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features, false);
//        Assert.assertFalse(vae.getVCFAnnotationDescriptions(false).contains(null));
//        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
//        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
//    }
//
//    @Test(expectedExceptions = CommandLineException.BadArgumentValue.class)
//    public void testBadGroup(){
//        final List<String> annotationsToExclude= Collections.emptyList();
//        final FeatureInput<VariantContext> dbSNPBinding = null;
//        final List<String> annotationGroupsToUse = Arrays.asList("fred");
//        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
//        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
//        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
//        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
//        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
//    }

    private VariantContext makeVC(final Allele refAllele, final Allele altAllele) {
        return makeVC(refAllele, altAllele, new SimpleInterval("1", 15, 15));
    }
    private VariantContext makeVC(final Allele refAllele, final Allele altAllele, final Locatable loc) {
        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final Genotype g = new GenotypeBuilder("sample1", alleles).make();

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr(loc.getContig()).start(loc.getStart()).stop(loc.getEnd()).genotypes(g).make();
    }

    private ReadLikelihoods<Allele> makeReadLikelihoods(final int ref, final int alt, final Allele refAllele, final Allele altAllele) {
        return makeReadLikelihoods(ref, alt, refAllele, altAllele, "1", 10000);
    }

    private ReadLikelihoods<Allele> makeReadLikelihoods(final int ref,
                                                        final int alt,
                                                        final Allele refAllele,
                                                        final Allele altAllele,
                                                        final String contig,
                                                        final int position) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < alt; i++) {
            final GATKRead read = ArtificialReadUtils.createUniqueArtificialRead("10M");
            read.setPosition(contig, position);
            reads.add(read);
        }
        for (int i = 0; i < ref; i++) {
            final GATKRead read = ArtificialReadUtils.createUniqueArtificialRead("10M");
            read.setPosition(contig, position);
            reads.add(read);
        }

        final Map<String, List<GATKRead>> readsBySample = ImmutableMap.of("sample1", reads);
        final org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList = new IndexedSampleList(Arrays.asList("sample1"));
        final AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(refAllele, altAllele));
        final ReadLikelihoods<Allele> likelihoods = new ReadLikelihoods<>(sampleList, alleleList, readsBySample);

        // modify likelihoods in-place
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(0);

        int n = 0;
        for (int i = 0; i < alt; i++) {
            matrix.set(0, n, -100.0);
            matrix.set(1, n, -1.0);
            n++;
        }
        for (int i = 0; i < ref; i++) {
            matrix.set(0, n, -1.0);
            matrix.set(1, n, -100.0);
            n++;
        }

        return likelihoods;
    }

    @Test
    public void testAllAnnotations() throws Exception {
        /**
         * exclude {@link ReferenceBases} until https://github.com/broadinstitute/gatk/issues/2799 is fixed
         * */
        final List<Class<? extends Annotation>> annotationsToExclude= Arrays.asList( ReferenceBases.class);
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features, false);
        Assert.assertFalse(vae.getVCFAnnotationDescriptions(false).contains(null));

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final ReadLikelihoods<Allele> likelihoods = makeReadLikelihoods(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(makeVC(refAllele, altAllele), new FeatureContext(), null, likelihoods, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_COUNT_KEY), 1);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY), 1.0/2);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_NUMBER_KEY), 2);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.FISHER_STRAND_KEY), FisherStrand.makeValueObjectForAnnotation(new int[][]{{ref,0},{alt,0}}));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.NOCALL_CHROM_KEY), 0);
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.RMS_MAPPING_QUALITY_KEY), RMSMappingQuality.formattedValue(0.0));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.MAPPING_QUALITY_ZERO_KEY), MappingQualityZero.formattedValue(8));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.SAMPLE_LIST_KEY), "sample1");
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.STRAND_ODDS_RATIO_KEY), StrandOddsRatio.formattedValue(StrandOddsRatio.calculateSOR(new int[][]{{ref,0},{alt,0}})));
    }

    @Test
    public void testMultipleAnnotations() throws Exception {
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<Annotation> annotationsToUse = Arrays.asList(new Coverage(), new FisherStrand());
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annotationsToUse, dbSNPBinding, features, false);

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final ReadLikelihoods<Allele> likelihoods = makeReadLikelihoods(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(makeVC(refAllele, altAllele), new FeatureContext(), null, likelihoods, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.FISHER_STRAND_KEY), FisherStrand.makeValueObjectForAnnotation(new int[][]{{ref,0},{alt,0}}));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_COUNT_KEY));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.NOCALL_CHROM_KEY));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(VCFConstants.RMS_MAPPING_QUALITY_KEY));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(VCFConstants.MAPPING_QUALITY_ZERO_KEY));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.SAMPLE_LIST_KEY));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.STRAND_ODDS_RATIO_KEY));
    }

    @Test
    public void testCoverageAnnotationOnDbSnpSite() throws Exception {
        final List<Annotation> annotationsToUse = Arrays.asList(new Coverage());//good one
        final String path = publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf";
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>(path, "dbsnp", Collections.emptyMap());

        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annotationsToUse, dbSNPBinding, features, false);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions(false);
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));

        final int alt = 5;
        final int ref = 3;

        final SimpleInterval loc= new SimpleInterval("1", 69428, 69428);
        final VariantContext vcRS = new FeatureDataSource<VariantContext>(path, null, 0, VariantContext.class).query(loc).next();

        final Allele refAllele = vcRS.getReference();
        final Allele altAllele = vcRS.getAlternateAllele(0);

        final VariantContext vcToAnnotate = makeVC(refAllele, altAllele, loc);
        final ReadLikelihoods<Allele> likelihoods = makeReadLikelihoods(ref, alt, refAllele, altAllele, loc.getContig(), loc.getStart() - 5);

        final FeatureContext featureContext= when(mock(FeatureContext.class).getValues(dbSNPBinding, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcRS)).getMock();
        final VariantContext resultVC = vae.annotateContext(vcToAnnotate, featureContext, null, likelihoods, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));
        Assert.assertEquals(resultVC.getID(), vcRS.getID());
    }

    @Test
    public void testCoverageAnnotationOnOverlapSite() throws Exception {
        final List<Annotation> annotationsToUse = Arrays.asList(new Coverage());//good one
        final String path = publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf";
        final FeatureInput<VariantContext> dbSNPBinding = null;

        final String featureSourceName = "fred";
        final FeatureInput<VariantContext> fredInput = new FeatureInput<>(path, featureSourceName, Collections.emptyMap());//we'll just reuse the DBSnp file under a different name
        final List<FeatureInput<VariantContext>> features = Arrays.asList(fredInput);

        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annotationsToUse, dbSNPBinding, features, false);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions(false);
        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));

        final VCFInfoHeaderLine fredHeaderLine = new VCFInfoHeaderLine(featureSourceName, 0, VCFHeaderLineType.Flag, featureSourceName + " Membership");
        Assert.assertTrue(vcfAnnotationDescriptions.contains(fredHeaderLine));
        final int alt = 5;
        final int ref = 3;

        final SimpleInterval loc= new SimpleInterval("1", 69428, 69428);
        final VariantContext vcRS = new FeatureDataSource<VariantContext>(path, null, 0, VariantContext.class).query(loc).next();

        final Allele refAllele = vcRS.getReference();
        final Allele altAllele = vcRS.getAlternateAllele(0);

        final VariantContext vcToAnnotate = makeVC(refAllele, altAllele, loc);
        final ReadLikelihoods<Allele> likelihoods = makeReadLikelihoods(ref, alt, refAllele, altAllele, loc.getContig(), loc.getStart() - 5);

        final FeatureContext featureContext = when(mock(FeatureContext.class).getValues(fredInput, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcRS)).getMock();

        final VariantContext resultVC = vae.annotateContext(vcToAnnotate, featureContext, null, likelihoods, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));

        Assert.assertEquals(resultVC.getID(), VCFConstants.EMPTY_ID_FIELD);
        Assert.assertTrue((boolean)resultVC.getCommonInfo().getAttribute(featureSourceName));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute("does not exist"));
    }

    @Test
    public void testCoverageAnnotationOnDBSNPAndOverlapSite() throws Exception {
        final List<Annotation> annotationsToUse = Arrays.asList(new Coverage());//good one
        final String dbSNPPath= publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf";
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>(dbSNPPath, "dbsnp", Collections.emptyMap());

        final File fredFile = getTestFile("one_entry_source.vcf");
        final String featureSourceName = "fred";
        final FeatureInput<VariantContext> fredInput = new FeatureInput<>(fredFile.getAbsolutePath(), featureSourceName, Collections.emptyMap());
        final List<FeatureInput<VariantContext>> features = Arrays.asList(fredInput);

        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annotationsToUse, dbSNPBinding, features, false);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions(false);
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));
        final VCFInfoHeaderLine fredHeaderLine = new VCFInfoHeaderLine(featureSourceName, 0, VCFHeaderLineType.Flag, featureSourceName + " Membership");
        Assert.assertTrue(vcfAnnotationDescriptions.contains(fredHeaderLine));

        final VCFInfoHeaderLine headerLine = new VCFInfoHeaderLine(featureSourceName, 0, VCFHeaderLineType.Flag, featureSourceName + " Membership");
        Assert.assertTrue(vcfAnnotationDescriptions.contains(headerLine));
        final int alt = 5;
        final int ref = 3;

        final SimpleInterval loc= new SimpleInterval("1", 69428, 69428);
        final VariantContext vcDbSNP = new FeatureDataSource<VariantContext>(dbSNPPath, null, 0, VariantContext.class).query(loc).next();
        final VariantContext vcFred = new FeatureDataSource<VariantContext>(fredFile.getAbsolutePath(), null, 0, VariantContext.class).query(loc).next();

        final Allele refAllele = vcDbSNP.getReference();
        final Allele altAllele = vcDbSNP.getAlternateAllele(0);

        final VariantContext vcToAnnotate = makeVC(refAllele, altAllele, loc);
        final ReadLikelihoods<Allele> likelihoods = makeReadLikelihoods(ref, alt, refAllele, altAllele, loc.getContig(), loc.getStart() - 5);

        //both features
        final FeatureContext featureContext0 = when(mock(FeatureContext.class).getValues(dbSNPBinding, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcDbSNP)).getMock();
        final FeatureContext featureContext = when(featureContext0.getValues(fredInput, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcFred)).getMock();

        final VariantContext resultVC = vae.annotateContext(vcToAnnotate, featureContext, null, likelihoods, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));

        //Check that if has both the DBSNP and Fred annotations
        Assert.assertEquals(resultVC.getID(), vcDbSNP.getID());
        Assert.assertTrue((boolean)resultVC.getCommonInfo().getAttribute(featureSourceName));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute("does not exist"));
    }

    @Test(expectedExceptions = GATKException.class)
    public void testDBSNPONlyViaSpecialArg() throws Exception {
        final List<Annotation> annotationsToUse = Arrays.asList(new Coverage());//good one
        final File dbSNPFile = new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>(dbSNPFile.getAbsolutePath(), VCFConstants.DBSNP_KEY, Collections.emptyMap());

        final File fredFile = getTestFile("one_entry_source.vcf");
        final FeatureInput<VariantContext> fredInput = new FeatureInput<>(fredFile.getAbsolutePath(), VCFConstants.DBSNP_KEY, Collections.emptyMap());
        final List<FeatureInput<VariantContext>> features = Arrays.asList(fredInput);

        new VariantAnnotatorEngine(annotationsToUse, dbSNPBinding, features, false);
    }

    @Test
    public void testCoverageAnnotationViaEngine() throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>(file.getAbsolutePath(), "dbsnp", Collections.emptyMap());

        final List<Annotation> annotationsToUse = Arrays.asList(new Coverage(),
                new DepthPerAlleleBySample(),
                new SampleList());
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = new VariantAnnotatorEngine(annotationsToUse, dbSNPBinding, features, false);

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final ReadLikelihoods<Allele> likelihoods = makeReadLikelihoods(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(makeVC(refAllele, altAllele), new FeatureContext(), null, likelihoods,
                ann -> ann instanceof Coverage || ann instanceof DepthPerAlleleBySample);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));

        Assert.assertEquals(resultVC.getGenotype(0).getAD(), new int[]{ref, alt});

        //skipped because we only asked for Coverage and DepthPerAlleleBySample
        Assert.assertFalse(resultVC.getCommonInfo().hasAttribute(GATKVCFConstants.SAMPLE_LIST_KEY));
    }

    @Test
    public void testAnnotationsHaveDescriptions() throws Exception {

        Set<String> sampleSet = Collections.singleton("FRED");
        final Set<VCFHeaderLine> headerInfo = new LinkedHashSet<>();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();

        headerInfo.addAll(VariantAnnotatorEngine.ofAllMinusExcluded(Collections.emptyList(), dbSNPBinding, features, false).getVCFAnnotationDescriptions(false));

        Assert.assertFalse(headerInfo.contains(null));
        new VCFHeader(headerInfo, sampleSet);//make sure this does not blow up: https://github.com/broadinstitute/gatk/issues/1713
    }

    @Test
    public void testNoNullInKeysAndDescriptions() throws Exception {
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();

        final VariantAnnotatorEngine variantAnnotatorEngine = VariantAnnotatorEngine.ofAllMinusExcluded(Collections.emptyList(), dbSNPBinding, features, false);
        for (GenotypeAnnotation ga : variantAnnotatorEngine.getGenotypeAnnotations()) {
            Assert.assertFalse(ga.getDescriptions().contains(null), "getDescriptions contains null:" + ga);
            Assert.assertFalse(ga.getKeyNames().contains(null), "getKeyNames contains null" + ga);
        }
        for (InfoFieldAnnotation ifa : variantAnnotatorEngine.getInfoAnnotations()) {
            Assert.assertFalse(ifa.getDescriptions().contains(null), "getDescriptions contains null:" + ifa);
            Assert.assertFalse(ifa.getKeyNames().contains(null), "getKeyNames contains null:" + ifa);
        }
        Assert.assertFalse(variantAnnotatorEngine.getVCFAnnotationDescriptions(false).contains(null));
    }

}
