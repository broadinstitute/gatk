package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public final class VariantAnnotatorEngineUnitTest extends BaseTest {
    @Test
    public void testEmpty(){
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse= Collections.emptyList();
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testExclude(){
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<String> annotationsToExclude= annotationsToUse;
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
        Assert.assertTrue(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertTrue(vae.getInfoAnnotations().isEmpty());
    }

    @Test
    public void testIncludeGroupExcludeIndividual(){
        final List<String> annotationGroupsToUse= Collections.singletonList(StandardAnnotation.class.getSimpleName());
        final List<String> annotationsToUse = Collections.emptyList();
        final List<String> annotationsToExclude= Arrays.asList(Coverage.class.getSimpleName());
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());

        //check that Coverage is out
        Assert.assertTrue(vae.getInfoAnnotations().stream().noneMatch(a -> a.getClass().getSimpleName().equals(Coverage.class.getSimpleName())));
    }

    @Test
    public void testAll(){
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features);
        Assert.assertFalse(vae.getVCFAnnotationDescriptions().contains(null));
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());

        final List<GenotypeAnnotation> knowGenoAnnos = ClassUtils.makeInstancesOfSubclasses(GenotypeAnnotation.class, Annotation.class.getPackage());
        final List<InfoFieldAnnotation> knowInfoAnnos = ClassUtils.makeInstancesOfSubclasses(InfoFieldAnnotation.class, Annotation.class.getPackage());
        Assert.assertEquals(vae.getGenotypeAnnotations().size(), knowGenoAnnos.size());
        Assert.assertEquals(vae.getInfoAnnotations().size(), knowInfoAnnos.size());

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions();
        Assert.assertFalse(vcfAnnotationDescriptions.isEmpty());

        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY)));          //yes DP
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY)));   //yes AC
    }

    @Test
    public void testAllMinusCoverage(){
        final List<String> annotationsToExclude= Arrays.asList(Coverage.class.getSimpleName());
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features);
        Assert.assertFalse(vae.getVCFAnnotationDescriptions().contains(null));
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions();
        Assert.assertFalse(vcfAnnotationDescriptions.isEmpty());

        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));
        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY))); //no DP
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY)));  //yes AC
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadAnnot(){
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<String> annotationGroupsToUse = Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList("fred");
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullFeatures(){
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<String> annotationGroupsToUse = Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());
        final List<FeatureInput<VariantContext>> features = null;
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadExcludeAnnot(){
        final List<String> annotationsToExclude= Arrays.asList("fred");
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features);
        Assert.assertFalse(vae.getVCFAnnotationDescriptions().contains(null));
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testBadGroup(){
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<String> annotationGroupsToUse = Arrays.asList("fred");
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
        Assert.assertFalse(vae.getGenotypeAnnotations().isEmpty());
        Assert.assertFalse(vae.getInfoAnnotations().isEmpty());
    }

    private VariantContext makeVC(final Allele refAllele, final Allele altAllele) {
        return makeVC(refAllele, altAllele, new SimpleInterval("1", 15, 15));
    }
    private VariantContext makeVC(final Allele refAllele, final Allele altAllele, final Locatable loc) {
        final List<Allele> alleles = Arrays.asList(refAllele, altAllele);
        final Genotype g = new GenotypeBuilder("sample1", alleles).make();

        return (new VariantContextBuilder())
                .alleles(Arrays.asList(refAllele, altAllele)).chr(loc.getContig()).start(loc.getStart()).stop(loc.getEnd()).genotypes(g).make();
    }

    private Map<String, PerReadAlleleLikelihoodMap> makeReadMap(final int ref, final int alt, final Allele refAllele, final Allele altAllele) {
        return makeReadMap(ref, alt, refAllele, altAllele, "1", 10000);
    }

    private Map<String, PerReadAlleleLikelihoodMap> makeReadMap(final int ref,
                                                                final int alt,
                                                                final Allele refAllele,
                                                                final Allele altAllele,
                                                                final String contig,
                                                                final int position) {
        final PerReadAlleleLikelihoodMap map= new PerReadAlleleLikelihoodMap();

        for (int i = 0; i < alt; i++) {
            final GATKRead read = ArtificialReadUtils.createUniqueArtificialRead("10M");
            read.setPosition(contig, position);
            map.add(read, altAllele, -1.0);
            map.add(read, refAllele, -100.0);
        }
        for (int i = 0; i < ref; i++) {
            final GATKRead read = ArtificialReadUtils.createUniqueArtificialRead("10M");
            read.setPosition(contig, position);
            map.add(read, altAllele, -100.0);
            map.add(read, refAllele, -1.0);
        }

        return Collections.singletonMap("sample1", map);
    }

    @Test
    public void testAllAnnotations() throws Exception {
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features);
        Assert.assertFalse(vae.getVCFAnnotationDescriptions().contains(null));

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(makeVC(refAllele, altAllele), new FeatureContext(), null, perReadAlleleLikelihoodMap, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_COUNT_KEY), 1);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY), 1.0/2);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.ALLELE_NUMBER_KEY), 2);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.FISHER_STRAND_KEY), FisherStrand.makeValueObjectForAnnotation(new int[][]{{ref,0},{alt,0}}));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.NOCALL_CHROM_KEY), 0);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.RMS_MAPPING_QUALITY_KEY), RMSMappingQuality.formatedValue(0.0));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.MAPPING_QUALITY_ZERO_KEY), MappingQualityZero.formattedValue(8));
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.SAMPLE_LIST_KEY), "sample1");
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(GATKVCFConstants.STRAND_ODDS_RATIO_KEY), StrandOddsRatio.formattedValue(StrandOddsRatio.calculateSOR(new int[][]{{ref,0},{alt,0}})));
    }

    @Test
    public void testMultipleAnnotations() throws Exception {
        final List<String> annotationsToExclude = Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<String> annotationGroupsToUse = Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName(), FisherStrand.class.getSimpleName());
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(makeVC(refAllele, altAllele), new FeatureContext(), null, perReadAlleleLikelihoodMap, a->true);
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
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<String> annotationsToExclude= Collections.emptyList();
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>("dbsnp", Collections.emptyMap(), file);

        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions();
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));

        final int alt = 5;
        final int ref = 3;

        final SimpleInterval loc= new SimpleInterval("1", 69428, 69428);
        final VariantContext vcRS = new FeatureDataSource<>(file, new VCFCodec(), null, 0).query(loc).next();

        final Allele refAllele = vcRS.getReference();
        final Allele altAllele = vcRS.getAlternateAllele(0);

        final VariantContext vcToAnnotate = makeVC(refAllele, altAllele, loc);
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele, loc.getContig(), loc.getStart() - 5);

        final FeatureContext featureContext= when(mock(FeatureContext.class).getValues(dbSNPBinding, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcRS)).getMock();
        final VariantContext resultVC = vae.annotateContext(vcToAnnotate, featureContext, null, perReadAlleleLikelihoodMap, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));
        Assert.assertEquals(resultVC.getID(), vcRS.getID());
    }

    @Test
    public void testCoverageAnnotationOnOverlapSite() throws Exception {
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<String> annotationsToExclude= Collections.emptyList();
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = null;

        final String featureSourceName = "fred";
        final FeatureInput<VariantContext> fredInput = new FeatureInput<>(featureSourceName, Collections.emptyMap(), file);//we'll just reuse the DBSnp file under a different name
        final List<FeatureInput<VariantContext>> features = Arrays.asList(fredInput);

        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions();
        Assert.assertFalse(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));

        final VCFInfoHeaderLine fredHeaderLine = new VCFInfoHeaderLine(featureSourceName, 0, VCFHeaderLineType.Flag, featureSourceName + " Membership");
        Assert.assertTrue(vcfAnnotationDescriptions.contains(fredHeaderLine));
        final int alt = 5;
        final int ref = 3;

        final SimpleInterval loc= new SimpleInterval("1", 69428, 69428);
        final VariantContext vcRS = new FeatureDataSource<>(file, new VCFCodec(), null, 0).query(loc).next();

        final Allele refAllele = vcRS.getReference();
        final Allele altAllele = vcRS.getAlternateAllele(0);

        final VariantContext vcToAnnotate = makeVC(refAllele, altAllele, loc);
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele, loc.getContig(), loc.getStart() - 5);

        final FeatureContext featureContext = when(mock(FeatureContext.class).getValues(fredInput, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcRS)).getMock();

        final VariantContext resultVC = vae.annotateContext(vcToAnnotate, featureContext, null, perReadAlleleLikelihoodMap, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));

        Assert.assertEquals(resultVC.getID(), VCFConstants.EMPTY_ID_FIELD);
        Assert.assertTrue((boolean)resultVC.getCommonInfo().getAttribute(featureSourceName));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute("does not exist"));
    }

    @Test
    public void testCoverageAnnotationOnDBSNPAndOverlapSite() throws Exception {
        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<String> annotationsToExclude= Collections.emptyList();
        final File dbSNPFile= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>("dbsnp", Collections.emptyMap(), dbSNPFile);

        final File fredFile = getTestFile("one_entry_source.vcf");
        final String featureSourceName = "fred";
        final FeatureInput<VariantContext> fredInput = new FeatureInput<>(featureSourceName, Collections.emptyMap(), fredFile);
        final List<FeatureInput<VariantContext>> features = Arrays.asList(fredInput);

        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);

        final Set<VCFHeaderLine> vcfAnnotationDescriptions = vae.getVCFAnnotationDescriptions();
        Assert.assertTrue(vcfAnnotationDescriptions.contains(VCFStandardHeaderLines.getInfoLine(VCFConstants.DBSNP_KEY)));
        final VCFInfoHeaderLine fredHeaderLine = new VCFInfoHeaderLine(featureSourceName, 0, VCFHeaderLineType.Flag, featureSourceName + " Membership");
        Assert.assertTrue(vcfAnnotationDescriptions.contains(fredHeaderLine));

        final VCFInfoHeaderLine headerLine = new VCFInfoHeaderLine(featureSourceName, 0, VCFHeaderLineType.Flag, featureSourceName + " Membership");
        Assert.assertTrue(vcfAnnotationDescriptions.contains(headerLine));
        final int alt = 5;
        final int ref = 3;

        final SimpleInterval loc= new SimpleInterval("1", 69428, 69428);
        final VariantContext vcDbSNP = new FeatureDataSource<>(dbSNPFile, new VCFCodec(), null, 0).query(loc).next();
        final VariantContext vcFred = new FeatureDataSource<>(fredFile, new VCFCodec(), null, 0).query(loc).next();

        final Allele refAllele = vcDbSNP.getReference();
        final Allele altAllele = vcDbSNP.getAlternateAllele(0);

        final VariantContext vcToAnnotate = makeVC(refAllele, altAllele, loc);
        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele, loc.getContig(), loc.getStart() - 5);

        //both features
        final FeatureContext featureContext0 = when(mock(FeatureContext.class).getValues(dbSNPBinding, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcDbSNP)).getMock();
        final FeatureContext featureContext = when(featureContext0.getValues(fredInput, loc.getStart())).thenReturn(Arrays.<VariantContext>asList(vcFred)).getMock();

        final VariantContext resultVC = vae.annotateContext(vcToAnnotate, featureContext, null, perReadAlleleLikelihoodMap, a->true);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));

        //Check that if has both the DBSNP and Fred annotations
        Assert.assertEquals(resultVC.getID(), vcDbSNP.getID());
        Assert.assertTrue((boolean)resultVC.getCommonInfo().getAttribute(featureSourceName));
        Assert.assertNull(resultVC.getCommonInfo().getAttribute("does not exist"));
    }

    @Test(expectedExceptions = GATKException.class)
    public void testDBSNPONlyViaSpecialArg() throws Exception {
        final List<String> annotationGroupsToUse = Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName());//good one
        final List<String> annotationsToExclude = Collections.emptyList();
        final File dbSNPFile = new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>(VCFConstants.DBSNP_KEY, Collections.emptyMap(), dbSNPFile);

        final File fredFile = getTestFile("one_entry_source.vcf");
        final FeatureInput<VariantContext> fredInput = new FeatureInput<>(VCFConstants.DBSNP_KEY, Collections.emptyMap(), fredFile);
        final List<FeatureInput<VariantContext>> features = Arrays.asList(fredInput);

        VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);
    }

    @Test
    public void testCoverageAnnotationViaEngine() throws Exception {
        final File file= new File(publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf");
        final FeatureInput<VariantContext> dbSNPBinding = new FeatureInput<>("dbsnp", Collections.emptyMap(), file);

        final List<String> annotationGroupsToUse= Collections.emptyList();
        final List<String> annotationsToUse = Arrays.asList(Coverage.class.getSimpleName(),
                                                            DepthPerAlleleBySample.class.getSimpleName(),
                                                            SampleList.class.getSimpleName());
        final List<String> annotationsToExclude= Collections.emptyList();
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();
        final VariantAnnotatorEngine vae = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbSNPBinding, features);

        final int alt = 5;
        final int ref = 3;
        final Allele refAllele = Allele.create("A", true);
        final Allele altAllele = Allele.create("T");

        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap = makeReadMap(ref, alt, refAllele, altAllele);
        final VariantContext resultVC = vae.annotateContext(makeVC(refAllele, altAllele), new FeatureContext(), null, perReadAlleleLikelihoodMap,
                                            ann -> ann instanceof Coverage || ann instanceof DepthPerAlleleBySample);
        Assert.assertEquals(resultVC.getCommonInfo().getAttribute(VCFConstants.DEPTH_KEY), String.valueOf(ref+alt));

        Assert.assertEquals(resultVC.getGenotype(0).getAD(), new int[]{ref, alt});

        //skipped because we only asked for Coverage and DepthPerAlleleBySample
        Assert.assertFalse(resultVC.getCommonInfo().hasAttribute(GATKVCFConstants.SAMPLE_LIST_KEY));
    }

    @Test
    public void testAnnotationsHaveDescriptions() throws Exception {

        Set<String> sampleSet = Collections.singleton("FRED");
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();

        headerInfo.addAll(VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features).getVCFAnnotationDescriptions());

        Assert.assertFalse(headerInfo.contains(null));
        new VCFHeader(headerInfo, sampleSet);//make sure this does not blow up: https://github.com/broadinstitute/gatk/issues/1713
    }

    @Test
    public void testNoNullInKeysAndDescriptions() throws Exception {
        final List<String> annotationsToExclude= Collections.emptyList();
        final FeatureInput<VariantContext> dbSNPBinding = null;
        final List<FeatureInput<VariantContext>> features = Collections.emptyList();

        final VariantAnnotatorEngine variantAnnotatorEngine = VariantAnnotatorEngine.ofAllMinusExcluded(annotationsToExclude, dbSNPBinding, features);
        for (GenotypeAnnotation ga : variantAnnotatorEngine.getGenotypeAnnotations()) {
            Assert.assertFalse(ga.getDescriptions().contains(null), "getDescriptions contains null:" + ga);
            Assert.assertFalse(ga.getKeyNames().contains(null), "getKeyNames contains null" + ga);
        }
        for (InfoFieldAnnotation ifa : variantAnnotatorEngine.getInfoAnnotations()) {
            Assert.assertFalse(ifa.getDescriptions().contains(null), "getDescriptions contains null:" + ifa);
            Assert.assertFalse(ifa.getKeyNames().contains(null), "getKeyNames contains null:" + ifa);
        }
        Assert.assertFalse(variantAnnotatorEngine.getVCFAnnotationDescriptions().contains(null));
    }
}
