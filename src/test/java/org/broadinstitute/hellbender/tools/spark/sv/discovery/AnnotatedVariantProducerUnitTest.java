package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.InsDelVariantDetector;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyReferenceLocations;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.mockito.Mockito;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple4;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.mockito.Mockito.when;


public class AnnotatedVariantProducerUnitTest extends GATKBaseTest {

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SVDiscoveryTestDataProvider.testDataInitialized) {
            new SVDiscoveryTestDataProvider();
        }
    }


    // -----------------------------------------------------------------------------------------------
    // Evidence summary annotation
    // -----------------------------------------------------------------------------------------------
    /**
     * Not an exhaustive test on all attributes, only tests:
     * MAPPING_QUALITIES, ALIGNMENT_LENGTH
     */
    private static void seeIfItWorks_evidenceAnnotation(final Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData,
                                                        final String[] expectedMappingQualitiesAsStrings,
                                                        final String[] expectedAlignmentLengthsAsStrings) {

        final AlignmentInterval region1 = testData._1();
        final AlignmentInterval region2 = testData._2();
        final byte[] contigSeq = null; // hack, as the contig sequence is really not necessary for this test purpose

        final Map<String, Object> attributeMap =
                AnnotatedVariantProducer.getEvidenceRelatedAnnotations(Collections.singletonList(new ChimericAlignment(region1, region2, Collections.emptyList(), testData._4(), SVDiscoveryTestDataProvider.seqDict)));

        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFConstants.MAPPING_QUALITIES)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedMappingQualitiesAsStrings);
        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFConstants.ALIGN_LENGTHS)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedAlignmentLengthsAsStrings);
    }

    @Test(groups = "sv")
    public void testGetEvidenceRelatedAnnotations() {

        // inversion
        Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData = SVDiscoveryTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(1984)});

        // simple deletion
        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple insertion
        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(100)});

        // long range substitution
        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple deletion with homology
        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple tandem dup contraction from 2 units to 1 unit
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple tandem dup expansion from 1 unit to 2 units
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(50)});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(137)});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(127)});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(31)});

        // tandem dup contraction from 3 units to 2 units
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(31)});

        // tandem dup expansion from 2 units to 3 units
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(127)});
    }

    // -----------------------------------------------------------------------------------------------
    // Integrative test
    // -----------------------------------------------------------------------------------------------
    private static void seeIfItWorks_integrative(final Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData,
                                                 final List<String> expectedAttributeKeys, final String sampleId) throws IOException {

        final AlignmentInterval region1 = testData._1();
        final AlignmentInterval region2 = testData._2();

        final Iterable<ChimericAlignment> evidence = Collections.singletonList(new ChimericAlignment(region1, region2, Collections.emptyList(), testData._4(), SVDiscoveryTestDataProvider.seqDict));

        final NovelAdjacencyReferenceLocations breakpoints = testData._3();

        final VariantContext variantContext =
                AnnotatedVariantProducer.produceAnnotatedVcFromInferredTypeAndRefLocations(breakpoints.leftJustifiedLeftRefLoc,
                        breakpoints.leftJustifiedRightRefLoc.getStart(), breakpoints.complication, InsDelVariantDetector.inferTypeFromNovelAdjacency(breakpoints),
                        null, evidence, SparkContextFactory.getTestSparkContext().broadcast(SVDiscoveryTestDataProvider.reference),
                        SparkContextFactory.getTestSparkContext().broadcast(SVDiscoveryTestDataProvider.seqDict), null, sampleId);

        final List<String> attributeKeys = variantContext.getAttributes().keySet().stream().sorted().collect(Collectors.toList());

        Assert.assertEquals(attributeKeys, expectedAttributeKeys);
    }

    @Test(groups = "sv")
    public void testIntegrative() throws IOException {

        final Set<String> commonAttributes = Sets.newHashSet(VCFConstants.END_KEY, GATKSVVCFConstants.SVLEN, GATKSVVCFConstants.SVTYPE,
                GATKSVVCFConstants.TOTAL_MAPPINGS, GATKSVVCFConstants.HQ_MAPPINGS, GATKSVVCFConstants.MAPPING_QUALITIES,
                GATKSVVCFConstants.ALIGN_LENGTHS, GATKSVVCFConstants.MAX_ALIGN_LENGTH, GATKSVVCFConstants.CONTIG_NAMES);

        final String sampleId = "sample";
        
        // inversion
        Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData = SVDiscoveryTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.INV33, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()), sampleId);

        // simple deletion
        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_minus;

        seeIfItWorks_integrative(testData, commonAttributes.stream().sorted().collect(Collectors.toList()), sampleId);

        // simple insertion
        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // long range substitution
        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // simple deletion with homology
        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // simple tandem dup contraction from 2 units to 1 unit
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH, GATKSVVCFConstants.DUP_INV_ORIENTATIONS).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // simple tandem dup expansion from 1 unit to 2 units
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUP_SEQ_CIGARS, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_INV_ORIENTATIONS).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUP_SEQ_CIGARS, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.INSERTED_SEQUENCE, GATKSVVCFConstants.DUP_INV_ORIENTATIONS).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH, GATKSVVCFConstants.DUP_INV_ORIENTATIONS).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH, GATKSVVCFConstants.DUP_INV_ORIENTATIONS).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // tandem dup contraction from 3 units to 2 units
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH, GATKSVVCFConstants.DUP_INV_ORIENTATIONS).stream())
                .sorted().collect(Collectors.toList()), sampleId);

        // tandem dup expansion from 2 units to 3 units
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH, GATKSVVCFConstants.DUP_INV_ORIENTATIONS).stream())
                .sorted().collect(Collectors.toList()), sampleId);
    }


    @DataProvider(name = "CIIntervals")
    public Object[][] getCIIntervalTests() {
        return new Object[][] {
                new Object[] { 200, new SVInterval(1, 190, 225), "-10,25", null},
                new Object[] { 200, new SVInterval(1, 200, 225), "0,25", null},
                new Object[] { 200, new SVInterval(1, 201, 225), null, new IllegalStateException("Interval must contain point")}
        };
    }

    @Test(dataProvider = "CIIntervals")
    public void testProduceCIInterval(final int point, final SVInterval interval, final String expected, final Exception expectedException) {
        if (expectedException == null) {
            Assert.assertEquals(AnnotatedVariantProducer.produceCIInterval(point, interval), expected);
        } else {
            try {
                AnnotatedVariantProducer.produceCIInterval(point, interval);
                Assert.fail("did not throw expected exception " + expectedException);
            } catch (Throwable e) {
                Assert.assertEquals(e.getClass(), expectedException.getClass());
                Assert.assertEquals(e.getMessage(), expectedException.getMessage());
            }
        }
    }


    // -----------------------------------------------------------------------------------------------
    // EvidenceTargetLink-based annotations
    // -----------------------------------------------------------------------------------------------
    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @DataProvider(name = "evidenceTargetLinksAndPreciseVariants")
    public Object[][] getEvidenceTargetLinksAndPreciseVariants() {

        final VariantContext unAnnotatedVC = new VariantContextBuilder()
                .id("TESTID")
                .chr("20").start(200).stop(300)
                .alleles("N", SimpleSVType.ImpreciseDeletion.createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL))
                .attribute(VCFConstants.END_KEY, 300)
                .attribute(GATKSVVCFConstants.SVTYPE, SimpleSVType.TYPES.DEL.toString())
                .make();

        final VariantContext annotatedVC = new VariantContextBuilder()
                .id("TESTID")
                .chr("20").start(200).stop(300)
                .alleles("N", SimpleSVType.ImpreciseDeletion.createBracketedSymbAlleleString(GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL))
                .attribute(VCFConstants.END_KEY, 300)
                .attribute(GATKSVVCFConstants.SVTYPE, SimpleSVType.TYPES.DEL.toString())
                .attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, 7)
                .attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, 5)
                .make();

        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[] {
                Arrays.asList(
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 190, 210), true),
                                new StrandedInterval(new SVInterval(0, 310, 320), false),
                                5, 7, new HashSet<>(), new HashSet<>())),
                Arrays.asList( unAnnotatedVC ),
                Arrays.asList( annotatedVC ) }
        );
        tests.add(new Object[] {
                Arrays.asList(
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 190, 210), true),
                                new StrandedInterval(new SVInterval(0, 310, 320), true),
                                5, 7, new HashSet<>(), new HashSet<>())),
                Arrays.asList( unAnnotatedVC ),
                Arrays.asList( unAnnotatedVC ) }
        );
        tests.add(new Object[] {
                Arrays.asList(
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 190, 210), true),
                                new StrandedInterval(new SVInterval(0, 310, 320), false),
                                3, 4, new HashSet<>(), new HashSet<>()),
                        new EvidenceTargetLink(
                                new StrandedInterval(new SVInterval(0, 192, 215), true),
                                new StrandedInterval(new SVInterval(0, 299, 303), false),
                                2, 3, new HashSet<>(), new HashSet<>())),
                Arrays.asList( unAnnotatedVC ),
                Arrays.asList( annotatedVC ) }
        );

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "evidenceTargetLinksAndPreciseVariants", groups = "sv")
    public void testProcessEvidenceTargetLinks(final List<EvidenceTargetLink> etls,
                                               final List<VariantContext> inputVariants,
                                               final List<VariantContext> expectedVariants) {

        final Logger localLogger = LogManager.getLogger(AnnotatedVariantProducer.class);
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection params =
                new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection();

        final ReferenceMultiSource referenceMultiSource = new ReferenceMultiSource(
                twoBitRefURL, ReferenceWindowFunctions.IDENTITY_FUNCTION);

        ReadMetadata metadata = Mockito.mock(ReadMetadata.class);
        when(metadata.getMaxMedianFragmentSize()).thenReturn(300);
        when(metadata.getContigName(0)).thenReturn("20");

        PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTree = new PairedStrandedIntervalTree<>();
        etls.forEach(e -> evidenceTree.put(e.getPairedStrandedIntervals(), e));

        final List<VariantContext> processedVariantContexts =
                AnnotatedVariantProducer.annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(inputVariants,
                        evidenceTree, metadata, referenceMultiSource, params, localLogger);

        VariantContextTestUtils.assertEqualVariants(processedVariantContexts, expectedVariants);
    }
}