package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyAndChimericAlignmentEvidence;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.mockito.Mockito;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;
import static org.mockito.Mockito.when;


public class AnnotatedVariantProducerUnitTest extends GATKBaseTest {

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SimpleSVDiscoveryTestDataProvider.testDataInitialized) {
            new SimpleSVDiscoveryTestDataProvider();
        }
    }


    // -----------------------------------------------------------------------------------------------
    // Evidence summary annotation
    // -----------------------------------------------------------------------------------------------
    /**
     * Not an exhaustive test on all attributes, only tests:
     * MAPPING_QUALITIES, ALIGNMENT_LENGTH
     */
    @Test(groups = "sv", dataProvider = "forEvidenceAnnotation")
    public void testEvidenceAnnotation(final TestDataForSimpleSVs testData,
                                       final String[] expectedMappingQualitiesAsStrings,
                                       final String[] expectedAlignmentLengthsAsStrings) {

        final List<SimpleNovelAdjacencyAndChimericAlignmentEvidence.SimpleChimeraAndNCAMstring> chimericAlignments = Collections.singletonList(
                new SimpleNovelAdjacencyAndChimericAlignmentEvidence.SimpleChimeraAndNCAMstring(
                        new SimpleChimera(testData.firstAlignment, testData.secondAlignment,
                                Collections.emptyList(), testData.evidenceAssemblyContigName, b37_seqDict),
                        AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME));

        final Map<String, Object> attributeMap =
                AnnotatedVariantProducer.getEvidenceRelatedAnnotations(chimericAlignments);

        Assert.assertEquals(((String)attributeMap.get(MAPPING_QUALITIES)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedMappingQualitiesAsStrings);
        Assert.assertEquals(((String)attributeMap.get(ALIGN_LENGTHS)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedAlignmentLengthsAsStrings);
    }

    @DataProvider(name = "forEvidenceAnnotation")
    private Object[][] dataForEvidenceAnnotation() {
        final List<Object[]> data = new ArrayList<>(20);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint, new String[]{"60"}, new String[]{String.valueOf(1984)}});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_minus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_plus, new String[]{"60"}, new String[]{String.valueOf(100)}});

        // long range substitution fudged del
        data.add(new Object[]{forLongRangeSubstitution_fudgedDel_minus, new String[]{"60"}, new String[]{String.valueOf(70)}});

        // long range substitution fat ins
        data.add(new Object[]{forLongRangeSubstitution_fatIns_plus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // long range substitution del+ins
        data.add(new Object[]{forLongRangeSubstitution_DelAndIns_plus, new String[]{"60"}, new String[]{String.valueOf(70)}});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_plus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_minus, new String[]{"60"}, new String[]{String.valueOf(40)}});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_ins_plus, new String[]{"60"}, new String[]{String.valueOf(50)}});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_dup_minus, new String[]{"60"}, new String[]{String.valueOf(137)}});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_plus, new String[]{"60"}, new String[]{String.valueOf(127)}});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_minus, new String[]{"60"}, new String[]{String.valueOf(31)}});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_plus, new String[]{"60"}, new String[]{String.valueOf(31)}});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_minus, new String[]{"60"}, new String[]{String.valueOf(127)}});

        return data.toArray(new Object[data.size()][]);
    }

    // -----------------------------------------------------------------------------------------------
    // Integrative test
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv", dataProvider = "forIntegrativeTest")
    public void testIntegrative(final TestDataForSimpleSVs testData,
                                final List<String> expectedAttributeKeys,
                                final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                final Broadcast<ReferenceMultiSource> referenceBroadcast,
                                final Broadcast<SAMSequenceDictionary> refSeqDictBroadcast) throws IOException {

        final NovelAdjacencyAndAltHaplotype breakpoints = testData.biPathBubble;
        final List<SimpleNovelAdjacencyAndChimericAlignmentEvidence.SimpleChimeraAndNCAMstring> evidence = Collections.singletonList(
                new SimpleNovelAdjacencyAndChimericAlignmentEvidence.SimpleChimeraAndNCAMstring(new SimpleChimera(testData.firstAlignment, testData.secondAlignment,
                        Collections.emptyList(), testData.evidenceAssemblyContigName, b37_seqDict),
                        AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME));
        final String sampleID = "testSample";

        final VariantContext variantContext =
                AnnotatedVariantProducer
                        .produceAnnotatedVcFromInferredTypeAndRefLocations(
                                breakpoints,
                                DiscoverVariantsFromContigAlignmentsSAMSpark.inferSimpleTypeFromNovelAdjacency(breakpoints),
                                evidence,
                                referenceBroadcast, refSeqDictBroadcast, broadcastCNVCalls, sampleID);

        final List<String> attributeKeys = variantContext.getAttributes().keySet().stream().sorted().collect(Collectors.toList());

        Assert.assertEquals(attributeKeys, expectedAttributeKeys);
    }

    @DataProvider(name = "forIntegrativeTest")
    private Object[][] dataForIntegrativeTest() {
        final List<Object[]> data = new ArrayList<>(20);

        final JavaSparkContext testSparkContext = SparkContextFactory.getTestSparkContext();
        final Broadcast<ReferenceMultiSource> referenceBroadcast = testSparkContext.broadcast(b37_reference);
        final Broadcast<SAMSequenceDictionary> refSeqDictBroadcast = testSparkContext.broadcast(b37_seqDict);

        final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls = null;

        final Set<String> commonAttributes = Sets.newHashSet(VCFConstants.END_KEY, SVLEN, SVTYPE, TOTAL_MAPPINGS,
                HQ_MAPPINGS, MAPPING_QUALITIES, ALIGN_LENGTHS, MAX_ALIGN_LENGTH, CONTIG_NAMES, SEQ_ALT_HAPLOTYPE);

        final Set<String> commAttributesWithoutAltSeq = Sets.newHashSet(VCFConstants.END_KEY, SVLEN, SVTYPE, TOTAL_MAPPINGS,
                HQ_MAPPINGS, MAPPING_QUALITIES, ALIGN_LENGTHS, MAX_ALIGN_LENGTH, CONTIG_NAMES);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint,
                Stream.concat( commAttributesWithoutAltSeq.stream(), Sets.newHashSet(INV33, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                 broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_minus,
                commAttributesWithoutAltSeq.stream().sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_plus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(INSERTED_SEQUENCE, INSERTED_SEQUENCE_LENGTH).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // long range substitution
        data.add(new Object[]{forLongRangeSubstitution_fudgedDel_minus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(INSERTED_SEQUENCE, INSERTED_SEQUENCE_LENGTH).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_plus,
                Stream.concat( commAttributesWithoutAltSeq.stream(), Sets.newHashSet(HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_minus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(DUP_TAN_CONTRACTION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_ins_plus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUP_SEQ_CIGARS, DUPLICATION_NUMBERS, DUP_ORIENTATIONS).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_dup_minus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUP_SEQ_CIGARS, DUPLICATION_NUMBERS, DUP_ORIENTATIONS, INSERTED_SEQUENCE, INSERTED_SEQUENCE_LENGTH).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_plus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH, DUP_IMPRECISE_AFFECTED_RANGE).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_minus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(DUP_TAN_CONTRACTION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH, DUP_IMPRECISE_AFFECTED_RANGE).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_plus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(DUP_TAN_CONTRACTION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH, DUP_IMPRECISE_AFFECTED_RANGE).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_minus,
                Stream.concat( commonAttributes.stream(), Sets.newHashSet(DUP_TAN_EXPANSION_STRING, DUP_REPEAT_UNIT_REF_SPAN, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_ORIENTATIONS, HOMOLOGY, HOMOLOGY_LENGTH, DUP_IMPRECISE_AFFECTED_RANGE).stream()).sorted().collect(Collectors.toList()),
                broadcastCNVCalls, referenceBroadcast, refSeqDictBroadcast});

        return data.toArray(new Object[data.size()][]);
    }

    // -----------------------------------------------------------------------------------------------
    // CI test
    // -----------------------------------------------------------------------------------------------
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

    @DataProvider(name = "evidenceTargetLinksAndPreciseVariants")
    public Object[][] getEvidenceTargetLinksAndPreciseVariants() {

        final VariantContext unAnnotatedVC = new VariantContextBuilder()
                .id("TESTID")
                .chr("20").start(200).stop(300)
                .alleles("N", SimpleSVType.ImpreciseDeletion.createBracketedSymbAlleleString(SYMB_ALT_ALLELE_DEL))
                .attribute(VCFConstants.END_KEY, 300)
                .attribute(SVTYPE, SimpleSVType.TYPES.DEL.toString())
                .make();

        final VariantContext annotatedVC = new VariantContextBuilder()
                .id("TESTID")
                .chr("20").start(200).stop(300)
                .alleles("N", SimpleSVType.ImpreciseDeletion.createBracketedSymbAlleleString(SYMB_ALT_ALLELE_DEL))
                .attribute(VCFConstants.END_KEY, 300)
                .attribute(SVTYPE, SimpleSVType.TYPES.DEL.toString())
                .attribute(READ_PAIR_SUPPORT, 7)
                .attribute(SPLIT_READ_SUPPORT, 5)
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

        ReadMetadata metadata = Mockito.mock(ReadMetadata.class);
        when(metadata.getMaxMedianFragmentSize()).thenReturn(300);
        when(metadata.getContigName(0)).thenReturn("20");

        PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTree = new PairedStrandedIntervalTree<>();
        etls.forEach(e -> evidenceTree.put(e.getPairedStrandedIntervals(), e));

        final List<VariantContext> processedVariantContexts =
                AnnotatedVariantProducer.annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(inputVariants,
                        evidenceTree, metadata, b37_reference, params, localLogger);

        VariantContextTestUtils.assertEqualVariants(processedVariantContexts, expectedVariants);
    }
}