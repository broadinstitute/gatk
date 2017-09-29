package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple4;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class AnnotatedVariantProducerUnitTest extends BaseTest {

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
                                                        final String[] expectedAlignmentLengthsAsStrings) throws IOException {

        final AlignmentInterval region1 = testData._1();
        final AlignmentInterval region2 = testData._2();
        final byte[] contigSeq = null; // hack, as the contig sequence is really not necessary for this test purpose

        final Map<String, Object> attributeMap =
                AnnotatedVariantProducer.getEvidenceRelatedAnnotations(Collections.singletonList(new ChimericAlignment(region1, region2, Collections.emptyList(), testData._4())));

        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFConstants.MAPPING_QUALITIES)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedMappingQualitiesAsStrings);
        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFConstants.ALIGN_LENGTHS)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                expectedAlignmentLengthsAsStrings);
    }

    @Test(groups = "sv")
    public void testGetEvidenceRelatedAnnotations() throws IOException {

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
                                                 final List<String> expectedAttributeKeys) throws IOException {

        final AlignmentInterval region1 = testData._1();
        final AlignmentInterval region2 = testData._2();

        final Iterable<ChimericAlignment> evidence = Collections.singletonList(new ChimericAlignment(region1, region2, Collections.emptyList(), testData._4()));

        final NovelAdjacencyReferenceLocations breakpoints = testData._3();

        final VariantContext variantContext =
                AnnotatedVariantProducer.produceAnnotatedVcFromInferredTypeAndRefLocations(breakpoints.leftJustifiedLeftRefLoc,
                        breakpoints.leftJustifiedRightRefLoc.getStart(), breakpoints.complication, SvTypeInference.inferFromNovelAdjacency(breakpoints),
                        evidence, SparkContextFactory.getTestSparkContext().broadcast(SVDiscoveryTestDataProvider.reference));

        final List<String> attributeKeys = variantContext.getAttributes().keySet().stream().sorted().collect(Collectors.toList());

        Assert.assertEquals(attributeKeys, expectedAttributeKeys);
    }

    @Test(groups = "sv")
    public void testIntegrative() throws IOException {

        final Set<String> commonAttributes = Sets.newHashSet(VCFConstants.END_KEY, GATKSVVCFConstants.SVLEN, GATKSVVCFConstants.SVTYPE,
                GATKSVVCFConstants.TOTAL_MAPPINGS, GATKSVVCFConstants.HQ_MAPPINGS, GATKSVVCFConstants.MAPPING_QUALITIES,
                GATKSVVCFConstants.ALIGN_LENGTHS, GATKSVVCFConstants.MAX_ALIGN_LENGTH, GATKSVVCFConstants.CONTIG_NAMES);

        // inversion
        Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData = SVDiscoveryTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.INV33, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()));

        // simple deletion
        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_minus;

        seeIfItWorks_integrative(testData, commonAttributes.stream().sorted().collect(Collectors.toList()));

        // simple insertion
        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()));

        // long range substitution
        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()));

        // simple deletion with homology
        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // simple tandem dup contraction from 2 units to 1 unit
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // simple tandem dup expansion from 1 unit to 2 units
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.TANDUP_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUP_SEQ_CIGARS, GATKSVVCFConstants.DUPLICATION_NUMBERS).stream())
                .sorted().collect(Collectors.toList()));

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.TANDUP_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUP_SEQ_CIGARS, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.TANDUP_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup contraction from 3 units to 2 units
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup expansion from 2 units to 3 units
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFConstants.TANDUP_EXPANSION_STRING, GATKSVVCFConstants.DUP_REPEAT_UNIT_REF_SPAN, GATKSVVCFConstants.DUPLICATION_NUMBERS, GATKSVVCFConstants.DUP_ANNOTATIONS_IMPRECISE, GATKSVVCFConstants.HOMOLOGY, GATKSVVCFConstants.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));
    }
}