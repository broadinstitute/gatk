package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple2;
import scala.Tuple3;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class SVVariantConsensusCallUnitTest extends BaseTest {


    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SVCallerTestDataProvider.testDataInitialized) {
            new SVCallerTestDataProvider();
        }
    }


    // -----------------------------------------------------------------------------------------------
    // Type inference
    // -----------------------------------------------------------------------------------------------
    private static void seeIfItWorks_typeInference(final NovelAdjacencyReferenceLocations breakpoints,
                                                   final String expectedTypeString,
                                                   final Set<String> expectedFlags) throws IOException {

        final SvType variant = SVVariantConsensusCall.getType(breakpoints);
        Assert.assertEquals(variant.toString(), expectedTypeString);

        final Set<String> flags = variant.getTypeSpecificAttributes().keySet();
        Assert.assertEquals(flags.size(), expectedFlags.size());
        if (!expectedFlags.isEmpty()) Assert.assertTrue(flags.containsAll(expectedFlags));
    }

    @Test
    public void testGetType() throws IOException {

        // inversion
        NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.INV.name(), ImmutableSet.of(GATKSVVCFHeaderLines.INV33));

        breakpoints = SVCallerTestDataProvider.forSimpleInversionWithHom_leftPlus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.INV.name(), ImmutableSet.of(GATKSVVCFHeaderLines.INV55));

        // simple deletion
        breakpoints = SVCallerTestDataProvider.forSimpleDeletion_plus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DEL.name(), Collections.emptySet());

        // simple insertion
        breakpoints = SVCallerTestDataProvider.forSimpleInsertion_minus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.INS.name(), Collections.emptySet());

        // long range substitution
        breakpoints = SVCallerTestDataProvider.forLongRangeSubstitution_plus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DEL.name(), Collections.emptySet());

        // simple deletion with homology
        breakpoints = SVCallerTestDataProvider.forDeletionWithHomology_minus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DEL.name(), Collections.emptySet());

        // simple tandem dup contraction from 2 units to 1 unit
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupContraction_plus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DEL.name(), ImmutableSet.of(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING));

        // simple tandem dup expansion from 1 unit to 2 units
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansion_minus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DUP.name(), ImmutableSet.of(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING));

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DUP.name(), ImmutableSet.of(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING));

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DUP.name(), ImmutableSet.of(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING));

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DEL.name(), ImmutableSet.of(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING));

        // tandem dup contraction from 3 units to 2 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DEL.name(), ImmutableSet.of(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING));

        // tandem dup expansion from 2 units to 3 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3();
        seeIfItWorks_typeInference(breakpoints, SvType.TYPES.DUP.name(), ImmutableSet.of(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING));
    }

    // -----------------------------------------------------------------------------------------------
    // Evidence summary annotation
    // -----------------------------------------------------------------------------------------------
    /**
     * Not an exhaustive test on all attributes, only tests:
     * MAPPING_QUALITIES, ALIGNMENT_LENGTH
     */
    private static void seeIfItWorks_evidenceAnnotation(final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData,
                                                        final String[] expectedMappingQualitiesAsStrings,
                                                        final String[] expectedAlignmentLengthsAsStrings) throws IOException {

        final AlignmentRegion region1 = testData._1();
        final AlignmentRegion region2 = testData._2();
        final byte[] contigSeq = null; // hack, as the contig sequence is really not necessary for this test purpose

        final Map<String, Object> attributeMap =
                SVVariantConsensusCall.getEvidenceRelatedAnnotations(Collections.singletonList(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList())));

        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFHeaderLines.MAPPING_QUALITIES)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                            expectedMappingQualitiesAsStrings);
        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFHeaderLines.ALIGN_LENGTHS)).split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR),
                            expectedAlignmentLengthsAsStrings);
    }

    @Test
    public void testGetEvidenceRelatedAnnotations() throws IOException {

        // inversion
        Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(1984)});

        // simple deletion
        testData = SVCallerTestDataProvider.forSimpleDeletion_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple insertion
        testData = SVCallerTestDataProvider.forSimpleInsertion_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(100)});

        // long range substitution
        testData = SVCallerTestDataProvider.forLongRangeSubstitution_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple deletion with homology
        testData = SVCallerTestDataProvider.forDeletionWithHomology_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple tandem dup contraction from 2 units to 1 unit
        testData = SVCallerTestDataProvider.forSimpleTanDupContraction_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(40)});

        // simple tandem dup expansion from 1 unit to 2 units
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansion_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(50)});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(137)});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(127)});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(31)});

        // tandem dup contraction from 3 units to 2 units
        testData = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(31)});

        // tandem dup expansion from 2 units to 3 units
        testData = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;

        seeIfItWorks_evidenceAnnotation(testData, new String[]{"60"}, new String[]{String.valueOf(127)});
    }

    // -----------------------------------------------------------------------------------------------
    // Integrative test
    // -----------------------------------------------------------------------------------------------
    private static void seeIfItWorks_integrative(final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData,
                                                 final List<String> expectedAttributeKeys) throws IOException {

        final AlignmentRegion region1 = testData._1();
        final AlignmentRegion region2 = testData._2();
        final byte[] contigSeq = null; // hack, as the contig sequence is really not necessary for this test purpose

        final Iterable<ChimericAlignment> evidence = Collections.singletonList(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));

        final NovelAdjacencyReferenceLocations breakpoints = testData._3();

        final VariantContext variantContext
                = SVVariantConsensusCall.callVariantsFromConsensus(new Tuple2<>(breakpoints, evidence), SparkContextFactory.getTestSparkContext().broadcast(SVCallerTestDataProvider.reference));

        final List<String> attributeKeys = variantContext.getAttributes().keySet().stream().sorted().collect(Collectors.toList());

        Assert.assertEquals(attributeKeys, expectedAttributeKeys);
    }

    @Test
    public void testIntegrative() throws IOException {

        final Set<String> commonAttributes = Sets.newHashSet(VCFConstants.END_KEY, GATKSVVCFHeaderLines.SVLEN, GATKSVVCFHeaderLines.SVTYPE,
                GATKSVVCFHeaderLines.TOTAL_MAPPINGS, GATKSVVCFHeaderLines.HQ_MAPPINGS, GATKSVVCFHeaderLines.MAPPING_QUALITIES,
                GATKSVVCFHeaderLines.ALIGN_LENGTHS, GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, GATKSVVCFHeaderLines.ASSEMBLY_IDS, GATKSVVCFHeaderLines.CONTIG_IDS);

        // inversion
        Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFHeaderLines.INV33, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream()).sorted().collect(Collectors.toList()));

        // simple deletion
        testData = SVCallerTestDataProvider.forSimpleDeletion_minus;

        seeIfItWorks_integrative(testData, commonAttributes.stream().sorted().collect(Collectors.toList()));

        // simple insertion
        testData = SVCallerTestDataProvider.forSimpleInsertion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFHeaderLines.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()));

        // long range substitution
        testData = SVCallerTestDataProvider.forLongRangeSubstitution_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFHeaderLines.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()));

        // simple deletion with homology
        testData = SVCallerTestDataProvider.forDeletionWithHomology_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // simple tandem dup contraction from 2 units to 1 unit
        testData = SVCallerTestDataProvider.forSimpleTanDupContraction_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // simple tandem dup expansion from 1 unit to 2 units
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS).stream())
                .sorted().collect(Collectors.toList()));

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.INSERTED_SEQUENCE).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup contraction from 3 units to 2 units
        testData = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));

        // tandem dup expansion from 2 units to 3 units
        testData = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .sorted().collect(Collectors.toList()));
    }
}