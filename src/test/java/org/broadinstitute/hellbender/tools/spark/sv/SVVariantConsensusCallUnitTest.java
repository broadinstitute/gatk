package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
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
import java.util.*;
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
    // Variant length
    // -----------------------------------------------------------------------------------------------
    private static void seeIfItWorks_SvLen(final NovelAdjacencyReferenceLocations breakpoints,
                                           final int expectedSvLength) {
        Assert.assertEquals(SVVariantConsensusCall.getSvLength(breakpoints), expectedSvLength);
    }

    @Test
    public void testInferSvLength() {
        // inversion
        NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3();
        seeIfItWorks_SvLen(breakpoints, 14644);

        breakpoints = SVCallerTestDataProvider.forSimpleInversionWithHom_leftPlus._3();
        seeIfItWorks_SvLen(breakpoints, 405);

        // simple deletion
        breakpoints = SVCallerTestDataProvider.forSimpleDeletion_plus._3();
        seeIfItWorks_SvLen(breakpoints, 20);

        // simple insertion
        breakpoints = SVCallerTestDataProvider.forSimpleInsertion_minus._3();
        seeIfItWorks_SvLen(breakpoints, 50);

        // long range substitution
        breakpoints = SVCallerTestDataProvider.forLongRangeSubstitution_plus._3();
        seeIfItWorks_SvLen(breakpoints, 20);

        // simple deletion with homology
        breakpoints = SVCallerTestDataProvider.forDeletionWithHomology_minus._3();
        seeIfItWorks_SvLen(breakpoints, 38);

        // simple tandem dup contraction from 2 units to 1 unit
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupContraction_plus._3();
        seeIfItWorks_SvLen(breakpoints, 10);

        // simple tandem dup expansion from 1 unit to 2 units
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansion_minus._3();
        seeIfItWorks_SvLen(breakpoints, 10);

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3();
        seeIfItWorks_SvLen(breakpoints, 99);

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3();
        seeIfItWorks_SvLen(breakpoints, 96);

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3();
        seeIfItWorks_SvLen(breakpoints, 96);

        // tandem dup contraction from 3 units to 2 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3();
        seeIfItWorks_SvLen(breakpoints, 96);

        // tandem dup expansion from 2 units to 3 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3();
        seeIfItWorks_SvLen(breakpoints, 96);
    }

    // -----------------------------------------------------------------------------------------------
    // Allele production
    // -----------------------------------------------------------------------------------------------
    private static void seeIfItWorks_alleleProduction(final List<Allele> producedAlleles,
                                                      final String expectedSymbolicAltAlleleString) {
        Assert.assertEquals(producedAlleles.size(), 2);
        Assert.assertTrue(producedAlleles.get(0).isReference() && producedAlleles.get(1).isNonReference());
        Assert.assertTrue(producedAlleles.get(1).isSymbolic() && producedAlleles.get(1).toString().equals(expectedSymbolicAltAlleleString));
    }

    @Test
    public void testProduceAlleles() throws IOException {

        // inversion
        List<Allele> alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.INV);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_INV);

        // simple deletion
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forSimpleDeletion_plus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.DEL);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DEL);

        // simple insertion
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forSimpleInsertion_minus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.INS);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_INS);

        // long range substitution (i.e. scarred deletion)
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forLongRangeSubstitution_plus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.DEL);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DEL);

        // simple deletion with homology
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forDeletionWithHomology_minus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.DEL);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DEL);

        // simple tandem dup contraction from 2 units to 1 unit
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forSimpleTanDupContraction_plus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.DEL);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DEL);

        // simple tandem dup expansion from 1 unit to 2 units
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forSimpleTanDupExpansion_minus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.INS);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_INS);

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.DUP);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DUP);

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.INS);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_INS);

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.DEL);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DEL);

        // tandem dup contraction from 3 units to 2 units
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.DEL);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_DEL);

        // tandem dup expansion from 2 units to 3 units
        alleles = SVVariantConsensusCall.produceAlleles(SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3(), SVCallerTestDataProvider.reference, GATKSVVCFHeaderLines.SVTYPES.INS);
        seeIfItWorks_alleleProduction(alleles, SVConstants.CallingStepConstants.VCF_ALT_ALLELE_STRING_INS);
    }

    // -----------------------------------------------------------------------------------------------
    // VariantID production
    // -----------------------------------------------------------------------------------------------
    private static void seeIfItWorks_idProduction(final String idCompactString,
                                                  final String expectedFirstField, final String expectedSecondField,
                                                  final String expectedThirdField, final String expectedFourthField){
        Assert.assertFalse(idCompactString.isEmpty());
        final String[] fields = idCompactString.split(SVConstants.CallingStepConstants.VARIANT_ID_FIELD_SEPARATOR);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(fields[0], expectedFirstField);
        Assert.assertEquals(fields[1], expectedSecondField);
        Assert.assertEquals(fields[2], expectedThirdField);
        Assert.assertEquals(fields[3], expectedFourthField);
    }

    @Test
    public void testProduceVariantId() {

        // inversion
        NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.INV), GATKSVVCFHeaderLines.INV_3_TO_5,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        breakpoints = SVCallerTestDataProvider.forSimpleInversionWithHom_leftMinus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.INV), GATKSVVCFHeaderLines.INV_5_TO_3,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // simple deletion
        breakpoints = SVCallerTestDataProvider.forSimpleDeletion_minus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.DEL), GATKSVVCFHeaderLines.SVTYPES.DEL.name(),
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // simple insertion
        breakpoints = SVCallerTestDataProvider.forSimpleInsertion_plus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.INS), GATKSVVCFHeaderLines.SVTYPES.INS.name(),
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // long range substitution
        breakpoints = SVCallerTestDataProvider.forLongRangeSubstitution_minus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.DEL), GATKSVVCFHeaderLines.SVTYPES.DEL.name(),
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // simple deletion with homology
        breakpoints = SVCallerTestDataProvider.forDeletionWithHomology_plus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.DEL), GATKSVVCFHeaderLines.SVTYPES.DEL.name(),
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // simple tandem dup contraction from 2 units to 1 unit
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupContraction_minus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.DEL), GATKSVVCFHeaderLines.SVTYPES.DEL.name()+"_"+GATKSVVCFHeaderLines.SVTYPES.DUP.name()+"_"+ SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // simple tandem dup expansion from 1 unit to 2 units
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansion_plus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.DUP), GATKSVVCFHeaderLines.SVTYPES.DUP.name()+ SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.INS), GATKSVVCFHeaderLines.SVTYPES.INS.name()+"_"+GATKSVVCFHeaderLines.SVTYPES.DUP.name()+"_"+ SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.INS), GATKSVVCFHeaderLines.SVTYPES.INS.name()+"_"+GATKSVVCFHeaderLines.SVTYPES.DUP.name()+"_"+ SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.DEL), GATKSVVCFHeaderLines.SVTYPES.DEL.name()+"_"+GATKSVVCFHeaderLines.SVTYPES.DUP.name()+"_"+ SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // tandem dup contraction from 3 units to 2 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.DEL), GATKSVVCFHeaderLines.SVTYPES.DEL.name()+"_"+GATKSVVCFHeaderLines.SVTYPES.DUP.name()+"_"+ SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));

        // tandem dup expansion from 2 units to 3 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus._3();
        seeIfItWorks_idProduction(SVVariantConsensusCall.produceVariantId(breakpoints, GATKSVVCFHeaderLines.SVTYPES.INS), GATKSVVCFHeaderLines.SVTYPES.INS.name()+"_"+GATKSVVCFHeaderLines.SVTYPES.DUP.name()+"_"+ SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING,
                breakpoints.leftJustifiedLeftRefLoc.getContig(), String.valueOf(breakpoints.leftJustifiedLeftRefLoc.getEnd()), String.valueOf(breakpoints.leftJustifiedRightRefLoc.getStart()));
    }

    // -----------------------------------------------------------------------------------------------
    // Type inference
    // -----------------------------------------------------------------------------------------------
    private static void seeIfItWorks_typeInference(final NovelAdjacencyReferenceLocations breakpoints,
                                                   final Map<String, String> expectedAttributesAsStrings) throws IOException{

        final Map<String, Object> attributeMap = SVVariantConsensusCall.getType(breakpoints);
        Assert.assertEquals(attributeMap.size(), expectedAttributesAsStrings.size());
        attributeMap.forEach((k,v) -> Assert.assertEquals(expectedAttributesAsStrings.get(k), (String)v));
    }

    @Test
    public void testGetType() throws IOException {

        // inversion
        NovelAdjacencyReferenceLocations breakpoints = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INV.name(), GATKSVVCFHeaderLines.INV_3_TO_5, ""));

        breakpoints = SVCallerTestDataProvider.forSimpleInversionWithHom_leftPlus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INV.name(), GATKSVVCFHeaderLines.INV_5_TO_3, ""));

        // simple deletion
        breakpoints = SVCallerTestDataProvider.forSimpleDeletion_plus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DEL.name()));

        // simple insertion
        breakpoints = SVCallerTestDataProvider.forSimpleInsertion_minus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.INS.name()));

        // long range substitution
        breakpoints = SVCallerTestDataProvider.forLongRangeSubstitution_plus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DEL.name()));

        // simple deletion with homology
        breakpoints = SVCallerTestDataProvider.forDeletionWithHomology_minus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DEL.name()));

        // simple tandem dup contraction from 2 units to 1 unit
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupContraction_plus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DEL.name(), SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, ""));

        // simple tandem dup expansion from 1 unit to 2 units
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansion_minus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DUP.name(), SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, ""));

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        breakpoints = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DUP.name(), SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, ""));

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DUP.name(), SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, ""));

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DEL.name(), SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, ""));

        // tandem dup contraction from 3 units to 2 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DEL.name(), SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, ""));

        // tandem dup expansion from 2 units to 3 units
        breakpoints = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3();
        seeIfItWorks_typeInference(breakpoints, ImmutableMap.of(GATKSVVCFHeaderLines.SVTYPE, GATKSVVCFHeaderLines.SVTYPES.DUP.name(), SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, ""));
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

        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFHeaderLines.MAPPING_QUALITIES)).split(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR),
                            expectedMappingQualitiesAsStrings);
        Assert.assertEquals(((String)attributeMap.get(GATKSVVCFHeaderLines.ALIGN_LENGTHS)).split(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR),
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
                                                 final Set<String> expectedAttributeKeys) throws IOException {

        final AlignmentRegion region1 = testData._1();
        final AlignmentRegion region2 = testData._2();
        final byte[] contigSeq = null; // hack, as the contig sequence is really not necessary for this test purpose

        final Iterable<ChimericAlignment> evidence = Collections.singletonList(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));

        final NovelAdjacencyReferenceLocations breakpoints = testData._3();

        final VariantContext variantContext
                = SVVariantConsensusCall.callVariantsFromConsensus(new Tuple2<>(breakpoints, evidence), SparkContextFactory.getTestSparkContext().broadcast(SVCallerTestDataProvider.reference));

        final Set<String> attributeKeys = variantContext.getAttributes().keySet();

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
                Sets.newHashSet(GATKSVVCFHeaderLines.INV_3_TO_5, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream()).collect(Collectors.toSet()));

        // simple deletion
        testData = SVCallerTestDataProvider.forSimpleDeletion_minus;

        seeIfItWorks_integrative(testData, commonAttributes);

        // simple insertion
        testData = SVCallerTestDataProvider.forSimpleInsertion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFHeaderLines.INSERTED_SEQUENCE).stream())
                .collect(Collectors.toSet()));

        // long range substitution
        testData = SVCallerTestDataProvider.forLongRangeSubstitution_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFHeaderLines.INSERTED_SEQUENCE).stream())
                .collect(Collectors.toSet()));

        // simple deletion with homology
        testData = SVCallerTestDataProvider.forDeletionWithHomology_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .collect(Collectors.toSet()));

        // simple tandem dup contraction from 2 units to 1 unit
        testData = SVCallerTestDataProvider.forSimpleTanDupContraction_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .collect(Collectors.toSet()));

        // simple tandem dup expansion from 1 unit to 2 units
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansion_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS).stream())
                .collect(Collectors.toSet()));

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.INSERTED_SEQUENCE).stream())
                .collect(Collectors.toSet()));

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .collect(Collectors.toSet()));

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .collect(Collectors.toSet()));

        // tandem dup contraction from 3 units to 2 units
        testData = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_CONTRACTION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .collect(Collectors.toSet()));

        // tandem dup expansion from 2 units to 3 units
        testData = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;

        seeIfItWorks_integrative(testData, Stream.concat( commonAttributes.stream(),
                Sets.newHashSet(SVConstants.CallingStepConstants.TANDUP_EXPANSION_STRING, GATKSVVCFHeaderLines.DUPLICATED_SEQUENCE, GATKSVVCFHeaderLines.DUPLICATION_NUMBERS, GATKSVVCFHeaderLines.HOMOLOGY, GATKSVVCFHeaderLines.HOMOLOGY_LENGTH).stream())
                .collect(Collectors.toSet()));
    }
}