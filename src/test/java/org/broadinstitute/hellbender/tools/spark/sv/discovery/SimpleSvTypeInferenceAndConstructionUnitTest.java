package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.ImmutableSet;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SVDiscoveryTestDataProvider.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.TYPES.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.createBracketedSymbAlleleString;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SvTypeInference.inferFromNovelAdjacency;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;


public class SimpleSvTypeInferenceAndConstructionUnitTest extends BaseTest {

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SVDiscoveryTestDataProvider.testDataInitialized) {
            new SVDiscoveryTestDataProvider();
        }
    }

    @Test(groups = "sv", dataProvider = "forTypeInference")
    public void testGetType(final NovelAdjacencyReferenceLocations breakpoints, final String expectedTypeString,
                            final Set<String> expectedAttributeIDs) {

        final SvType variant = SvTypeInference.inferFromNovelAdjacency(breakpoints);
        Assert.assertEquals(variant.toString(), expectedTypeString);

        final Set<String> attributeIDs = variant.getTypeSpecificAttributes().keySet();
        Assert.assertEquals(attributeIDs, expectedAttributeIDs);
    }

    @DataProvider(name = "forTypeInference")
    private Object[][] forTypeInference() {
        final List<Object[]> data = new ArrayList<>(20);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3(), INV.name(), ImmutableSet.of(INV33)});

        data.add(new Object[]{forSimpleInversionWithHom_leftPlus._3(), INV.name(), ImmutableSet.of(INV55)});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_plus._3(), DEL.name(), Collections.emptySet()});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_minus._3(), INS.name(), Collections.emptySet()});

        // long range substitution
        data.add(new Object[]{forLongRangeSubstitution_plus._3(), DEL.name(), Collections.emptySet()});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_minus._3(), DEL.name(), Collections.emptySet()});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_plus._3(), DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_minus._3(), DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_plus._3(), DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus._3(), DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus._3(), DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus._3(), DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus._3(), DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(groups = "sv", dataProvider = "forAltAlleleSvLenAndIdProductions")
    public static void testAltAlleleSvLenAndIdProductions(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                                                          final SvType simpleType,
                                                          final String expectedSymbolicAltAlleleStringWithoutBracket,
                                                          final int expectedSvLen,
                                                          final String expectedTypeInfoInIdString) throws IOException {

        final List<Allele> producedAlleles = AnnotatedVariantProducer.produceAlleles(novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc, reference, simpleType);

        Assert.assertEquals(producedAlleles.size(), 2);
        Assert.assertTrue(producedAlleles.get(0).isReference() && producedAlleles.get(1).isNonReference() && producedAlleles.get(1).isSymbolic());
        Assert.assertEquals(producedAlleles.get(1).toString(), createBracketedSymbAlleleString(expectedSymbolicAltAlleleStringWithoutBracket));

        Assert.assertEquals(simpleType.getSVLength(), expectedSvLen);

        final String variantId = simpleType.getInternalVariantId();
        Assert.assertFalse(variantId.isEmpty());
        final String[] fields = variantId.split(INTERVAL_VARIANT_ID_FIELD_SEPARATOR);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(fields[0], expectedTypeInfoInIdString);
        final String expectedRefContigNameInIdString = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig(),
                     expectedPOSInfoInIdString = String.valueOf(novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd()),
                     expectedENDInfoInIdString = String.valueOf(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart());
        Assert.assertEquals(fields[1], expectedRefContigNameInIdString);
        Assert.assertEquals(fields[2], expectedPOSInfoInIdString);
        Assert.assertEquals(fields[3], expectedENDInfoInIdString);
    }

    @DataProvider(name = "forAltAlleleSvLenAndIdProductions")
    private Object[][] forAltAlleleSvLenAndIdProductions(){
        final List<Object[]> data = new ArrayList<>(20);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3(),
                inferFromNovelAdjacency(forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3()),
                SYMB_ALT_ALLELE_INV, 14644, INV33});
        data.add(new Object[]{forSimpleInversionWithHom_leftPlus._3(),
                inferFromNovelAdjacency(forSimpleInversionWithHom_leftPlus._3()),
                SYMB_ALT_ALLELE_INV, 405, INV55});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_plus._3(), inferFromNovelAdjacency(forSimpleDeletion_plus._3()),
                SYMB_ALT_ALLELE_DEL, -20, SimpleSVType.TYPES.DEL.name()});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_minus._3(), inferFromNovelAdjacency(forSimpleInsertion_minus._3()),
                SYMB_ALT_ALLELE_INS, 50, SimpleSVType.TYPES.INS.name()});

        // long range substitution (i.e. scarred deletion)
        data.add(new Object[]{forLongRangeSubstitution_plus._3(), inferFromNovelAdjacency(forLongRangeSubstitution_plus._3()),
                SYMB_ALT_ALLELE_DEL, -20, SimpleSVType.TYPES.DEL.name()});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_minus._3(), inferFromNovelAdjacency(forDeletionWithHomology_minus._3()),
                SYMB_ALT_ALLELE_DEL, -38, SimpleSVType.TYPES.DEL.name()});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_plus._3(), inferFromNovelAdjacency(forSimpleTanDupContraction_plus._3()),
                SYMB_ALT_ALLELE_DEL, -10,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_minus._3(), inferFromNovelAdjacency(forSimpleTanDupExpansion_minus._3()),
                SYMB_ALT_ALLELE_DUP, 10,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_plus._3(), inferFromNovelAdjacency(forSimpleTanDupExpansionWithNovelIns_plus._3()),
                SYMB_ALT_ALLELE_DUP, 99,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus._3(), inferFromNovelAdjacency(forComplexTanDup_1to2_pseudoHom_minus._3()),
                SYMB_ALT_ALLELE_DUP, 96,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus._3(), inferFromNovelAdjacency(forComplexTanDup_2to1_pseudoHom_plus._3()),
                SYMB_ALT_ALLELE_DEL, -96,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});


        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus._3(), inferFromNovelAdjacency(forComplexTanDup_3to2_noPseudoHom_minus._3()),
                SYMB_ALT_ALLELE_DEL, -96,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus._3(), inferFromNovelAdjacency(forComplexTanDup_2to3_noPseudoHom_plus._3()),
                SYMB_ALT_ALLELE_DUP, 96,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        return data.toArray(new Object[data.size()][]);
    }
}
