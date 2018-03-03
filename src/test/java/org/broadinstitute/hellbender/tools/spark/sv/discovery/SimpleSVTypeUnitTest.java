package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark.inferSimpleTypeFromNovelAdjacency;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.createBracketedSymbAlleleString;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public class SimpleSVTypeUnitTest extends GATKBaseTest {

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
    // Allele production
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv", dataProvider = "forAltAlleleSvLenAndIdProductions")
    public static void testAltAlleleSvLenAndIdProductions(final NovelAdjacencyAndAltHaplotype novelAdjacencyReferenceLocations,
                                                          final SvType simpleType,
                                                          final String expectedSymbolicAltAlleleStringWithoutBracket,
                                                          final int expectedSvLen,
                                                          final String expectedTypeInfoInIdString) throws IOException {

        final List<Allele> producedAlleles = AnnotatedVariantProducer.produceAlleles(novelAdjacencyReferenceLocations.getLeftJustifiedLeftRefLoc(), b37_reference, simpleType);

        Assert.assertEquals(producedAlleles.size(), 2);
        Assert.assertTrue(producedAlleles.get(0).isReference() && producedAlleles.get(1).isNonReference() && producedAlleles.get(1).isSymbolic());
        Assert.assertEquals(producedAlleles.get(1).toString(), createBracketedSymbAlleleString(expectedSymbolicAltAlleleStringWithoutBracket));

        Assert.assertEquals(simpleType.getSVLength(), expectedSvLen);

        final String variantId = simpleType.getInternalVariantId();
        Assert.assertFalse(variantId.isEmpty());
        final String[] fields = variantId.split(INTERVAL_VARIANT_ID_FIELD_SEPARATOR);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(fields[0], expectedTypeInfoInIdString);
        final String expectedRefContigNameInIdString = novelAdjacencyReferenceLocations.getLeftJustifiedLeftRefLoc().getContig(),
                expectedPOSInfoInIdString = String.valueOf(novelAdjacencyReferenceLocations.getLeftJustifiedLeftRefLoc().getEnd()),
                expectedENDInfoInIdString = String.valueOf(novelAdjacencyReferenceLocations.getLeftJustifiedRightRefLoc().getStart());
        Assert.assertEquals(fields[1], expectedRefContigNameInIdString);
        Assert.assertEquals(fields[2], expectedPOSInfoInIdString);
        Assert.assertEquals(fields[3], expectedENDInfoInIdString);
    }

    @DataProvider(name = "forAltAlleleSvLenAndIdProductions")
    private Object[][] forAltAlleleSvLenAndIdProductions() {
        final List<Object[]> data = new ArrayList<>(20);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint.biPathBubble,
                inferSimpleTypeFromNovelAdjacency(forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint.biPathBubble),
                SYMB_ALT_ALLELE_INV, 14644, INV33});
        data.add(new Object[]{forSimpleInversionWithHom_leftPlus.biPathBubble,
                inferSimpleTypeFromNovelAdjacency(forSimpleInversionWithHom_leftPlus.biPathBubble),
                SYMB_ALT_ALLELE_INV, 405, INV55});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_plus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleDeletion_plus.biPathBubble),
                SYMB_ALT_ALLELE_DEL, -20, SimpleSVType.TYPES.DEL.name()});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_minus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleInsertion_minus.biPathBubble),
                SYMB_ALT_ALLELE_INS, 50, SimpleSVType.TYPES.INS.name()});

        // long range substitution (i.e. scarred deletion)
        data.add(new Object[]{forLongRangeSubstitution_plus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forLongRangeSubstitution_plus.biPathBubble),
                SYMB_ALT_ALLELE_DEL, -20, SimpleSVType.TYPES.DEL.name()});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_minus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forDeletionWithHomology_minus.biPathBubble),
                SYMB_ALT_ALLELE_DEL, -38, SimpleSVType.TYPES.DEL.name()});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_plus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleTanDupContraction_plus.biPathBubble),
                SYMB_ALT_ALLELE_DEL, -10,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_minus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleTanDupExpansion_minus.biPathBubble),
                SYMB_ALT_ALLELE_DUP, 10,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_plus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forSimpleTanDupExpansionWithNovelIns_plus.biPathBubble),
                SYMB_ALT_ALLELE_DUP, 99,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_1to2_pseudoHom_minus.biPathBubble),
                SYMB_ALT_ALLELE_DUP, 96,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_2to1_pseudoHom_plus.biPathBubble),
                SYMB_ALT_ALLELE_DEL, -96,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});


        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_3to2_noPseudoHom_minus.biPathBubble),
                SYMB_ALT_ALLELE_DEL, -96,
                DUP_TAN_CONTRACTION_INTERNAL_ID_START_STRING});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus.biPathBubble, inferSimpleTypeFromNovelAdjacency(forComplexTanDup_2to3_noPseudoHom_plus.biPathBubble),
                SYMB_ALT_ALLELE_DUP, 96,
                DUP_TAN_EXPANSION_INTERNAL_ID_START_STRING});

        return data.toArray(new Object[data.size()][]);
    }
}
