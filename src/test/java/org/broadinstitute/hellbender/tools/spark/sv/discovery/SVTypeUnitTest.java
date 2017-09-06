package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.List;

import static  org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.createBracketedSymbAlleleString;

public class SVTypeUnitTest {

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
    // Allele production
    // -----------------------------------------------------------------------------------------------
    private static void worker(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                               final String expectedSymbolicAltAlleleStringWithoutBracket,
                               final int expectedSvLen,
                               final String expectedFirstFieldInIdString) throws IOException {

        final SvType SvType = SvTypeInference.inferFromNovelAdjacency(novelAdjacencyReferenceLocations);
        final List<Allele> producedAlleles = AnnotatedVariantProducer.produceAlleles(novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc, SVDiscoveryTestDataProvider.reference, SvType);

        Assert.assertEquals(producedAlleles.size(), 2);
        Assert.assertTrue(producedAlleles.get(0).isReference() && producedAlleles.get(1).isNonReference());
        Assert.assertTrue(producedAlleles.get(1).isSymbolic() && producedAlleles.get(1).toString().equals(createBracketedSymbAlleleString(expectedSymbolicAltAlleleStringWithoutBracket)));

        Assert.assertEquals(SvType.getSVLength(), expectedSvLen);

        final String variantId = SvType.getInternalVariantId();
        Assert.assertFalse(variantId.isEmpty());
        final String[] fields = variantId.split(GATKSVVCFConstants.INTERVAL_VARIANT_ID_FIELD_SEPARATOR);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(fields[0], expectedFirstFieldInIdString);
        final String expectedSecondField = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig(),
                expectedThirdField  = String.valueOf(novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd()),
                expectedFourthField = String.valueOf(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart());
        Assert.assertEquals(fields[1], expectedSecondField);
        Assert.assertEquals(fields[2], expectedThirdField);
        Assert.assertEquals(fields[3], expectedFourthField);
    }

    @Test(groups = "sv")
    public void testAltAlleleSvLenAndIdProductions() throws IOException {

        // inversion
        worker(SVDiscoveryTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_INV_IN_HEADER, 14644,
                GATKSVVCFConstants.INV33);
        worker(SVDiscoveryTestDataProvider.forSimpleInversionWithHom_leftPlus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_INV_IN_HEADER, 405,
                GATKSVVCFConstants.INV55);

        // simple deletion
        worker(SVDiscoveryTestDataProvider.forSimpleDeletion_plus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL_IN_HEADER, -20,
                SimpleSVType.TYPES.DEL.name());

        // simple insertion
        worker(SVDiscoveryTestDataProvider.forSimpleInsertion_minus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_INS_IN_HEADER, 50,
                SimpleSVType.TYPES.INS.name());

        // long range substitution (i.e. scarred deletion)
        worker(SVDiscoveryTestDataProvider.forLongRangeSubstitution_plus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL_IN_HEADER, -20,
                SimpleSVType.TYPES.DEL.name());

        // simple deletion with homology
        worker(SVDiscoveryTestDataProvider.forDeletionWithHomology_minus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL_IN_HEADER, -38,
                SimpleSVType.TYPES.DEL.name());

        // simple tandem dup contraction from 2 units to 1 unit
        worker(SVDiscoveryTestDataProvider.forSimpleTanDupContraction_plus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL_IN_HEADER, -10,
                GATKSVVCFConstants.TANDUP_CONTRACTION_INTERNAL_ID_START_STRING);

        // simple tandem dup expansion from 1 unit to 2 units
        worker(SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_minus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP_IN_HEADER, 10,
                GATKSVVCFConstants.TANDUP_EXPANSION_INTERNAL_ID_START_STRING);

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        worker(SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP_IN_HEADER, 99,
                GATKSVVCFConstants.TANDUP_EXPANSION_INTERNAL_ID_START_STRING);

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        worker(SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP_IN_HEADER, 96,
                GATKSVVCFConstants.TANDUP_EXPANSION_INTERNAL_ID_START_STRING);

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        worker(SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL_IN_HEADER, -96,
                GATKSVVCFConstants.TANDUP_CONTRACTION_INTERNAL_ID_START_STRING);

        // tandem dup contraction from 3 units to 2 units
        worker(SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DEL_IN_HEADER, -96,
                GATKSVVCFConstants.TANDUP_CONTRACTION_INTERNAL_ID_START_STRING);

        // tandem dup expansion from 2 units to 3 units
        worker(SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3(), GATKSVVCFConstants.SYMB_ALT_ALLELE_DUP_IN_HEADER, 96,
                GATKSVVCFConstants.TANDUP_EXPANSION_INTERNAL_ID_START_STRING);
    }
}
