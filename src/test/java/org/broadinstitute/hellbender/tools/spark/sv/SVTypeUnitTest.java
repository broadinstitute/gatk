package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.variant.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.List;


public class SVTypeUnitTest {

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
    // Allele production
    // -----------------------------------------------------------------------------------------------
    private static void worker(final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations,
                               final String expectedSymbolicAltAlleleString,
                               final int expectedSvLen,
                               final String expectedFirstFieldInIdString) throws IOException {

        final SvType SvType = SVVariantConsensusDiscovery.getType(novelAdjacencyReferenceLocations);
        final List<Allele> producedAlleles = SVVariantConsensusDiscovery.produceAlleles(novelAdjacencyReferenceLocations, SVCallerTestDataProvider.reference, SvType);

        Assert.assertEquals(producedAlleles.size(), 2);
        Assert.assertTrue(producedAlleles.get(0).isReference() && producedAlleles.get(1).isNonReference());
        Assert.assertTrue(producedAlleles.get(1).isSymbolic() && producedAlleles.get(1).toString().equals(expectedSymbolicAltAlleleString));

        Assert.assertEquals(SvType.getSVLength(), expectedSvLen);

        final String variantId = SvType.getVariantId();
        Assert.assertFalse(variantId.isEmpty());
        final String[] fields = variantId.split(SVConstants.DiscoveryStepConstants.VARIANT_ID_FIELD_SEPARATOR);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(fields[0], expectedFirstFieldInIdString);
        final String expectedSecondField = novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getContig(),
                     expectedThirdField  = String.valueOf(novelAdjacencyReferenceLocations.leftJustifiedLeftRefLoc.getEnd()),
                     expectedFourthField = String.valueOf(novelAdjacencyReferenceLocations.leftJustifiedRightRefLoc.getStart());
        Assert.assertEquals(fields[1], expectedSecondField);
        Assert.assertEquals(fields[2], expectedThirdField);
        Assert.assertEquals(fields[3], expectedFourthField);
    }

    @Test
    public void testAltAlleleSvLenAndIdProductions() throws IOException {

        // inversion
        worker(SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_INV, 14644,
                GATKSVVCFHeaderLines.INV33);
        worker(SVCallerTestDataProvider.forSimpleInversionWithHom_leftPlus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_INV, 405,
                GATKSVVCFHeaderLines.INV55);

        // simple deletion
        worker(SVCallerTestDataProvider.forSimpleDeletion_plus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DEL, -20,
                SvType.TYPES.DEL.name());

        // simple insertion
        worker(SVCallerTestDataProvider.forSimpleInsertion_minus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_INS, 50,
                SvType.TYPES.INS.name());

        // long range substitution (i.e. scarred deletion)
        worker(SVCallerTestDataProvider.forLongRangeSubstitution_plus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DEL, -20,
                SvType.TYPES.DEL.name());

        // simple deletion with homology
        worker(SVCallerTestDataProvider.forDeletionWithHomology_minus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DEL, -38,
                SvType.TYPES.DEL.name());

        // simple tandem dup contraction from 2 units to 1 unit
        worker(SVCallerTestDataProvider.forSimpleTanDupContraction_plus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DEL, -10,
                SvType.TANDEMUPLICATION_CONTRACTION_ID_START_STRING);

        // simple tandem dup expansion from 1 unit to 2 units
        worker(SVCallerTestDataProvider.forSimpleTanDupExpansion_minus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DUP, 10,
                SvType.TANDEMUPLICATION_EXPANSION_ID_START_STRING);

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        worker(SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DUP, 99,
                SvType.TANDEMUPLICATION_EXPANSION_ID_START_STRING);

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        worker(SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DUP, 96,
                SvType.TANDEMUPLICATION_EXPANSION_ID_START_STRING);

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        worker(SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DEL, -96,
                SvType.TANDEMUPLICATION_CONTRACTION_ID_START_STRING);

        // tandem dup contraction from 3 units to 2 units
        worker(SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DEL, -96,
                SvType.TANDEMUPLICATION_CONTRACTION_ID_START_STRING);

        // tandem dup expansion from 2 units to 3 units
        worker(SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3(), SVConstants.DiscoveryStepConstants.VCF_ALT_ALLELE_STRING_DUP, 96,
                SvType.TANDEMUPLICATION_EXPANSION_ID_START_STRING);
    }
}
