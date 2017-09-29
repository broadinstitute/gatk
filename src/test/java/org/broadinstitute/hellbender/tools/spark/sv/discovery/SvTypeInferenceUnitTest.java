package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.ImmutableSet;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;
import java.util.Set;

public class SvTypeInferenceUnitTest extends BaseTest {

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SVDiscoveryTestDataProvider.testDataInitialized) {
            new SVDiscoveryTestDataProvider();
        }
    }


    private static void seeIfItWorks_typeInference(final NovelAdjacencyReferenceLocations breakpoints,
                                                   final String expectedTypeString,
                                                   final Set<String> expectedFlags) throws IOException {

        final SvType variant = SvTypeInference.inferFromNovelAdjacency(breakpoints);
        Assert.assertEquals(variant.toString(), expectedTypeString);

        final Set<String> flags = variant.getTypeSpecificAttributes().keySet();
        Assert.assertEquals(flags.size(), expectedFlags.size());
        if (!expectedFlags.isEmpty()) Assert.assertTrue(flags.containsAll(expectedFlags));
    }

    @Test(groups = "sv")
    public void testGetType() throws IOException {

        // inversion
        NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.INV.name(), ImmutableSet.of(GATKSVVCFConstants.INV33));

        breakpoints = SVDiscoveryTestDataProvider.forSimpleInversionWithHom_leftPlus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.INV.name(), ImmutableSet.of(GATKSVVCFConstants.INV55));

        // simple deletion
        breakpoints = SVDiscoveryTestDataProvider.forSimpleDeletion_plus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DEL.name(), Collections.emptySet());

        // simple insertion
        breakpoints = SVDiscoveryTestDataProvider.forSimpleInsertion_minus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.INS.name(), Collections.emptySet());

        // long range substitution
        breakpoints = SVDiscoveryTestDataProvider.forLongRangeSubstitution_plus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DEL.name(), Collections.emptySet());

        // simple deletion with homology
        breakpoints = SVDiscoveryTestDataProvider.forDeletionWithHomology_minus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DEL.name(), Collections.emptySet());

        // simple tandem dup contraction from 2 units to 1 unit
        breakpoints = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_plus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DEL.name(), ImmutableSet.of(GATKSVVCFConstants.TANDUP_CONTRACTION_STRING));

        // simple tandem dup expansion from 1 unit to 2 units
        breakpoints = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_minus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DUP.name(), ImmutableSet.of(GATKSVVCFConstants.TANDUP_EXPANSION_STRING));

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        breakpoints = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DUP.name(), ImmutableSet.of(GATKSVVCFConstants.TANDUP_EXPANSION_STRING));

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DUP.name(), ImmutableSet.of(GATKSVVCFConstants.TANDUP_EXPANSION_STRING));

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DEL.name(), ImmutableSet.of(GATKSVVCFConstants.TANDUP_CONTRACTION_STRING));

        // tandem dup contraction from 3 units to 2 units
        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DEL.name(), ImmutableSet.of(GATKSVVCFConstants.TANDUP_CONTRACTION_STRING));

        // tandem dup expansion from 2 units to 3 units
        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3();
        seeIfItWorks_typeInference(breakpoints, SimpleSVType.TYPES.DUP.name(), ImmutableSet.of(GATKSVVCFConstants.TANDUP_EXPANSION_STRING));
    }
}
