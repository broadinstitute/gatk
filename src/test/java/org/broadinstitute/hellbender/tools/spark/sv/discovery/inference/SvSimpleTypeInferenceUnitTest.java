package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.collect.ImmutableSet;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.TYPES.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public class SvSimpleTypeInferenceUnitTest extends GATKBaseTest {

    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SimpleSVDiscoveryTestDataProvider.testDataInitialized) {
            new SimpleSVDiscoveryTestDataProvider();
        }
    }


    @Test(groups = "sv", dataProvider = "forTypeInference")
    public void testGetType(final NovelAdjacencyAndAltHaplotype breakpoints, final String expectedTypeString,
                            final Set<String> expectedAttributeIDs) {

        final SvType variant = DiscoverVariantsFromContigAlignmentsSAMSpark.inferSimpleTypeFromNovelAdjacency(breakpoints);
        Assert.assertEquals(variant.toString(), expectedTypeString);

        final Set<String> attributeIDs = variant.getTypeSpecificAttributes().keySet();
        Assert.assertEquals(attributeIDs, expectedAttributeIDs);
    }

    @DataProvider(name = "forTypeInference")
    private Object[][] forTypeInference() {
        final List<Object[]> data = new ArrayList<>(20);

        // inversion
        data.add(new Object[]{forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint.biPathBubble, INV.name(), ImmutableSet.of(INV33)});

        data.add(new Object[]{forSimpleInversionWithHom_leftPlus.biPathBubble, INV.name(), ImmutableSet.of(INV55)});

        // simple deletion
        data.add(new Object[]{forSimpleDeletion_plus.biPathBubble, DEL.name(), Collections.emptySet()});

        // simple insertion
        data.add(new Object[]{forSimpleInsertion_minus.biPathBubble, INS.name(), Collections.emptySet()});

        // long range substitution
        data.add(new Object[]{forLongRangeSubstitution_plus.biPathBubble, DEL.name(), Collections.emptySet()});

        // simple deletion with homology
        data.add(new Object[]{forDeletionWithHomology_minus.biPathBubble, DEL.name(), Collections.emptySet()});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{forSimpleTanDupContraction_plus.biPathBubble, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{forSimpleTanDupExpansion_minus.biPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{forSimpleTanDupExpansionWithNovelIns_plus.biPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{forComplexTanDup_1to2_pseudoHom_minus.biPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{forComplexTanDup_2to1_pseudoHom_plus.biPathBubble, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{forComplexTanDup_3to2_noPseudoHom_minus.biPathBubble, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING)});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{forComplexTanDup_2to3_noPseudoHom_plus.biPathBubble, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING)});

        return data.toArray(new Object[data.size()][]);
    }
}
