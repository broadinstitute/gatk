package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.TYPES.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.TYPES.DEL;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.TYPES.DUP;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.DUP_TAN_CONTRACTION_STRING;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.DUP_TAN_EXPANSION_STRING;

public class DiscoverVariantsFromContigAlignmentsSAMSparkUnitTest extends GATKBaseTest {

    //  for methods used in this CLI are used only in this class

    @Test(groups = "sv")
    public void testFilterByRegionTooSmall() {
        final byte[] contigSequence = SimpleSVDiscoveryTestDataProvider.LONG_CONTIG1.getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval(SimpleSVDiscoveryTestDataProvider.chrForLongContig1, 20138007, 20142231), 1, contigSequence.length - 1986, TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval(SimpleSVDiscoveryTestDataProvider.chrForLongContig1, 20152030, 20154634), 3604, contigSequence.length, TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertFalse( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region1, region2, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region2, region1, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( DiscoverVariantsFromContigAlignmentsSAMSpark.firstAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test(groups = "sv")
    public void testFilterByNextAlignmentMayBeInsertion() {
        final AlignmentInterval overlappingRegion1 = new AlignmentInterval(new SimpleInterval("19", 48699881, 48700034), 1, 154, TextCigarCodec.decode("47S154M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval overlappingRegion2 = new AlignmentInterval(new SimpleInterval("19", 48700584, 48700668), 117, 201, TextCigarCodec.decode("116H85M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertTrue(DiscoverVariantsFromContigAlignmentsSAMSpark.nextAlignmentMayBeInsertion(overlappingRegion1, overlappingRegion2,  CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, 50,true));
    }

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
