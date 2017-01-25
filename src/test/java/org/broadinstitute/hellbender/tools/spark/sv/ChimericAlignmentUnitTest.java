package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple3;

public class ChimericAlignmentUnitTest extends BaseTest {

    @Test
    public void testFilterByRegionTooSmall() {
        final byte[] contigSequence = SVCallerTestDataProvider.LONG_CONTIG1.getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("702700", "702700", new SimpleInterval(SVCallerTestDataProvider.chrForLongContig1, 20138007, 20142231), TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 1, contigSequence.length - 1986);
        final AlignmentRegion region2 = new AlignmentRegion("702700", "702700", new SimpleInterval(SVCallerTestDataProvider.chrForLongContig1, 20152030, 20154634), TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 3604, contigSequence.length);
        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, SVConstants.CallingStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, SVConstants.CallingStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test
    public void testFilterByNextAlignmentMayBeNovelInsertion() throws Exception {
        AlignmentRegion overlappingRegion1 = new AlignmentRegion("overlap", "22", new SimpleInterval("19", 48699881, 48700035), TextCigarCodec.decode("47S154M"), false, 60, 0, 1, 154);
        AlignmentRegion overlappingRegion2 = new AlignmentRegion("overlap", "22", new SimpleInterval("19", 48700584, 48700669), TextCigarCodec.decode("116H85M"), true, 60, 0, 117, 201);

        Assert.assertTrue(ChimericAlignment.nextAlignmentMayBeNovelInsertion(overlappingRegion1, overlappingRegion2, 50));
    }

    @Test
    public void testBooleanStates_inversion() {

        Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData = SVCallerTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.REVERSE_TO_FORWARD, false, true, true);

        testData = SVCallerTestDataProvider.forSimpleInversionWithHom_leftPlus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.FORWARD_TO_REVERSE, false, true, true);

        testData = SVCallerTestDataProvider.forSimpleInversionWithHom_leftMinus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.FORWARD_TO_REVERSE, true, true, false);

        testData = SVCallerTestDataProvider.forSimpleInversionWithHom_rightPlus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.REVERSE_TO_FORWARD, false, true, true);

        testData = SVCallerTestDataProvider.forSimpleInversionWithHom_rightMinus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.REVERSE_TO_FORWARD, true, true, false);
    }

    @Test
    public void testBooleanStates_simpleInsertionAndDeletion() {

        // simple deletion
        Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData = SVCallerTestDataProvider.forSimpleDeletion_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forSimpleDeletion_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // simple insertion
        testData = SVCallerTestDataProvider.forSimpleInsertion_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forSimpleInsertion_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // long range substitution
        testData = SVCallerTestDataProvider.forLongRangeSubstitution_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forLongRangeSubstitution_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // simple deletion with homology
        testData = SVCallerTestDataProvider.forDeletionWithHomology_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forDeletionWithHomology_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);
    }

    @Test
    public void testBooleanStates_tandemDuplication_simple() {

        // tandem duplication simple contraction
        Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData = SVCallerTestDataProvider.forSimpleTanDupContraction_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forSimpleTanDupContraction_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // tandem duplication simple expansion
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansion_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansion_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // tandem duplication simple expansion with novel insertion
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);
    }

    @Test
    public void testBooleanStates_tandemDuplication_complex() {

        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> testData = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // third test: contraction from 3 units to 2 units without pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        testData = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, false, true, true);
        testData = SVCallerTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;
        testBooleanSeries(testData._1(), testData._2(), ChimericAlignment.StrandSwitch.NO_SWITCH, true, true, false);
    }

    private static void testBooleanSeries(final AlignmentRegion region1, final AlignmentRegion region2,
                                          final ChimericAlignment.StrandSwitch expectedStrandSwitch,
                                          final boolean expectedRefPositionSwitch,
                                          final boolean expectedIsNotSimpleTranslocation,
                                          final boolean expectedIsForwardStrandRepresentation) {

        Assert.assertEquals(ChimericAlignment.determineStrandSwitch(region1, region2), expectedStrandSwitch);
        Assert.assertEquals(ChimericAlignment.involvesRefPositionSwitch(region1, region2), expectedRefPositionSwitch);
        Assert.assertEquals(ChimericAlignment.isForwardStrandRepresentation(region1, region2, expectedStrandSwitch, expectedRefPositionSwitch), expectedIsForwardStrandRepresentation);
        Assert.assertEquals(ChimericAlignment.isNotSimpleTranslocation(region1, region2, expectedStrandSwitch, expectedRefPositionSwitch), expectedIsNotSimpleTranslocation);
    }
}