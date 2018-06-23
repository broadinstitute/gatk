package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME;

public class BreakpointComplicationsAndLocationAndAltSeqInferenceUnitTest extends AssemblyBasedSVDiscoveryBaseTest {

    // -----------------------------------------------------------------------------------------------
    // Tests for utility functions
    // -----------------------------------------------------------------------------------------------
    @DataProvider
    private Object[][] forValidateInferredLocations() {
        final List<Object[]> data = new ArrayList<>(20);
        final SAMSequenceDictionary bareBoneHg38SAMSeqDict = TestUtilsForAssemblyBasedSVDiscovery.bareBoneHg38SAMSeqDict;

        SimpleInterval inferredLeftBreakpoint = new SimpleInterval("chr20", 10000, 10000);
        SimpleInterval inferredRightBreakpoint = new SimpleInterval("chr20", 10000-1, 10000-1);
        data.add(new Object[]{inferredLeftBreakpoint, inferredRightBreakpoint, bareBoneHg38SAMSeqDict});

        inferredLeftBreakpoint = new SimpleInterval("chr20", 64444167+1, 64444167+1);
        inferredRightBreakpoint = new SimpleInterval("chr20", 64444167+1, 64444167+1);
        data.add(new Object[]{inferredLeftBreakpoint, inferredRightBreakpoint, bareBoneHg38SAMSeqDict});

        inferredLeftBreakpoint = new SimpleInterval("chr20", 10000, 10000);
        inferredRightBreakpoint = new SimpleInterval("chr21", 46709983+1, 46709983+1);

        data.add(new Object[]{inferredLeftBreakpoint, inferredRightBreakpoint, bareBoneHg38SAMSeqDict});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(expectedExceptions = GATKException.ShouldNeverReachHereException.class, dataProvider = "forValidateInferredLocations")
    public void testValidateInferredLocations(final SimpleInterval leftBreakpoint,
                                              final SimpleInterval rightBreakpoint,
                                              final SAMSequenceDictionary referenceSequenceDictionary) {
        BreakpointsInference.validateInferredLocations(leftBreakpoint, rightBreakpoint, referenceSequenceDictionary, "");
    }

    @DataProvider
    private Object[][] forGetInferenceClass() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera assemblyBasedSVDiscoveryTestDataForSimpleChimera : getAllTestData()) {
            final SAMSequenceDictionary appropriateDictionary = assemblyBasedSVDiscoveryTestDataForSimpleChimera.getAppropriateDictionary();
            SimpleChimera simpleChimera = new SimpleChimera(assemblyBasedSVDiscoveryTestDataForSimpleChimera.firstAlignment, assemblyBasedSVDiscoveryTestDataForSimpleChimera.secondAlignment,
                    Collections.emptyList(), assemblyBasedSVDiscoveryTestDataForSimpleChimera.evidenceAssemblyContigName, NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME,
                    appropriateDictionary);
            data.add(new Object[]{simpleChimera, assemblyBasedSVDiscoveryTestDataForSimpleChimera.evidenceContigSeq,
                    appropriateDictionary, assemblyBasedSVDiscoveryTestDataForSimpleChimera.getAppropriateBreakpointInferencer()
            });
        }
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forGetInferenceClass")
    public void testGetInferenceClass(final SimpleChimera simpleChimera, final byte[] contigSeq,
                                      final SAMSequenceDictionary samSequenceDictionary,
                                      final Class<? extends BreakpointsInference> expectedInferencerClass) {
        Assert.assertEquals(BreakpointsInference.getInferenceClass(simpleChimera, contigSeq, samSequenceDictionary).getClass(),
                expectedInferencerClass);
    }

    /**
     * To add complexity in test data.
     */
    @Test(groups = "sv")
    public void testExtractCigarAndAltSeqForSimpleTandup_complicated() {

        final int contigTotalLength = 355;

        // forward strand
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("1", 1000001, 1000125), 16, 90,
                TextCigarCodec.decode("5H10S15M20D25M30D35M260S5H"),
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("1", 1000041, 1000145), 191, 345,
                TextCigarCodec.decode("5H185S45M30I55M20I5M10S5H"),
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);


        final Cigar cigar1 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region1, 1000125, 1000041);
        Assert.assertEquals(cigar1, TextCigarCodec.decode("20M30D35M"));
        final Cigar cigar2 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region2, 1000125, 1000041);
        Assert.assertEquals(cigar2, TextCigarCodec.decode("45M30I40M"));

        // reverse strand
        final AlignmentInterval region3 = new AlignmentInterval(region2.referenceSpan, contigTotalLength-region2.endInAssembledContig+1, contigTotalLength-region2.startInAssembledContig+1,
                CigarUtils.invertCigar(region2.cigarAlong5to3DirectionOfContig),
                false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region4 = new AlignmentInterval(region1.referenceSpan, contigTotalLength-region1.endInAssembledContig+1, contigTotalLength-region1.startInAssembledContig+1,
                CigarUtils.invertCigar(region1.cigarAlong5to3DirectionOfContig),
                false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final Cigar cigar3 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region3, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar3), cigar2);
        final Cigar cigar4 = BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications.extractCigarForTandupExpansion(region4, 1000125, 1000041);
        Assert.assertEquals(CigarUtils.invertCigar(cigar4), cigar1);

        // modifying the real event below
        //  "asm030282:tig00005     1_342_chrX:1294785-1295169_-_168M43D174M1321H_60_65_173_O       463_1663_chrX:1293941-1295139_-_462S1098M2I101M_60_17_1106_O";
        final AlignmentInterval region5 = new AlignmentInterval(new SimpleInterval("chr20", 1294785, 1295169),
                1, 342, TextCigarCodec.decode("168M43D174M1321H"),
                false, 60, 65, 173, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region6 = new AlignmentInterval(new SimpleInterval("chr20", 1293941, 1295139),
                463, 1663, TextCigarCodec.decode("462S1098M2I101M"),
                false, 60, 17, 1106, ContigAlignmentsModifier.AlnModType.NONE);
        final SimpleChimera simpleChimera = new SimpleChimera(region5, region6, Collections.emptyList(), "asm030282:tig00005",
                AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME,
                TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);
        final BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications smallDuplicationWithPreciseDupRangeBreakpointComplications =
                new BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications(simpleChimera, "CTCCTGTCGTCCAGGTAGACATGGAGCAGGCTCTCCTTCTTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGATCTCCTCTTGTCTACACTGGATCGAGGTGGATGTGAGGCAGTCTCTCTTCCTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCCAGGTAGCTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCTCCATGGAGGTGGGGCAGGGTCTCCTTCTGTCTACACTGCGTGTAGTTGGAGGTGGGGCAGGGTCTCCTCCTGTCTACACTGGGTCCAGGTAGACATGGGGCAGTCTCTCCTTCTTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCCAGGTAGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCTCCATGGAGGTGGGGCAGTCTCTCCTTCTTGTCTACACTGGGTCCAGGTGGTTGTGGGACAAGATCTCCTCTTGTCTACACTGGCTCGAGGTGGACATGGGGCAGGGTCTCTTCTTGTCTACACTGGGTCCAGGAGGTTGTGAGACAAGATCTCCTCTTGTCTACACTGGATCGAGGTGGACGTGAGGCAGTCTCTCTCCCTGTCTACACTGGGTCCAGGTAGTTGTGGGACAAGAGCTCCTCCTGTCTACACTGGGTCTCCATGGAGGTGGAGCAGGGTCTCCTCCTGTCTACACTGGGTGTAGGTGGAGGTGGGGTCGGGTGTCCTCCTATCTACACTGGGTCCAGGTAGACATGGGGCAGGGTCTCCTTCTCTCTACACTGCGTCCAGCTGGAGGTGGAGCAGAGGCTCTCCTTGCTTGTGGCATCGTCCCCCACACCTCCCGGTCCACTTCCTGGTTCCATGGTTGCAGGATCATCCTTGTCCACCCTCCCTGCAACCTCTTTCAAGGTGGCTCCACAGGCCACAGACCCTTCACCTCTTCCTCCGCTACCCGAAGTGTGTTCACCCCAGAGTCACCGCTCACCACCCACCCATCCTTCCCCCAGGCCACTTCCCCGGGATTCCCAGGCTCCTGTGCGGGCGTGTCCCGTACGCCTCCCTCCTGGTGCCCAGCCCCGGGGAGCTCTCACCGACCTTTCTGTGGACGTCCAGCTGGTACTGAAAGTCCAGGTACGACAGCTTCTGATAGGTCCTGGGCTGTTTCCACCGTACGAGGCAGTGCGTCGTGTTGCAACGTACGGTGACATTGCTGGGAGGGTTGAATCGTTCTGTAACGAGGGCGCAGGACACACCCCTGAACCCGAGAGGTCCTGTCTACACTGGGTCCAGGTGGAGGTGGTGCAGAGTCTCCTCCTGTCTACACTGGGTCCAGGTGGAGGTGGAGTAGGTCCTGTCTACACTGGGTCCAGGTGGAGGTGGAGTAGGGACACCTCTTTGACTACACTGGGTCCAGGTGGAGATGGGGCAGGGTCTCCTCCTGTCTACACTGGGTCCAGGTGGAGGTGGGGCAGGGTCTCCTCCTGTCTACACTGGGTCTAGGTGGAGGTGGTGCAGAGTCTTCTCCTGTCTACACTGGGTCCAGGTGGAGGTGGGGCAGGGTCTCCTCCTGTCTACACTCGGTCCAGGTGGATGTGGACTAGGGACACCTCTTTGTCTA".getBytes(), true);
        Assert.assertEquals(smallDuplicationWithPreciseDupRangeBreakpointComplications.getCigarStringsForDupSeqOnCtgForwardStrandRep(),
                Arrays.asList("355M", "174M43D138M"));
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for breakpoint location inference and complication resolving
    // -----------------------------------------------------------------------------------------------
    @DataProvider
    private Object[][] forBreakpointComplicationsAndLocationAndAltSeqInference() {
        final List<Object[]> data = new ArrayList<>(20);
        for (final AssemblyBasedSVDiscoveryTestDataProvider.AssemblyBasedSVDiscoveryTestDataForSimpleChimera assemblyBasedSVDiscoveryTestDataForSimpleChimera : getAllTestData()) {
            final SAMSequenceDictionary appropriateDictionary = assemblyBasedSVDiscoveryTestDataForSimpleChimera.getAppropriateDictionary();
            SimpleChimera simpleChimera = new SimpleChimera(assemblyBasedSVDiscoveryTestDataForSimpleChimera.firstAlignment, assemblyBasedSVDiscoveryTestDataForSimpleChimera.secondAlignment,
                    Collections.emptyList(), assemblyBasedSVDiscoveryTestDataForSimpleChimera.evidenceAssemblyContigName, NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME,
                    appropriateDictionary);
            data.add(new Object[]{simpleChimera, assemblyBasedSVDiscoveryTestDataForSimpleChimera.evidenceContigSeq, appropriateDictionary,
                                  assemblyBasedSVDiscoveryTestDataForSimpleChimera.expectedNovelAdjacencyAndAltSeq
            });
        }
        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forBreakpointComplicationsAndLocationAndAltSeqInference")
    public void testBreakpointComplicationsAndLocationAndAltSeqInference(final SimpleChimera simpleChimera, final byte[] contigSeq,
                                                                         final SAMSequenceDictionary samSequenceDictionary,
                                                                         final NovelAdjacencyAndAltHaplotype expectedNAAAH) {
        final NovelAdjacencyAndAltHaplotype actual = new NovelAdjacencyAndAltHaplotype(simpleChimera, contigSeq, samSequenceDictionary);
        try { // if equality test failed, it is because of slight difference, usually due to slight difference in contig sequence used in opposite strands
            Assert.assertEquals(actual, expectedNAAAH);
        } catch (final AssertionError ex) {
            testEqualLocations(actual, expectedNAAAH);
            testEqualComplicationsAllowingSlightDiff(actual, expectedNAAAH);
        }
    }

    private static void testEqualLocations(final NovelAdjacencyAndAltHaplotype actual, final NovelAdjacencyAndAltHaplotype expected) {
        Assert.assertEquals(actual.getLeftJustifiedLeftRefLoc(), expected.getLeftJustifiedLeftRefLoc());
        Assert.assertEquals(actual.getLeftJustifiedRightRefLoc(), expected.getLeftJustifiedRightRefLoc());
        Assert.assertEquals(actual.getStrandSwitch(), expected.getStrandSwitch());
        Assert.assertEquals(actual.getTypeInferredFromSimpleChimera(), expected.getTypeInferredFromSimpleChimera());
    }
    // doing it this way because some complications will be different slightly but doesn't matter
    private static void testEqualComplicationsAllowingSlightDiff(final NovelAdjacencyAndAltHaplotype actual, final NovelAdjacencyAndAltHaplotype expected) {


        final BreakpointComplications actualComplication = actual.getComplication();
        final BreakpointComplications expectedComplication = expected.getComplication();
        Assert.assertEquals(actualComplication.getClass(), expectedComplication.getClass());
        int levenshteinDistance = StringUtils.getLevenshteinDistance(actualComplication.getHomologyForwardStrandRep(),
                                                                     expectedComplication.getHomologyForwardStrandRep());
        Assert.assertTrue(levenshteinDistance <= 2);
        levenshteinDistance = StringUtils.getLevenshteinDistance(actualComplication.getInsertedSequenceForwardStrandRep(),
                                                                 expectedComplication.getInsertedSequenceForwardStrandRep());
        Assert.assertTrue(levenshteinDistance <= 2);
        levenshteinDistance = StringUtils.getLevenshteinDistance(new String(actual.getAltHaplotypeSequence()),
                                                                 new String(expected.getAltHaplotypeSequence()));
        Assert.assertTrue(levenshteinDistance <= 2);

        // more tests if more complications available
        if (actualComplication instanceof BreakpointComplications.SmallDuplicationBreakpointComplications) {
            BreakpointComplications.SmallDuplicationBreakpointComplications actualSmallDupComplication = (BreakpointComplications.SmallDuplicationBreakpointComplications) actualComplication;
            BreakpointComplications.SmallDuplicationBreakpointComplications expectedSmallDupComplication = (BreakpointComplications.SmallDuplicationBreakpointComplications) expectedComplication;

            Assert.assertEquals(actualSmallDupComplication.getDupSeqOrientationsOnCtg(), expectedSmallDupComplication.getDupSeqOrientationsOnCtg());
            Assert.assertEquals(actualSmallDupComplication.getDupSeqRepeatNumOnCtg(), expectedSmallDupComplication.getDupSeqRepeatNumOnCtg());
            Assert.assertEquals(actualSmallDupComplication.getDupSeqRepeatNumOnRef(), expectedSmallDupComplication.getDupSeqRepeatNumOnRef());
            Assert.assertEquals(actualSmallDupComplication.getDupSeqRepeatUnitRefSpan(), expectedSmallDupComplication.getDupSeqRepeatUnitRefSpan());
            if (actualComplication instanceof BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications) {
                Assert.assertEquals(((BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications) actualComplication).getCigarStringsForDupSeqOnCtgForwardStrandRep(),
                                    ((BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications) expectedComplication).getCigarStringsForDupSeqOnCtgForwardStrandRep());
            } else {
                Assert.assertEquals(((BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications) actualComplication).getImpreciseDupAffectedRefRange(),
                                    ((BreakpointComplications.SmallDuplicationWithImpreciseDupRangeBreakpointComplications) expectedComplication).getImpreciseDupAffectedRefRange());
            }
        } else if (actualComplication instanceof BreakpointComplications.InvertedDuplicationBreakpointComplications) {
            BreakpointComplications.InvertedDuplicationBreakpointComplications actualIntraChrSSComplication = (BreakpointComplications.InvertedDuplicationBreakpointComplications) actualComplication;
            BreakpointComplications.InvertedDuplicationBreakpointComplications expectedIntraChrSSComplication = (BreakpointComplications.InvertedDuplicationBreakpointComplications) expectedComplication;

            Assert.assertEquals(actualIntraChrSSComplication.getDupSeqRepeatUnitRefSpan(), expectedIntraChrSSComplication.getDupSeqRepeatUnitRefSpan());
            Assert.assertEquals(actualIntraChrSSComplication.getDupSeqOrientationsOnCtg(), expectedIntraChrSSComplication.getDupSeqOrientationsOnCtg());
            Assert.assertEquals(actualIntraChrSSComplication.getDupSeqRepeatNumOnRef(), expectedIntraChrSSComplication.getDupSeqRepeatNumOnRef());
            Assert.assertEquals(actualIntraChrSSComplication.getDupSeqRepeatNumOnCtg(), expectedIntraChrSSComplication.getDupSeqRepeatNumOnCtg());
            Assert.assertEquals(actualIntraChrSSComplication.getCigarStringsForDupSeqOnCtg(), expectedIntraChrSSComplication.getCigarStringsForDupSeqOnCtg());
            Assert.assertEquals(actualIntraChrSSComplication.getInvertedTransInsertionRefSpan(), expectedIntraChrSSComplication.getInvertedTransInsertionRefSpan());
            Assert.assertEquals(actualIntraChrSSComplication.isDupAnnotIsFromOptimization(), expectedIntraChrSSComplication.isDupAnnotIsFromOptimization());
        }
    }

    // -----------------------------------------------------------------------------------------------
    // legacy tests (case covered by other more comprehensive tests, but kept in case others need test data)
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv", enabled = false)
    public void testRefOrderSwitch() {
        AlignmentInterval region1 = new AlignmentInterval(
                // assigned from chr18 to chr21 to use the dict
                new SimpleInterval("chr21", 39477098, 39477363),
                1 ,268,
                TextCigarCodec.decode("236M2I30M108S"), true, 32, 25, 133, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(
                new SimpleInterval("chr21", 39192594, 39192692),
                252 ,350,
                TextCigarCodec.decode("251S99M26S"), true, 32, 1, 94, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleChimera simpleChimera = new SimpleChimera(region1, region2, Collections.emptyList(), "testContig", NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME, TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);
        NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(simpleChimera,
                "TTCCTTAAAATGCAGGTGAATACAAGAATTAGGTTTCAGGTTTTATATATATATTCTGATATATATATATAATATAACCTGAGATATATATATAAATATATATATTAATATATATTAATATATATAAATATATATATATTAATATATATTTATATATAAATATATATATATTAATATATATAAATATATATAAATATATATATATTAATATATATTAATATATAAATATATATATATTAATATATATTAATATATATAAATATATATATTAATATATATAAATATATATATAAATATATATAAATATATAAATATATATATAAATATATATAAATATATATAAATATATATACACACATACATACACATATACATT".getBytes(),
                TestUtilsForAssemblyBasedSVDiscovery.b38_seqDict_chr20_chr21);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("chr21", 39192594, 39192594));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("chr21", 39477346, 39477346));
        Assert.assertEquals(breakpoints.getComplication().getHomologyForwardStrandRep(), "ATATATAAATATATATA");
        Assert.assertTrue(breakpoints.getComplication().getInsertedSequenceForwardStrandRep().isEmpty());
    }

}
