package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.TestUtilsForAssemblyBasedSVDiscovery;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType.SupportedType.*;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ContigChimericAlignmentIterativeInterpreter.firstAlignmentIsTooShort;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ContigChimericAlignmentIterativeInterpreter.nextAlignmentMayBeInsertion;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public class ContigChimericAlignmentIterativeInterpreterUnitTest extends GATKBaseTest {

    private final AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV assemblyBasedSVDiscoveryTestDataProviderForSimpleSV = new AssemblyBasedSVDiscoveryTestDataProviderForSimpleSV();
    private final AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints assemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints = new AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints();

    // -----------------------------------------------------------------------------------------------
    // Test on step 1: chimeric alignments extraction
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv")
    public void testFilterByRegionTooSmall() {
        final byte[] contigSequence = AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.LONG_CONTIG1.getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval(AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.chrForLongContig1, 20138007, 20142231), 1, contigSequence.length - 1986, TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval(AssemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.chrForLongContig1, 20152030, 20154634), 3604, contigSequence.length, TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertFalse( firstAlignmentIsTooShort(region1, region2, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( firstAlignmentIsTooShort(region2, region1, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( firstAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( firstAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test(groups = "sv")
    public void testFilterByNextAlignmentMayBeInsertion() {
        final AlignmentInterval overlappingRegion1 = new AlignmentInterval(new SimpleInterval("19", 48699881, 48700034), 1, 154, TextCigarCodec.decode("47S154M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval overlappingRegion2 = new AlignmentInterval(new SimpleInterval("19", 48700584, 48700668), 117, 201, TextCigarCodec.decode("116H85M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertTrue(nextAlignmentMayBeInsertion(overlappingRegion1, overlappingRegion2,  CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, 50,true));
    }

    // -----------------------------------------------------------------------------------------------
    // Test on step 3: turn into variant context from novel adjacency (no tests on step 2 because that was tested by NovelAdjacencyAndAltHaplotypeUnitTest)
    // -----------------------------------------------------------------------------------------------
    @DataProvider
    private Object[][] forInferSimpleTypeFromNovelAdjacency() {
        final ReferenceMultiSparkSource referenceMultiSource = TestUtilsForAssemblyBasedSVDiscovery.b37_reference;
        final List<Object[]> data = new ArrayList<>(20);
        // inversion
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHomology_RightBreakpoint_minus.expectedNovelAdjacencyAndAltSeq, INV.name(), ImmutableSet.of(INV33), referenceMultiSource});

        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForInversionBreakpoints.forSimpleInversionWithHom_leftPlus.expectedNovelAdjacencyAndAltSeq, INV.name(), ImmutableSet.of(INV55), referenceMultiSource});

        // simple deletion
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleDeletion_plus.expectedNovelAdjacencyAndAltSeq, DEL.name(), Collections.emptySet(), referenceMultiSource});

        // simple deletion with homology
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forDeletionWithHomology_minus.expectedNovelAdjacencyAndAltSeq, DEL.name(), Collections.emptySet(), referenceMultiSource});

        // simple insertion
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleInsertion_minus.expectedNovelAdjacencyAndAltSeq, INS.name(), Collections.emptySet(), referenceMultiSource});

        // long range substitution
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forLongRangeSubstitution_fudgedDel_plus.expectedNovelAdjacencyAndAltSeq, DEL.name(), Collections.emptySet(), referenceMultiSource});

        // simple tandem dup contraction from 2 units to 1 unit
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupContraction_plus.expectedNovelAdjacencyAndAltSeq, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING), referenceMultiSource});

        // simple tandem dup expansion from 1 unit to 2 units
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupExpansion_ins_minus.expectedNovelAdjacencyAndAltSeq, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING), referenceMultiSource});

        // simple tandem dup expansion from 1 unit to 2 units and novel insertion
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forSimpleTanDupExpansionWithNovelIns_dup_plus.expectedNovelAdjacencyAndAltSeq, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING), referenceMultiSource});

        // tandem dup expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_1to2_pseudoHom_minus.expectedNovelAdjacencyAndAltSeq, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING), referenceMultiSource});

        // tandem dup contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to1_pseudoHom_plus.expectedNovelAdjacencyAndAltSeq, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING), referenceMultiSource});

        // tandem dup contraction from 3 units to 2 units
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_3to2_noPseudoHom_minus.expectedNovelAdjacencyAndAltSeq, DEL.name(), ImmutableSet.of(DUP_TAN_CONTRACTION_STRING), referenceMultiSource});

        // tandem dup expansion from 2 units to 3 units
        data.add(new Object[]{assemblyBasedSVDiscoveryTestDataProviderForSimpleSV.forComplexTanDup_2to3_noPseudoHom_plus.expectedNovelAdjacencyAndAltSeq, DUP.name(), ImmutableSet.of(DUP_TAN_EXPANSION_STRING), referenceMultiSource});

        return data.toArray(new Object[data.size()][]);
    }
    @Test(groups = "sv", dataProvider = "forInferSimpleTypeFromNovelAdjacency")
    public void testInferSimpleTypeFromNovelAdjacency(final NovelAdjacencyAndAltHaplotype biPathBubble,
                                                      final String expectedTypeString,
                                                      final Set<String> expectedAttributeIDs,
                                                      final ReferenceMultiSparkSource reference) {

        final SvType variant = ContigChimericAlignmentIterativeInterpreter.inferSimpleTypeFromNovelAdjacency(biPathBubble, reference);
        Assert.assertEquals(variant.toString(), expectedTypeString);

        final Set<String> attributeIDs = variant.getTypeSpecificAttributes().keySet();
        Assert.assertEquals(attributeIDs, expectedAttributeIDs);
    }

    // -----------------------------------------------------------------------------------------------
    // following might be legacy tests that could be removed but needs time to investigate (Dec.13/2016)
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv", enabled = false)
    public void testGetAssembledBreakpointsFromAlignmentIntervalsWithOverlappingAlignmentInterval() {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 23102785, 23103303), 1, 519, TextCigarCodec.decode("519M1006S"), true, 60, 1, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 23103196, 23103237), 516, 557, TextCigarCodec.decode("515S42M968S"), false, 60, 2, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region3 = new AlignmentInterval(new SimpleInterval("20", 23103633, 23104602), 556, 1525, TextCigarCodec.decode("555S970M"), true, 60, 3, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final AlignedContig alignedContig = new AlignedContig("asm00001:tig0001", contigSequence, Arrays.asList(region1, region2, region3));
        final List<SimpleChimera> assembledBreakpointsFromAlignmentIntervals = ContigChimericAlignmentIterativeInterpreter.parseOneContig(alignedContig, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict, true, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, true);
        Assert.assertEquals(assembledBreakpointsFromAlignmentIntervals.size(), 1);
        final SimpleChimera simpleChimera = assembledBreakpointsFromAlignmentIntervals.get(0);
        Assert.assertEquals(simpleChimera.sourceContigName, "asm00001:tig0001");
        Assert.assertEquals(simpleChimera.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(simpleChimera.regionWithHigherCoordOnContig, region3);
        Assert.assertEquals(simpleChimera.insertionMappings.size(), 1);
        final String expectedInsertionMappingsString = String.join(AlignmentInterval.PACKED_STRING_REP_SEPARATOR, "516", "557", "20:23103196-23103237", "-", "515S42M968S", "60", "2", "100", "O");
        Assert.assertEquals(simpleChimera.insertionMappings.get(0), expectedInsertionMappingsString);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(simpleChimera, contigSequence, TestUtilsForAssemblyBasedSVDiscovery.b37_seqDict);
        Assert.assertTrue(breakpoints.getComplication().getHomologyForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.getComplication().getInsertedSequenceForwardStrandRep().getBytes(), Arrays.copyOfRange(contigSequence, 519, 555));
        Assert.assertEquals(breakpoints.getAltHaplotypeSequence(), Arrays.copyOfRange(contigSequence, 519, 555));
    }
}
