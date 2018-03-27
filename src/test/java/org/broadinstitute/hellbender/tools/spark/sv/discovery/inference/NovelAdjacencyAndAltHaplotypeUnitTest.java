package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SVTestUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.*;


public class NovelAdjacencyAndAltHaplotypeUnitTest extends GATKBaseTest {


    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!testDataInitialized) {
            new SimpleSVDiscoveryTestDataProvider();
        }
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for generic functions on the base class (uses inversion subclass for testing)
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv")
    public void testEqualsAndHashCode() {

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype1 = getBreakpoints("asm00001:tig0001", "foo");

        final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype2 = getBreakpoints("asm00002:tig0002", "bar");

        Assert.assertEquals(novelAdjacencyAndAltHaplotype1, novelAdjacencyAndAltHaplotype2);
        Assert.assertEquals(novelAdjacencyAndAltHaplotype1.hashCode(), novelAdjacencyAndAltHaplotype2.hashCode());
    }

    @Test(groups = "sv")
    void testKryoSerializer() throws IOException {
        try (final ByteArrayOutputStream bos = new ByteArrayOutputStream()) {
            final Output out = new Output(bos);
            final Kryo kryo = new Kryo();
            final NovelAdjacencyAndAltHaplotype novelAdjacencyAndAltHaplotype =
                    getBreakpoints("asm00001:tig0001", "foo");
            kryo.writeClassAndObject(out, novelAdjacencyAndAltHaplotype);
            out.flush();

            try ( final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray()) ) {
                final Input in = new Input(bis);
                @SuppressWarnings("unchecked")
                final NovelAdjacencyAndAltHaplotype roundTrip = (NovelAdjacencyAndAltHaplotype) kryo.readClassAndObject(in);
                Assert.assertEquals(roundTrip, novelAdjacencyAndAltHaplotype);
            }
        }
    }

    private static NovelAdjacencyAndAltHaplotype getBreakpoints(final String contigName, final String insertionMapping) {
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 10000, 10100), 1, 100, TextCigarCodec.decode("100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 20100, 20200), 101, 200, TextCigarCodec.decode("100M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final ArrayList<String> insertionMappings = new ArrayList<>();
        insertionMappings.add(insertionMapping);
        final ChimericAlignment breakpoint = new ChimericAlignment(region1, region2, insertionMappings, contigName, b37_seqDict);
        return new NovelAdjacencyAndAltHaplotype(breakpoint, SVTestUtils.makeDummySequence(200, (byte)'A'), b37_seqDict);
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for inversion
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv")
    public void testGetBreakpoints_5to3Inversion_withSimpleInsertion() {

        final NovelAdjacencyAndAltHaplotype breakpoints = forSimpleInversionWithNovelInsertion.biPathBubble;
        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.FORWARD_TO_REVERSE,
                new SimpleInterval("21", 69294, 69294), new SimpleInterval("21", 69364, 69364),
                null, "", "T", 0, 0, null,
                new byte[]{'T'});
    }

    @Test(groups = "sv")
    public void testGetAssembledBreakpointFromAlignmentIntervalsStrangeLeftBreakpoint() {

        final NovelAdjacencyAndAltHaplotype breakpoints = forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint.biPathBubble;
        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.REVERSE_TO_FORWARD,
                new SimpleInterval(chrForLongContig1, 20138007, 20138007),
                new SimpleInterval(chrForLongContig1, 20152651, 20152651),
                null, homologyForLongContig1,
                "", 0, 0, null,
                new byte[0]);
    }

    /**
     *  @see SimpleSVDiscoveryTestDataProvider#forSimpleInversionWithHomology(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testSimpleInversionWithHomologyBreakpointsIdentification_allFourRepresentations() {

        // left flanking forward strand representation
        final NovelAdjacencyAndAltHaplotype breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand =
                forSimpleInversionWithHom_leftPlus.biPathBubble;
        seeIfItWorksForNonSimpleTranslocations(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand, StrandSwitch.FORWARD_TO_REVERSE,
                new SimpleInterval("20", 200, 200), new SimpleInterval("20", 605, 605),
                null, "ACACA", "", 0, 0, null,
                new byte[0]);
        // see if reverse strand changes anything
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand,
                forSimpleInversionWithHom_leftMinus.biPathBubble);

        // right flanking forward strand representation, and reverse strand representation
        final NovelAdjacencyAndAltHaplotype breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand =
                forSimpleInversionWithHom_rightPlus.biPathBubble;
        Assert.assertEquals(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.getStrandSwitch(), StrandSwitch.REVERSE_TO_FORWARD);
        Assert.assertEquals(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand,
                forSimpleInversionWithHom_rightMinus.biPathBubble);

        // two novel adjacencies detected by left and right assemblies differ by exactly one base
        Assert.assertEquals(
                shiftRightOneBase(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.getLeftJustifiedLeftRefLoc()),
                breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.getLeftJustifiedLeftRefLoc());
        Assert.assertEquals(
                shiftRightOneBase(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.getLeftJustifiedRightRefLoc()),
                breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.getLeftJustifiedRightRefLoc());
        // and in this case should have the same homology
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.getComplication().getHomologyForwardStrandRep(),
                breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.getComplication().getHomologyForwardStrandRep());
    }

    private static SimpleInterval shiftRightOneBase(final SimpleInterval toBeShifted) {
        return new SimpleInterval(toBeShifted.getContig(), toBeShifted.getStart()+1, toBeShifted.getEnd()+1);
    }

    @Test(groups = "sv")
    public void testGetAssembledBreakpointsFromAlignmentIntervalsWithOverlappingAlignmentInterval() {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 23102785, 23103304), 1, 519, TextCigarCodec.decode("519M1006S"), true, 60, 1, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 23103196, 23103238), 516, 557, TextCigarCodec.decode("515S42M968S"), false, 60, 2, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region3 = new AlignmentInterval(new SimpleInterval("20", 23103633, 23104603), 556, 1525, TextCigarCodec.decode("555S970M"), true, 60, 3, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final AlignedContig alignedContig = new AlignedContig("asm00001:tig0001", contigSequence, Arrays.asList(region1, region2, region3));
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentIntervals = ChimericAlignment.parseOneContig(alignedContig, b37_seqDict, true, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, true);
        Assert.assertEquals(assembledBreakpointsFromAlignmentIntervals.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentIntervals.get(0);
        Assert.assertEquals(chimericAlignment.sourceContigName, "asm00001:tig0001");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region3);
        Assert.assertEquals(chimericAlignment.insertionMappings.size(), 1);
        final String expectedInsertionMappingsString = String.join(AlignmentInterval.PACKED_STRING_REP_SEPARATOR, "516", "557", "20:23103196-23103238", "-", "515S42M968S", "60", "2", "100", "O");
        Assert.assertEquals(chimericAlignment.insertionMappings.get(0), expectedInsertionMappingsString);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(chimericAlignment, contigSequence, b37_seqDict);
        Assert.assertTrue(breakpoints.getComplication().getHomologyForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.getComplication().getInsertedSequenceForwardStrandRep().getBytes(), Arrays.copyOfRange(contigSequence, 519, 555));
        Assert.assertEquals(breakpoints.getAltHaplotypeSequence(), Arrays.copyOfRange(contigSequence, 519, 555));
    }

    // following might be legacy tests that could be removed but needs time to investigate (Dec.13/2016)
    @Test(groups = "sv")
    public void testGetBreakpoints_5to3Inversion_simple() {
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 101, 200), 1, 100, TextCigarCodec.decode("100M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 501, 600), 101, 200, TextCigarCodec.decode("100S100M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, Collections.emptyList(), "1", b37_seqDict);
        final Tuple2<SimpleInterval, SimpleInterval> breakpoints =
                BreakpointsInference.getInferenceClass(chimericAlignment, null, b37_seqDict)
                .getLeftJustifiedBreakpoints();
        Assert.assertEquals(breakpoints._1(), new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpoints._2(), new SimpleInterval("20", 600, 600));
    }

    @Test(groups = "sv")
    public void testGetBreakpoints_5to3Inversion_withSimpleHomology() {
        final byte[] contigSeq = StringUtils.repeat("C", 50).concat("ATATAT").concat(StringUtils.repeat("C", 50)).getBytes();

        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 101, 156), 1, 56, TextCigarCodec.decode("56M50S"), true, 60, 0, 56, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 501, 556), 51, 106, TextCigarCodec.decode("56M50S"), false, 60, 0, 56, ContigAlignmentsModifier.AlnModType.NONE);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, Collections.emptyList(), "1", b37_seqDict);

        final Tuple2<SimpleInterval, SimpleInterval> breakpoints =
                BreakpointsInference.getInferenceClass(chimericAlignment, contigSeq, b37_seqDict)
                        .getLeftJustifiedBreakpoints();
        Assert.assertEquals(breakpoints._1(), new SimpleInterval("20", 150, 150));
        Assert.assertEquals(breakpoints._2(), new SimpleInterval("20", 556, 556));
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for complication resolving and breakpoint justification with the inferred complications for insertion and deletion
    // -----------------------------------------------------------------------------------------------

    /**
     * @see SimpleSVDiscoveryTestDataProvider#forSimpleDeletion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_simpleDeletion() {

        final NovelAdjacencyAndAltHaplotype breakpoints = forSimpleDeletion_plus.biPathBubble;
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = forSimpleDeletion_minus.biPathBubble;

        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.NO_SWITCH, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100060, 100060),
                null, "", "", 0, 0, null,
                new byte[0]);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SimpleSVDiscoveryTestDataProvider#forSimpleInsertion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_simpleInsertion() {

        final byte[] insertedSeq  = SVTestUtils.makeDummySequence(50, (byte)'C');

        final NovelAdjacencyAndAltHaplotype breakpoints = forSimpleInsertion_plus.biPathBubble;
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = forSimpleInsertion_minus.biPathBubble;

        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.NO_SWITCH, new SimpleInterval("21", 100100, 100100), new SimpleInterval("21", 100100, 100100),
                null, "", new String(insertedSeq), 0, 0, null,
                insertedSeq);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SimpleSVDiscoveryTestDataProvider#forLongRangeSubstitution()
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_longRangeSubstitution() {

        final byte[] substitution = SVTestUtils.makeDummySequence(10, (byte)'C');
        final NovelAdjacencyAndAltHaplotype breakpoints = forLongRangeSubstitution_plus.biPathBubble;
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = forLongRangeSubstitution_minus.biPathBubble;

        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.NO_SWITCH, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100060, 100060),
                null, "", new String(substitution), 0, 0, null,
                substitution);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SimpleSVDiscoveryTestDataProvider#forDeletionWithHomology(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_simpleDeletion_withHomology() {

        final byte[] homology = "ATCG".getBytes();

        final NovelAdjacencyAndAltHaplotype breakpoints = forDeletionWithHomology_plus.biPathBubble;
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = forDeletionWithHomology_minus.biPathBubble;

        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.NO_SWITCH, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100078, 100078),
                null, new String(homology), "", 0, 0, null,
                new byte[0]);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SimpleSVDiscoveryTestDataProvider#forSimpleTandemDuplicationContraction()
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_contraction_simple() {

        final NovelAdjacencyAndAltHaplotype breakpoints = forSimpleTanDupContraction_plus.biPathBubble;
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = forSimpleTanDupContraction_minus.biPathBubble;

        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.NO_SWITCH, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100050, 100050),
                new SimpleInterval("21", 100041, 100050),
                new String(SVTestUtils.makeDummySequence(10, (byte)'C')), "",
                2, 1, Collections.emptyList(),
                SVTestUtils.makeDummySequence(10, (byte)'C'));
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SimpleSVDiscoveryTestDataProvider#forSimpleTandemDuplicationExpansion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_expansion_simple() {

        final NovelAdjacencyAndAltHaplotype breakpoints = forSimpleTanDupExpansion_plus.biPathBubble;
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = forSimpleTanDupExpansion_minus.biPathBubble;

        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.NO_SWITCH, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100040, 100040),
                new SimpleInterval("21", 100041, 100050),
                "", "",
                1, 2, Arrays.asList("10M", "10M"),
                SVTestUtils.makeDummySequence(20, (byte)'C'));
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SimpleSVDiscoveryTestDataProvider#forSimpleTandemDuplicationExpansionWithNovelInsertion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_expansion_andNovelInsertion() {

        final String insertedSeq = "CTCTCTCTCT";                                                                           //10
        final String dup = "AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG";    //89
        final String alt = dup + insertedSeq + dup;

        final NovelAdjacencyAndAltHaplotype breakpoints = forSimpleTanDupExpansionWithNovelIns_plus.biPathBubble;
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = forSimpleTanDupExpansionWithNovelIns_minus.biPathBubble;

        seeIfItWorksForNonSimpleTranslocations(breakpoints, StrandSwitch.NO_SWITCH, new SimpleInterval("21", 25297163, 25297163), new SimpleInterval("21", 25297163, 25297163),
                new SimpleInterval("21", 25297164,25297252),
                "", insertedSeq,
                1, 2, Arrays.asList("89M", "89M"),
                alt.getBytes());
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    private static void seeIfItWorksForNonSimpleTranslocations(final NovelAdjacencyAndAltHaplotype breakpoints, final StrandSwitch expectedStrandSwitch,
                                                               final SimpleInterval expectedLeftBreakpoint, final SimpleInterval expectedRightBreakpoint,
                                                               final SimpleInterval expectedRepeatUnitRefSpan, final String expectedHomology, final String expectedInsertion,
                                                               final int expectedRefDupNum, final int expectedCtgDupNum,
                                                               final List<String> expectedTandupCigarStrings,
                                                               final byte[] expectedAltHaplotypeSequence) {

        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), expectedLeftBreakpoint);
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), expectedRightBreakpoint);
        Assert.assertEquals(breakpoints.getStrandSwitch(), expectedStrandSwitch);
        Assert.assertEquals(breakpoints.getAltHaplotypeSequence(), expectedAltHaplotypeSequence);
        if (expectedStrandSwitch.equals(StrandSwitch.NO_SWITCH)) {
            if (expectedRepeatUnitRefSpan == null) {
                final BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications complication =
                        (BreakpointComplications.SimpleInsDelOrReplacementBreakpointComplications) breakpoints.getComplication();
                Assert.assertEquals(complication.getHomologyForwardStrandRep(), expectedHomology);
                Assert.assertEquals(complication.getInsertedSequenceForwardStrandRep(), expectedInsertion);
            } else {
                final BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications complication =
                        (BreakpointComplications.SmallDuplicationWithPreciseDupRangeBreakpointComplications) breakpoints.getComplication();
                Assert.assertEquals(complication.getHomologyForwardStrandRep(), expectedHomology);
                Assert.assertEquals(complication.getInsertedSequenceForwardStrandRep(), expectedInsertion);
                Assert.assertEquals(complication.getDupSeqRepeatUnitRefSpan(), expectedRepeatUnitRefSpan);
                Assert.assertEquals(complication.getDupSeqRepeatNumOnRef(), expectedRefDupNum);
                Assert.assertEquals(complication.getDupSeqRepeatNumOnCtg(), expectedCtgDupNum);
                Assert.assertEquals(complication.getCigarStringsForDupSeqOnCtgForwardStrandRep(), expectedTandupCigarStrings);
            }
        } else {
            final BreakpointComplications.IntraChrStrandSwitchBreakpointComplications complication =
                    (BreakpointComplications.IntraChrStrandSwitchBreakpointComplications) breakpoints.getComplication();
            Assert.assertEquals(complication.getHomologyForwardStrandRep(), expectedHomology);
            Assert.assertEquals(complication.getInsertedSequenceForwardStrandRep(), expectedInsertion);
            Assert.assertEquals(complication.getDupSeqRepeatUnitRefSpan(), expectedRepeatUnitRefSpan);
            Assert.assertEquals(complication.getDupSeqRepeatNumOnRef(), expectedRefDupNum);
            Assert.assertEquals(complication.getDupSeqRepeatNumOnCtg(), expectedCtgDupNum);
            Assert.assertEquals(complication.getCigarStringsForDupSeqOnCtg(), expectedTandupCigarStrings);
        }

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, breakpoints);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final NovelAdjacencyAndAltHaplotype roundTrip = (NovelAdjacencyAndAltHaplotype) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, breakpoints);
    }


    /**
     * @see SimpleSVDiscoveryTestDataProvider#forComplexTandemDuplication()
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_complex() {

        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";                                                                    // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";                                                            // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA";   // 96
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA";   // 96
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                                                                      // 13


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        NovelAdjacencyAndAltHaplotype breakpoints = forComplexTanDup_1to2_pseudoHom_plus.biPathBubble;
        BreakpointComplications.SmallDuplicationBreakpointComplications complication =
                (BreakpointComplications.SmallDuplicationBreakpointComplications) breakpoints.getComplication();

        Assert.assertTrue(StringUtils.getLevenshteinDistance(complication.getHomologyForwardStrandRep(), pseudoHomology)<=2);
        Assert.assertTrue(complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(complication.getDupSeqRepeatNumOnRef(), 1);
        Assert.assertEquals(complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.getStrandSwitch(), StrandSwitch.NO_SWITCH);
        Assert.assertEquals(breakpoints.getAltHaplotypeSequence(), (firstRepeat+secondRepeat+pseudoHomology).getBytes());
        NovelAdjacencyAndAltHaplotype breakpointsRev = forComplexTanDup_1to2_pseudoHom_minus.biPathBubble;
        Assert.assertEquals(breakpointsRev, breakpoints); // different representation, should lead to same result

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = forComplexTanDup_2to1_pseudoHom_plus.biPathBubble;
        complication = (BreakpointComplications.SmallDuplicationBreakpointComplications) breakpoints.getComplication();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(complication.getHomologyForwardStrandRep(), firstRepeat+pseudoHomology)<=2);
        Assert.assertTrue(complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(complication.getDupSeqRepeatNumOnCtg(), 1);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.getStrandSwitch(), StrandSwitch.NO_SWITCH);
        Assert.assertEquals(breakpoints.getAltHaplotypeSequence(), (firstRepeat+pseudoHomology).getBytes());
        breakpointsRev = forComplexTanDup_2to1_pseudoHom_minus.biPathBubble;
        Assert.assertEquals(breakpointsRev, breakpoints); // different representation, should lead to same result

        // third test: contraction from 3 units to 2 units without pseudo-homology
        breakpoints = forComplexTanDup_3to2_noPseudoHom_plus.biPathBubble;
        complication = (BreakpointComplications.SmallDuplicationBreakpointComplications) breakpoints.getComplication();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(complication.getHomologyForwardStrandRep(), firstRepeat+secondRepeat)<=2);
        Assert.assertTrue(complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(complication.getDupSeqRepeatNumOnRef(), 3);
        Assert.assertEquals(complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.getStrandSwitch(), StrandSwitch.NO_SWITCH);
        Assert.assertEquals(breakpoints.getAltHaplotypeSequence(), (firstRepeat+secondRepeat).getBytes());
        breakpointsRev = forComplexTanDup_3to2_noPseudoHom_minus.biPathBubble;
        Assert.assertEquals(breakpointsRev, breakpoints); // different representation, should lead to same result

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        breakpoints = forComplexTanDup_2to3_noPseudoHom_plus.biPathBubble;
        complication = (BreakpointComplications.SmallDuplicationBreakpointComplications) breakpoints.getComplication();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(complication.getHomologyForwardStrandRep(), firstRepeat)<=2);
        Assert.assertTrue(complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(complication.getDupSeqRepeatNumOnCtg(), 3);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.getStrandSwitch(), StrandSwitch.NO_SWITCH);
        Assert.assertEquals(breakpoints.getAltHaplotypeSequence(), (firstRepeat+secondRepeat+firstRepeat).getBytes());
        breakpointsRev = forComplexTanDup_2to3_noPseudoHom_minus.biPathBubble;
        Assert.assertEquals(breakpointsRev, breakpoints); // different representation, should lead to same result

    }

    @Test(groups = "sv")
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
        ChimericAlignment simpleChimera = new ChimericAlignment(region1, region2, Collections.emptyList(), "testContig", b38_seqDict);
        NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(simpleChimera,
                "TTCCTTAAAATGCAGGTGAATACAAGAATTAGGTTTCAGGTTTTATATATATATTCTGATATATATATATAATATAACCTGAGATATATATATAAATATATATATTAATATATATTAATATATATAAATATATATATATTAATATATATTTATATATAAATATATATATATTAATATATATAAATATATATAAATATATATATATTAATATATATTAATATATAAATATATATATATTAATATATATTAATATATATAAATATATATATTAATATATATAAATATATATATAAATATATATAAATATATAAATATATATATAAATATATATAAATATATATAAATATATATACACACATACATACACATATACATT".getBytes(),
                b38_seqDict);
        Assert.assertEquals(breakpoints.getLeftJustifiedLeftRefLoc(), new SimpleInterval("chr21", 39192594, 39192594));
        Assert.assertEquals(breakpoints.getLeftJustifiedRightRefLoc(), new SimpleInterval("chr21", 39477346, 39477346));
        Assert.assertEquals(breakpoints.getComplication().getHomologyForwardStrandRep(), "ATATATAAATATATATA");
        Assert.assertTrue(breakpoints.getComplication().getInsertedSequenceForwardStrandRep().isEmpty());
    }
}