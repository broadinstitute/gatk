package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
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

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.NovelAdjacencyReferenceLocations.EndConnectionType.*;

public class NovelAdjacencyReferenceLocationsUnitTest extends BaseTest{

        static final BreakpointComplications DEFAULT_BREAKPOINT_COMPLICATIONS = new BreakpointComplications();


    /**
     * Hack to force trigger test data generation.
     */
    @BeforeClass
    private void makeSureDataIsAvailable() {
        if(!SVDiscoveryTestDataProvider.testDataInitialized) {
            new SVDiscoveryTestDataProvider();
        }
    }

    private static void seeIfItWorks(final NovelAdjacencyReferenceLocations breakpoints, final NovelAdjacencyReferenceLocations.EndConnectionType expectedEndConnectionType,
                                     final SimpleInterval expectedLeftBreakpoint, final SimpleInterval expectedRightBreakpoint,
                                     final SimpleInterval expectedRepeatUnitRefSpan, final String expectedHomology, final String expectedInsertion,
                                     final int expectedRefDupNum, final int expectedCtgDupNum,
                                     final List<String> expectedTandupCigarStrings) {

        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, expectedLeftBreakpoint);
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, expectedRightBreakpoint);
        Assert.assertEquals(breakpoints.endConnectionType, expectedEndConnectionType);
        Assert.assertEquals(breakpoints.complication.getHomologyForwardStrandRep(), expectedHomology);
        Assert.assertEquals(breakpoints.complication.getInsertedSequenceForwardStrandRep(), expectedInsertion);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), expectedRepeatUnitRefSpan);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), expectedRefDupNum);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), expectedCtgDupNum);
        Assert.assertEquals(breakpoints.complication.getCigarStringsForDupSeqOnCtg(), expectedTandupCigarStrings);

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, breakpoints);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final NovelAdjacencyReferenceLocations roundTrip = (NovelAdjacencyReferenceLocations) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, breakpoints);
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for generic functions on the base class
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv")
    public void testEqualsAndHashCode() throws Exception {

        final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations1 = getBreakpoints("asm00001:tig0001", "foo");

        final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations2 = getBreakpoints("asm00002:tig0002", "bar");

        Assert.assertTrue(novelAdjacencyReferenceLocations1.equals(novelAdjacencyReferenceLocations2));
        Assert.assertEquals(novelAdjacencyReferenceLocations1.hashCode(), novelAdjacencyReferenceLocations2.hashCode());
    }

    @Test(groups = "sv")
    void testKryoSerializer() throws IOException {
        // uses inversion subclass for testing
        final NovelAdjacencyReferenceLocations novelAdjacencyReferenceLocations1 = getBreakpoints("asm00001:tig0001", "foo");
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, novelAdjacencyReferenceLocations1);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final NovelAdjacencyReferenceLocations roundTrip = (NovelAdjacencyReferenceLocations) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, novelAdjacencyReferenceLocations1);
    }

    private static NovelAdjacencyReferenceLocations getBreakpoints(final String contigName, final String insertionMapping) throws IOException{
        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 10000, 10100), 1, 100, TextCigarCodec.decode("100M"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 20100, 20200), 101, 200, TextCigarCodec.decode("100M"), false, 60, 0);
        final ArrayList<String> insertionMappings = new ArrayList<>();
        insertionMappings.add(insertionMapping);
        final ChimericAlignment breakpoint = new ChimericAlignment(region1, region2, insertionMappings, contigName);
        return new NovelAdjacencyReferenceLocations(breakpoint, SVDiscoveryTestDataProvider.makeDummySequence(200, (byte)'A'));
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for inversion
    // -----------------------------------------------------------------------------------------------
    @Test(groups = "sv")
    public void testGetBreakpoints_5to3Inversion_withSimpleInsertion() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleInversionWithNovelInsertion._3();
        seeIfItWorks(breakpoints, FIVE_TO_FIVE,
                new SimpleInterval("21", 108569294, 108569294), new SimpleInterval("21", 108569364, 108569364),
                null, "", "T", 0, 0, null);
    }

    @Test(groups = "sv")
    public void testGetAssembledBreakpointFromAlignmentIntervalsStrangeLeftBreakpoint() throws Exception {

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint._3();
        seeIfItWorks(breakpoints, THREE_TO_THREE,
                new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20138006, 20138006),
                new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20152650, 20152650),
                null, "TGAGGTCAGGAGTTCCTGATCCCATCTTTACTAAAAATACAAAACTTACCCAGGGTGGTTGTGCACACTTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACAAGGGAGGAAATGGTTGCAGTGAGCCATGATCATGCCACTGAACCCCAGCCTGGGCAAGAGAGTGAGACTGTCTCAAAAAAAAAAAAAACTGTTTAATTTTTATGAATGCAGGTTTTCTGCAAACACTACACATAACTATGCTAATTGTTCTGAAGTAATAAATAGAAAGCAAGGCACAACTACAGACTCCACTGTTCAGTTTATGCACTGAACTGTTCTTGCTTTTGCAGTGTAAGTATTTCTGCCTGCAAATACTGGATAATTACCTTGGATCATCAGATTTCTATCAAAGGAATTTAGTATCTTTTAGTCTTTATCATTTTGTATTGCTAAATTTATCTGTGTGTTAAGCTTCTGTGTGCTCTTAAAATGAGGTTTTATCTAAACAAACCTGTGTCTACTTTAAAAGACTAAACATGAAAAAACTAAACTTTTCAGAACCAAAAACAAAGCAATAAATCTGAAGTACTAGATAGTCTGGAGTGAGATTTATTTAGCTTTTTT",
                "", 0, 0, null);
    }

    /**
     *  @see SVDiscoveryTestDataProvider#forSimpleInversionWithHomology(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testSimpleInversionWithHomologyBreakpointsIdentification_allFourRepresentations() throws IOException {

        final byte[] homology = "ACACA".getBytes();

        // left flanking forward strand
        final NovelAdjacencyReferenceLocations breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand = SVDiscoveryTestDataProvider.forSimpleInversionWithHom_leftPlus._3();
        seeIfItWorks(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand, FIVE_TO_FIVE,
                new SimpleInterval("20", 200, 200), new SimpleInterval("20", 605, 605),
                null, new String(homology), "", 0, 0, null);

        // see if reverse strand changes anything
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand, SVDiscoveryTestDataProvider.forSimpleInversionWithHom_leftMinus._3());

        // see if right flanking evidence give the same breakpoint location and homology (up to RC)
        // and see if the two strands give the same result
        final NovelAdjacencyReferenceLocations breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand = SVDiscoveryTestDataProvider.forSimpleInversionWithHom_rightPlus._3();
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.leftJustifiedLeftRefLoc,
                breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.leftJustifiedLeftRefLoc);
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.leftJustifiedRightRefLoc,
                breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.leftJustifiedRightRefLoc);
        Assert.assertEquals(breakpointsIdentifiedFromLeftFlankingEvidenceAndForwardStrand.complication.getHomologyForwardStrandRep(),
                new String(SVDiscoveryTestDataProvider.getReverseComplimentCopy(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.complication.getHomologyForwardStrandRep().getBytes())));
        Assert.assertEquals(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand.endConnectionType, THREE_TO_THREE);
        Assert.assertEquals(breakpointsIdentifiedFromRightFlankingEvidenceAndForwardStrand, SVDiscoveryTestDataProvider.forSimpleInversionWithHom_rightMinus._3());
    }

    @Test(groups = "sv")
    public void testGetAssembledBreakpointsFromAlignmentIntervalsWithOverlappingAlignmentInterval() throws Exception {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 23102817, 23103304), 1, 487, TextCigarCodec.decode("487M1006S"), true, 60, 1);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 23103196, 23103238), 484, 525, TextCigarCodec.decode("483S42M968S"), false, 60, 2);
        final AlignedAssembly.AlignmentInterval region3 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 23103633, 23104603), 524, 1493, TextCigarCodec.decode("523S970M"), true, 60, 3);

        final AlignedContig alignedContig = new AlignedContig("asm00001:tig0001", contigSequence, Arrays.asList(region1, region2, region3));
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentIntervals = ChimericAlignment.fromSplitAlignments(alignedContig, SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH);
        Assert.assertEquals(assembledBreakpointsFromAlignmentIntervals.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentIntervals.get(0);
        Assert.assertEquals(chimericAlignment.sourceContigName, "asm00001:tig0001");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region3);
        Assert.assertEquals(chimericAlignment.insertionMappings.size(), 1);
        final String expectedInsertionMappingsString = String.join(AlignedAssembly.AlignmentInterval.PACKED_STRING_REP_SEPARATOR, "484", "525", "20:23103196-23103238", "-", "483S42M968S", "60", "2");
        Assert.assertEquals(chimericAlignment.insertionMappings.get(0), expectedInsertionMappingsString);//"1_contig-1_484_525_20_23103196_-_483S42M968S_60_2");
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(chimericAlignment, contigSequence);
        Assert.assertTrue(breakpoints.complication.getHomologyForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getInsertedSequenceForwardStrandRep(), "TGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCA");
    }

    // following might be legacy tests that could be removed but needs time to investigate (Dec.13/2016)
    @Test(groups = "sv")
    public void testStrandedness_Inversion() throws IOException {

        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 10000, 10100), 1, 100, TextCigarCodec.decode("100M"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 20100, 20200), 101, 200, TextCigarCodec.decode("100M"), false, 60, 0);
        final ChimericAlignment breakpoint1 = new ChimericAlignment(region1, region2, new ArrayList<>(), "1");

        Assert.assertNotEquals(NovelAdjacencyReferenceLocations.determineEndConnectionType(breakpoint1), NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_THREE);

        final AlignedAssembly.AlignmentInterval region3 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("21", 38343346, 38343483), 1, 137, TextCigarCodec.decode("137M141S"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region4 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 38342908, 38343049), 138, 278, TextCigarCodec.decode("137S141M"), false, 60, 0);
        final ChimericAlignment breakpoint2 = new ChimericAlignment(region3, region4, new ArrayList<>(), "4");

        Assert.assertEquals(NovelAdjacencyReferenceLocations.determineEndConnectionType(breakpoint2), NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_FIVE);

        final AlignedAssembly.AlignmentInterval region5 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("21", 38343346, 38343483), 1, 137, TextCigarCodec.decode("137M141S"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region6 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("21", 38342908, 38343049), 138, 278, TextCigarCodec.decode("137S141M"), false, 60, 0);
        final ChimericAlignment breakpoint3 = new ChimericAlignment(region5, region6, new ArrayList<>(), "3");

        Assert.assertEquals(NovelAdjacencyReferenceLocations.determineEndConnectionType(breakpoint3), NovelAdjacencyReferenceLocations.EndConnectionType.FIVE_TO_FIVE);
    }

    @Test(groups = "sv")
    public void testGetBreakpoints_5to3Inversion_simple() throws IOException {
        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 101, 200), 1, 100, TextCigarCodec.decode("100M100S"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 501, 600), 101, 200, TextCigarCodec.decode("100S100M"), false, 60, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, Collections.emptyList(), "1");
        final Tuple2<SimpleInterval, SimpleInterval> breakpoints = NovelAdjacencyReferenceLocations.leftJustifyBreakpoints(chimericAlignment, DEFAULT_BREAKPOINT_COMPLICATIONS);
        Assert.assertEquals(breakpoints._1(), new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpoints._2(), new SimpleInterval("20", 600, 600));
    }

    @Test(groups = "sv")
    public void testGetBreakpoints_5to3Inversion_withSimpleHomology() throws IOException {
        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 101, 205), 1, 105, TextCigarCodec.decode("105M100S"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 501, 605), 96, 200, TextCigarCodec.decode("105M100S"), false, 60, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, Collections.emptyList(), "1");
        final BreakpointComplications homologyComplications = new BreakpointComplications("ACACA", "",
                DEFAULT_BREAKPOINT_COMPLICATIONS.hasDuplicationAnnotation(), DEFAULT_BREAKPOINT_COMPLICATIONS.getDupSeqRepeatUnitRefSpan(), DEFAULT_BREAKPOINT_COMPLICATIONS.getDupSeqRepeatNumOnRef(), DEFAULT_BREAKPOINT_COMPLICATIONS.getDupSeqRepeatNumOnCtg(), DEFAULT_BREAKPOINT_COMPLICATIONS.getCigarStringsForDupSeqOnCtg(), DEFAULT_BREAKPOINT_COMPLICATIONS.isDupAnnotIsFromOptimization());
        final Tuple2<SimpleInterval, SimpleInterval> breakpoints = NovelAdjacencyReferenceLocations.leftJustifyBreakpoints(chimericAlignment, homologyComplications);
        Assert.assertEquals(breakpoints._1(), new SimpleInterval("20", 200, 200));
        Assert.assertEquals(breakpoints._2(), new SimpleInterval("20", 605, 605));
    }

    // -----------------------------------------------------------------------------------------------
    // Tests for complication resolving and breakpoint justification with the inferred complications for insertion and deletion
    // -----------------------------------------------------------------------------------------------
    @Test(expectedExceptions = GATKException.class)
    public void testGetBreakpoints_ExpectException() throws IOException {
        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("21", 100001, 100100), 1 ,100, TextCigarCodec.decode("100M"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("21", 100101, 100200), 101 ,200, TextCigarCodec.decode("100M"), true, 60, 0);
        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, Collections.emptyList(), "1");
        new BreakpointComplications(chimericAlignment, SVDiscoveryTestDataProvider.makeDummySequence(200, (byte)'A'));
    }

    /**
     * @see SVDiscoveryTestDataProvider#forSimpleDeletion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_simpleDeletion() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleDeletion_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVDiscoveryTestDataProvider.forSimpleDeletion_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100060, 100060),
                null, "", "", 0, 0, null);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVDiscoveryTestDataProvider#forSimpleInsertion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_simpleInsertion() throws IOException {

        byte[] insertedSeq  = SVDiscoveryTestDataProvider.makeDummySequence(50, (byte)'C');

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleInsertion_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVDiscoveryTestDataProvider.forSimpleInsertion_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100100, 100100), new SimpleInterval("21", 100100, 100100),
                null, "", new String(insertedSeq), 0, 0, null);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVDiscoveryTestDataProvider#forLongRangeSubstitution()
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_longRangeSubstitution() throws IOException {

        final byte[] substitution = SVDiscoveryTestDataProvider.makeDummySequence(10, (byte)'C');
        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forLongRangeSubstitution_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVDiscoveryTestDataProvider.forLongRangeSubstitution_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100060, 100060),
                null, "", new String(substitution), 0, 0, null);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVDiscoveryTestDataProvider#forDeletionWithHomology(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_simpleDeletion_withHomology() throws IOException {

        final byte[] homology = "ATCG".getBytes();

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forDeletionWithHomology_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVDiscoveryTestDataProvider.forDeletionWithHomology_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100078, 100078),
                null, new String(homology), "", 0, 0, null);
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVDiscoveryTestDataProvider#forSimpleTandemDuplicationContraction()
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_contraction_simple() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100050, 100050),
                new SimpleInterval("21", 100041, 100050),
                new String(SVDiscoveryTestDataProvider.makeDummySequence(10, (byte)'C')), "",
                2, 1, Collections.emptyList());
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVDiscoveryTestDataProvider#forSimpleTandemDuplicationExpansion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_expansion_simple() throws IOException {

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 100040, 100040), new SimpleInterval("21", 100040, 100040),
                new SimpleInterval("21", 100041, 100050),
                "", "",
                1, 2, Arrays.asList("10M", "10M"));
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVDiscoveryTestDataProvider#forSimpleTandemDuplicationExpansionWithNovelInsertion(ByteArrayOutputStream)
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_expansion_andNovelInsertion() throws IOException {

        final byte[] insertedSeq = "CTCTCTCTCT".getBytes();                                                                           //10
        final byte[] dup = "AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG".getBytes();    //89

        final NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus._3();
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus._3();

        seeIfItWorks(breakpoints, FIVE_TO_THREE, new SimpleInterval("21", 25297163, 25297163), new SimpleInterval("21", 25297163, 25297163),
                new SimpleInterval("21", 25297164,25297252),
                "", new String(insertedSeq),
                1, 2, Arrays.asList("89M", "89M"));
        Assert.assertEquals(breakpointsDetectedFromReverseStrand, breakpoints);
    }

    /**
     * @see SVDiscoveryTestDataProvider#forComplexTandemDuplication()
     */
    @Test(groups = "sv")
    public void testGetBreakpoints_tandemDuplication_complex() throws IOException {

        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";                                                                    // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";                                                            // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA";   // 96
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA";   // 96
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                                                                      // 13


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        NovelAdjacencyReferenceLocations breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus._3();

        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 1);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 1);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 1);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+pseudoHomology)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 1);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        // third test: contraction from 3 units to 2 units without pseudo-homology
        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+secondRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 3);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat+secondRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 3);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 2);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312705, 312705));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 3);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);

        breakpoints = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus._3();
        Assert.assertTrue(StringUtils.getLevenshteinDistance(breakpoints.complication.getHomologyForwardStrandRep(), firstRepeat)<=2);
        Assert.assertTrue(breakpoints.complication.getInsertedSequenceForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatUnitRefSpan(), new SimpleInterval("20", 312610, 312705));
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnRef(), 2);
        Assert.assertEquals(breakpoints.complication.getDupSeqRepeatNumOnCtg(), 3);
        Assert.assertEquals(breakpoints.leftJustifiedLeftRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.leftJustifiedRightRefLoc, new SimpleInterval("20", 312609, 312609));
        Assert.assertEquals(breakpoints.endConnectionType, FIVE_TO_THREE);
    }

    /**
     * todo : see below
     * Tests a simple inter-chromosomal translocation event, which can not be handled yet [4/4/2017];
     *   enable (and definitely expand test cases) when such variants can be handled
     * @throws Exception
     */
    @Test(enabled = false)
    public void testGetAssembledBreakpointsFromAlignmentIntervals() throws Exception {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignedAssembly.AlignmentInterval region1 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("21", 118873207, 118873739), 1, 532, TextCigarCodec.decode("532M87S"), true, 60, 0);
        final AlignedAssembly.AlignmentInterval region2 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 175705642, 175705671), 519, 547, TextCigarCodec.decode("518S29M72S"), false, 3, 0);
        final AlignedAssembly.AlignmentInterval region3 = new AlignedAssembly.AlignmentInterval(new SimpleInterval("20", 118875262, 118875338), 544, 619, TextCigarCodec.decode("543S76M"), false, 60, 0);
        final AlignedContig alignedContig = new AlignedContig("asm00001:tig0001", contigSequence, Arrays.asList(region1, region2, region3));
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentIntervals = ChimericAlignment.fromSplitAlignments(alignedContig, SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH);
        Assert.assertEquals(assembledBreakpointsFromAlignmentIntervals.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentIntervals.get(0);
        Assert.assertEquals(chimericAlignment.sourceContigName, "asm00001:tig0001");
        Assert.assertEquals(chimericAlignment.regionWithLowerCoordOnContig, region1);
        Assert.assertEquals(chimericAlignment.regionWithHigherCoordOnContig, region3);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(chimericAlignment, contigSequence);
        Assert.assertTrue(breakpoints.complication.getHomologyForwardStrandRep().isEmpty());
        Assert.assertEquals(breakpoints.complication.getInsertedSequenceForwardStrandRep(), "GAGATAGAGTC");
    }
}