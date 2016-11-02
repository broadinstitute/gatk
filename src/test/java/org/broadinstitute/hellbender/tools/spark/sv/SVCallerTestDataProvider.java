package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import scala.Tuple2;
import scala.Tuple3;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Provides test data for testing several methods involved in the SV variant caller.
 * NO TESTS ARE RUN IN THIS PARTICULAR CLASS.
 */
final class SVCallerTestDataProvider {

    static final ReferenceMultiSource reference = new ReferenceMultiSource((PipelineOptions)null, BaseTest.b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    static final SAMSequenceDictionary seqDict = reference.getReferenceSequenceDictionary(null);

    static byte[] getReverseComplimentCopy(final byte[] sequence) {
        final byte[] sequenceCopy = Arrays.copyOf(sequence, sequence.length);
        SequenceUtil.reverseComplement(sequenceCopy);
        return sequenceCopy;
    }

    // the chromosome that the long contig1 is supposed to be mapped to is actually chr19, but to make tests runnable, we could only use "20" or "21"
    // todo: this should be fixed, but since the exact mapped to chromosome is not important now, we push it to later
    static final String chrForLongContig1 = "20";
    static final String LONG_CONTIG1 =
            "TTTTTTTTTTTTTTTCTGAGACCGAATCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACGATCTTGGCTTACTGCAAGCTCTGCCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCCCACCCCCCCACCCCCCCAGGTAGCTG" +
            "GGACTACAGGTGTCTGCCACCACACCTGGCTAAATTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAAGATGGTTTCGCTCTCCTGACCTCGCGATCCGCCCACCTCGGCCTCTCAAAGTGCTGGGATTACAGGCCTGAGCCACTGCGCCC" +
            "TGCCTGACCCCTTCTTTTAAAACAAATCTTTTGTCTTTGTCTTCATTTCTGCATTCGTCCCCTTCGTTCAATCCTGTAGGAATTGACAGTGATATTGGGGAGTTCCCATCCCTGGATTTGGGATTTCCTCGAGTTTCCAGCCCTGTCCTTGTGGCCAAAAAGT" +
            "GGCACAGTTTGGTTTTTGAATCCTGTTCGTGACGCCAAAAATATTCTGTGTGGGAAACATTCAAGGGGAGAAGAAAAGACACACACACAATACCTTTAAGGGTAAATAAGCTTTGTGCCACGTAAATGGAAATGCAGATATAATTAGCAAATTATATAATAAG" +
            "CAAATCAATGTAATAACCAGATTGATATAATAAGCATGTTGATATAATAAGCAAATTGCAATGGGAAGAGGAGAAAGGAAAAGAGATATATATATTTACACTCACCAGACTATGGAGGATTCACCACCAGACTGGGAAGCAACAGCTTGGGCTCCAGAGTCAG" +
            "CTACTCATTCCTGCACAGATGAGGAGGGTCTAATGAAGCTTCAGCACAATCTGATACCCTAGCTCTTTTGTAACGAGTTGTTTGGCATAAGGCCCAGTCATGAGGGCCATTTGCAACTGGGGTCAAGGAACACAAAACTGTCAACTTGTTTTTGCGATTGTCT" +
            "ATTGTTTTTCAACAACTAATATATAGAAATAGATTGAAATAGAGATTTCTCTGAAACAGCGCTGAATGAATGCCTCAAGGGCCTCACACAACCTGTTTGAGAACTTGGTGACTACTGAGTTTGTCCACGTTCAATTAAGTTCAAATTTAGTATTTAACTTTTC" +
            "CTCCACAAATTACCCAGTCTTGATGAATCAGCTCTGTCTAGACATAGGGCAAGGTGAACCCCTTGGGCAGTTACACAACCTCCGCCTTCTGGGTTTAAGCAATTCTCCTGCCTCAGCCTCCGGACTAGCTGGGTCTACAGGTGTGCAGCACCACACCCAGCTA" +
            "GTTATTTGTACTTTTAGTAGAAATGGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAAGTGATCCACCCACCTTGGCCTCCCAAAGTGTTGCGATTACAGGCACTTGCCAGTGAACCTGGCCCTAAATGACTTCTTTCTATCTCCTATAA" +
            "CAATTTGAAATTACTTAAAGGTGGTTTCAAATTGAAAAAATAAAAAGAATTTGGATAAAAATAAAATGTAAACAGTTTTTAAAAATTACAAGAGATTACAAAATATACATGTAAAACCGGAGTGGTCAAAAATGACAAATTTGATTTATTTATAAGGTTTATT" +
            "AAAATTAGCTTTAGTATTGATAATACACTATTACCAAAGTAAAAGCTGATTTTCTCTTGAAAAAAATTTTATGTATTATTAATATGACAGCAAAATACTTCTGTTCACCTTTTGAATATATTCAAAAAGAGAGAGAGTAAAAAAAAAGATAGAATTTTCCCAT" +
            "GATCTGGGTTGGGCCTGGCTTAGCTCAGGGAGGAAGCCCTGCCTGAAAAACGCTGCAGCTTAGGCTGTGACTTTTTCTTCACTCAGTCCTGATCACATTTTCTGTTACTCAGGGCCTGAGGGGGCGGGGGCCTTAAGCATTATCCAATCAGAAACGCTGGGCT" +
            "GACGCCCCGTCCGGGAGGGAGGTGGGGGGGGTCAGCCCCCCGCCCAGCCAGCCGCCCCGTCCGGGAGGTGGGGGGTGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCACCCCGTCTGGGAGGTGTACCCAACAGCTCATT" +
            "GAGAACGGGCCATGATGACAATGGCAGTTTTGTGGAATAGAAACGGGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCTGTGTCTGTGTAGAAAGAGGTAGACATGGGAGACTTCATTTTGTTCTGTACTAAGAAAAATTCTTCTGCCTTGGGATC" +
            "CTGTTGACCTATGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCACTCCCTAATCTCAAGTA" +
            "CCCGGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTGTTGTCCTATGACCCTGCCAAATCCCCCTCTGCGAGAAACACCCAAGAATGATCAACTAAAAAAAA" +
            "AAAAGAAAAAAGAAAAAAGAAAAAATACACCCTGGCAGAAAACTTCTCGAATGTAAGTAAAGTGGGAGTTTTCAACTTGCCCACCTCCTTACAAAAACATGCAATAAGAATGTGCTAGATATAAAGCTTACAAATGTGAGAAATACGGGAAGTTTTTCAATTT" +
            "TAACAAATACTTTCAAAGTCATGTGAAAATTCTCACTAGAAAGAAATCCAATAAAGGCAAGTTGAGGCATAAAAAAAAAAAAAAGAAAAAGCTCAGAATACATATAAGATTGACTAAGTCAGCCAGAAAATAATCCCCTAAAAGAAATTTCTCTCTAAACACC" +
            "CAATGTGCACAGCTACTCTCAGCATGAGAAACATGAGCTTTATGAAGAAAGGTGGCAGATTTTCAGAAGAATTTTATAAAAGTTTCTTTTCCATCTCTGCTGCTCTCTCATCTCCTAGCCATTGAATGGGGGTTCTATATTGAAATACATCTGACAACTTCCA" +
            "ACAACACTTTTTGATCAAGAAATAGAATTTGACTATGTTCGTATAGTGGAATATATTAGAACTTGTAACACAGCTAACTGAATAGCTATTATGGTGTTTGGGTGGCCACATCACCTGTCTTTATTTGTCCGGTAATAGCAGCATTCCAATTTAAAGAAATAAA" +
            "AGATACCAAAATTGTGTTTACTTTTAATTATTCCTATTGAATAAAGTAATAAGCATGTCAGACTGATATCTATCATAACAATAAATTTTGTTTGGATATTATATTAGATATAAATATTTAAGTATGAATAATTTTAATGAACTAGTCATAATGTATGTAGCAT" +
            "TTTTAAAAATTGTAACTATACTTCAGTTAAAACACTTTATATTTCAAAAGCATAAATAACAATATTAAAATAACAATTTAGGTGATTCATTCAAAGTAAGTATTGGGGCTTTATGTTCATACTATTGTAGAAAATACTGCTTATGGCTCATGCCTGTAATCCC" +
            "AGCACATTGGGAGGCTGAGGTGGGTAGATCACCTGAGGTCAGGAGTTCCTGATCCCATCTTTACTAAAAATACAAAACTTACCCAGGGTGGTTGTGCACACTTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACAAGGGAGGAAATGG" +
            "TTGCAGTGAGCCATGATCATGCCACTGAACCCCAGCCTGGGCAAGAGAGTGAGACTGTCTCAAAAAAAAAAAAAACTGTTTAATTTTTATGAATGCAGGTTTTCTGCAAACACTACACATAACTATGCTAATTGTTCTGAAGTAATAAATAGAAAGCAAGGCA" +
            "CAACTACAGACTCCACTGTTCAGTTTATGCACTGAACTGTTCTTGCTTTTGCAGTGTAAGTATTTCTGCCTGCAAATACTGGATAATTACCTTGGATCATCAGATTTCTATCAAAGGAATTTAGTATCTTTTAGTCTTTATCATTTTGTATTGCTAAATTTAT" +
            "CTGTGTGTTAAGCTTCTGTGTGCTCTTAAAATGAGGTTTTATCTAAACAAACCTGTGTCTACTTTAAAAGACTAAACATGAAAAAACTAAACTTTTCAGAACCAAAAACAAAGCAATAAATCTGAAGTACTAGATAGTCTGGAGTGAGATTTATTTAGCTTTT" +
            "TTTTTTTTTTGAGATGGAGTCTCGCTCTGTCACCGAGGCTGGAGTGCAGTGGCACGAACTCGGCTCACTGCAAAAGCTCTGCCTCCCAGCTTCATGCCATTCTCCTACCTCAGCCTCCCAAGTAGCTGGGATTACAGGCAACTGCCACCACGCCCAGCTAATT" +
            "TTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTTGATCTCCTGACCTCATGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCATGCAAAACCGCGCCCTGCCCTTATTTAGCTCTTAATATGATTTACATATAT" +
            "TTCAAAAAAGCAGAGAAAAATATCTACATATAATCTAAATCCCTTAAGAGAGAATAGCAAAAAAATTTTTGGAACTTTTTAAAGAGTTTTGAACTCTTGGACATCTGAATTTTGCACACTATATGCACTTGAAAGAATGTTAATGGGGAAAAAGCAGAAGAGA" +
            "AAAAGACGTTATAAAAAAATCCATGAGTGCACAACACCTATGAAACAGAATAGAGAGCCCAGAAATAATGCCTTTCACCTGCAACCATCAGATTTCTGACAAAGCTGACAAGAGGAATGTGGGAAGAATTCTCTCTTTCATAAATGGTGCTGGAATAACTATC" +
            "TACCACTATGTAGAAGACTGAAGTGGACCCCTTCATTACACCATATAAAAAAATCAACTGAAGATAAATTAAGGACTTAAATGTAAAACTTAAAATTATAAGAAACCCTGCAAGATAACCTAGGAAATAGCATTCTAGACACAGAAACAGGTAAAGACTTCAT" +
            "GATGAAGCTACCAAAAGCAACTGCAACAGAAGTAAATTGACAAATGGGATGTATTTAAACTTAAGAGCTTCTTCACAGCAAAGGAAACTATCAACAGAGTAAACAGACAAACTAGAGAATAAAAGAATATATTTGTAAATTTTGCCTCTGAAAAAGGTCTAAT" +
            "ATACAGAATTTATTAGGAACTTAAACAAGTTTACAAGAAAAAAACACACTCATTAAAAAGTATGCAAAAAACATGAACAGATGCCTTTCAATAGAAGATACACATAGGGCTAACAAGCATATGAAAAAAAAAATGCTCATTGCTAATCATTAGAGAAATTAGA" +
            "AGAGAAATCACAGGCTGGGTGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGCAGATTACAAGGTCAGGAGATCAAGATCATCCTGGCTAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAAATCAGCCAGGTGTGGCAGTGT" +
            "GCACCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATTGCTTGAATCTGGTAGGCAGAGGTTGCAGTGAGCTGAGATCACACCACTGCACTCCTGCCTGGGCAACAGAGCAAGACTCCGTCTCAAACACACACACACAGACACACACACACACACAC" +
            "ACACACACACACACACACACGCAGAGAAACCACAATGAGATACCACCTCATACCAGTCAGAATGGCTATTTTTAAAAAGTCAAAAGATAACAGATGCTGACAAAGTTGCAGAGAAAAGGGAATGCTTATACTCCTCTGGTGGGAGTGTAAATTATTTCAACAA" +
            "CTGTAAAAATCAGTGTGACAATTCCTCACAGAATGAAAAACAGAATTATCATTCGACTGGGAAACCTCATAATTGAGTATATACCCAAAGAAATATGAAATATTATAAAGACACATCCACATGCATGTTCACTGCAGCACTATTCACAATAGCAAAGACACGG" +
            "ACAGACTAAATGCCTATCAATGGCAGACTGGATCAAGAAAATATGGTATGGTCAGATGCGGTGGCTCATGCCTGTAATTCCAGCCCTTTGGGAGGCTGAGGCAGGTGGATTGCCTGAGCTTAGAAGTTTGAGACCACTCTGGGCAACATGGCAAAATTTTGTC" +
            "TCCACAGAAGATACAAAAAAAAAAAAAAAAAA";

    static final boolean testDataInitialized;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInversionWithNovelInsertion;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInversionWithHom_leftPlus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInversionWithHom_leftMinus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInversionWithHom_rightPlus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInversionWithHom_rightMinus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleDeletion_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleDeletion_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInsertion_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleInsertion_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forLongRangeSubstitution_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forLongRangeSubstitution_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forDeletionWithHomology_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forDeletionWithHomology_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleTanDupContraction_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleTanDupContraction_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleTanDupExpansion_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleTanDupExpansion_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleTanDupExpansionWithNovelIns_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forSimpleTanDupExpansionWithNovelIns_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_1to2_pseudoHom_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_1to2_pseudoHom_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_2to1_pseudoHom_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_2to1_pseudoHom_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_3to2_noPseudoHom_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_3to2_noPseudoHom_minus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_2to3_noPseudoHom_plus;
    static final Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations> forComplexTanDup_2to3_noPseudoHom_minus;

    static {
        try{

            final ByteArrayOutputStream outputStream = new ByteArrayOutputStream();

            forSimpleInversionWithNovelInsertion = forSimpleInversionWithNovelInsertion_leftFlankingForwardStrandOnly();
            forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint = forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint();
            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> inversion4 = forSimpleInversionWithHomology(outputStream);
            forSimpleInversionWithHom_leftPlus = inversion4.get(0);
            forSimpleInversionWithHom_leftMinus = inversion4.get(1);
            forSimpleInversionWithHom_rightPlus = inversion4.get(2);
            forSimpleInversionWithHom_rightMinus = inversion4.get(3);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> simpleDeletion = forSimpleDeletion(outputStream);
            forSimpleDeletion_plus = simpleDeletion.get(0);
            forSimpleDeletion_minus = simpleDeletion.get(1);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> simpleInsertion = forSimpleInsertion(outputStream);
            forSimpleInsertion_plus = simpleInsertion.get(0);
            forSimpleInsertion_minus = simpleInsertion.get(1);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> longRangeSubstitution = forLongRangeSubstitution();
            forLongRangeSubstitution_plus = longRangeSubstitution.get(0);
            forLongRangeSubstitution_minus = longRangeSubstitution.get(1);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> deletionWithHomology = forDeletionWithHomology(outputStream);
            forDeletionWithHomology_plus = deletionWithHomology.get(0);
            forDeletionWithHomology_minus = deletionWithHomology.get(1);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> simpleTandemDuplicationContraction = forSimpleTandemDuplicationContraction();
            forSimpleTanDupContraction_plus = simpleTandemDuplicationContraction.get(0);
            forSimpleTanDupContraction_minus = simpleTandemDuplicationContraction.get(1);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> simpleTandemDuplicationExpansion = forSimpleTandemDuplicationExpansion(outputStream);
            forSimpleTanDupExpansion_plus = simpleTandemDuplicationExpansion.get(0);
            forSimpleTanDupExpansion_minus = simpleTandemDuplicationExpansion.get(1);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> simpleTandemDuplicationExpansionWithNovelInsertion = forSimpleTandemDuplicationExpansionWithNovelInsertion(outputStream);
            forSimpleTanDupExpansionWithNovelIns_plus = simpleTandemDuplicationExpansionWithNovelInsertion.get(0);
            forSimpleTanDupExpansionWithNovelIns_minus = simpleTandemDuplicationExpansionWithNovelInsertion.get(1);

            final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> complexTandemDuplication = forComplexTandemDuplication();
            forComplexTanDup_1to2_pseudoHom_plus = complexTandemDuplication.get(0);
            forComplexTanDup_1to2_pseudoHom_minus = complexTandemDuplication.get(1);
            forComplexTanDup_2to1_pseudoHom_plus = complexTandemDuplication.get(2);
            forComplexTanDup_2to1_pseudoHom_minus = complexTandemDuplication.get(3);
            forComplexTanDup_3to2_noPseudoHom_plus = complexTandemDuplication.get(4);
            forComplexTanDup_3to2_noPseudoHom_minus = complexTandemDuplication.get(5);
            forComplexTanDup_2to3_noPseudoHom_plus = complexTandemDuplication.get(6);
            forComplexTanDup_2to3_noPseudoHom_minus = complexTandemDuplication.get(7);

            outputStream.close();

            testDataInitialized = true;
        } catch (final IOException ioex) {
            throw new GATKException("Failed to create test data " + SVCallerTestDataProvider.class);
        }
    }

    private static Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>
    forSimpleInversionWithNovelInsertion_leftFlankingForwardStrandOnly() throws IOException {
        // inversion with inserted sequence
        final byte[] leftFlank = makeDummySequence(146, (byte)'A');
        final byte[] rightFlankRC = makeDummySequence(50, (byte)'C');
        final byte[] contigSeq = new byte[leftFlank.length+1+rightFlankRC.length];
        System.arraycopy(leftFlank, 0, contigSeq, 0, leftFlank.length);
        contigSeq[leftFlank.length] = (byte) 'T';
        System.arraycopy(rightFlankRC, 0, contigSeq, leftFlank.length+1, rightFlankRC.length);

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 108569149, 108569294), TextCigarCodec.decode("146M51S"), true, 60, 0, 1, 146);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 108569315, 108569364), TextCigarCodec.decode("147S50M"), false, 60, 0, 148, 197);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        return new Tuple3<>(region1, region2, breakpoints);
    }

    private static Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>
    forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint() throws IOException {
        // inversion with strange left breakpoint
        final byte[] contigSequence = LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", new SimpleInterval(chrForLongContig1, 20138007, 20142231), TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 1, contigSequence.length - 1986);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", new SimpleInterval(chrForLongContig1, 20152030, 20154634), TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 3604, contigSequence.length);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(ChimericAlignment.fromSplitAlignments(new Tuple2<>(Arrays.asList(region1, region2), contigSequence)).get(0));
        return new Tuple3<>(region1, region2, breakpoints);
    }

    /**
     * The following four tests are all going to be for the same inversion, testing if implementations are correct for
     * identifying the breakpoints by looking at different representations of evidence.
     * The inversion we are looking at is
     *
     * '+' strand representation: G....100....G|ACACA|C....100....C               A....100....A|TGTGT|T....100....T
     *
     * 100-bases of 'G' is the left flanking before the homologyForwardStrandRep |ACACA| and the region starting with 100-bases of 'C' and
     * ending with 100-bases of 'A' and maybe (homologyForwardStrandRep uncertainty) the homologyForwardStrandRep |TGTGT| is inverted.
     * 100-bases of 'T' is the right flanking region.
     *
     * Returns a list of four Tuple3's with left flanking evidence '+'/'-' strand representation and right flanking side.
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forSimpleInversionWithHomology(final ByteArrayOutputStream outputStream) throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();

        final byte[] leftLeftPlus = makeDummySequence(100, (byte)'G');
        final byte[] leftLeftMinus = makeDummySequence(100, (byte)'C');
        final byte[] leftRightPlus = makeDummySequence(100, (byte)'C');
        final byte[] leftRightMinus = makeDummySequence(100, (byte)'G');
        final byte[] rightLeftPlus = makeDummySequence(100, (byte)'A');
        final byte[] rightLeftMinus = makeDummySequence(100, (byte)'T');
        final byte[] rightRightPlus = makeDummySequence(100, (byte)'T');
        final byte[] rightRightMinus = makeDummySequence(100, (byte)'A');
        final byte[] leftHomology = "ACACA".getBytes();
        final byte[] rightHomology = "TGTGT".getBytes();
        {// left flanking evidence '+'/'-' strand representation
            outputStream.reset();
            outputStream.write(leftLeftPlus);outputStream.write(leftHomology);outputStream.write(rightLeftMinus);
            byte[] contigSeq = outputStream.toByteArray();
            AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 101, 205), TextCigarCodec.decode("105M100S"), true, 60, 0, 1, 105);
            AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 501, 605), TextCigarCodec.decode("100S105M"), false, 60, 0, 101, 205);
            final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, new ArrayList<>()));
            result.add(new Tuple3<>(region1, region2, breakpoints));

            outputStream.reset();
            outputStream.write(rightLeftPlus);outputStream.write(rightHomology);outputStream.write(leftLeftMinus);
            contigSeq = outputStream.toByteArray();
            region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 501, 605), TextCigarCodec.decode("105M100S"), true, 60, 0, 1, 105);
            region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 101, 205), TextCigarCodec.decode("100S105M"), false, 60, 0, 101, 205);
            final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, new ArrayList<>()));
            result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));
        }
        {// right flanking evidence '+'/'-' strand representation
            outputStream.reset();
            outputStream.write(leftRightMinus);outputStream.write(rightHomology);outputStream.write(rightRightPlus);
            byte[] contigSeq = outputStream.toByteArray();

            AlignmentRegion region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 201, 305), TextCigarCodec.decode("105M100S"), false, 60, 0, 1, 105);
            AlignmentRegion region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 601, 705), TextCigarCodec.decode("100S105M"), true, 60, 0, 101, 205);
            final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, new ArrayList<>()));
            result.add(new Tuple3<>(region1, region2, breakpoints));

            outputStream.reset();
            outputStream.write(rightRightMinus);outputStream.write(leftHomology);outputStream.write(leftRightPlus);
            contigSeq = outputStream.toByteArray();

            region1 = new AlignmentRegion("1","1", new SimpleInterval("20", 601, 705), TextCigarCodec.decode("105M100S"), false, 60, 0, 1, 105);
            region2 = new AlignmentRegion("1","1", new SimpleInterval("20", 201, 305), TextCigarCodec.decode("100S105M"), true, 60, 0, 101, 205);
            final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, new ArrayList<>()));
            result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));
        }
        return result;
    }

    /**
     * 40-'A' + 10-'C'+10-'T' + 40-'G' where the segment 10-'C'+10-'T' is deleted (forward strand representation description).
     *
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forSimpleDeletion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();
        // simple deletion '+' strand representation
        final byte[] leftRefFlank = makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = makeDummySequence(40, (byte)'G');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100040), TextCigarCodec.decode("40M40S"), true, 60, 0, 1 ,40);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100061, 100100), TextCigarCodec.decode("40S40M"), true, 60, 0, 41 ,80);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // simple deletion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100061, 100100), TextCigarCodec.decode("40M40S"), false, 60, 0, 1 ,40);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100040), TextCigarCodec.decode("40S40M"), false, 60, 0, 41 ,80);
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));

        return result;
    }

    /**
     * 100-'A' + 100-'T' and a 50 bases of 'C' is inserted at the A->T junction point (forward strand description)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forSimpleInsertion(final ByteArrayOutputStream outputStream) throws IOException {
        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();

        // simple insertion '+' strand representation
        final byte[] leftRefFlank = makeDummySequence(100, (byte)'A');
        final byte[] insertedSeq  =makeDummySequence(50, (byte)'C');
        final byte[] rightRefFlank = makeDummySequence(100, (byte)'T');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(insertedSeq);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100100), TextCigarCodec.decode("100M100S"), true, 60, 0, 1 ,100);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100101, 100200), TextCigarCodec.decode("100S100M"), true, 60, 0, 151 ,250);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // simple insertion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(insertedSeq);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(insertedSeq);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100101, 100200), TextCigarCodec.decode("100M100S"), false, 60, 0, 1 ,100);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100100), TextCigarCodec.decode("100S100M"), false, 60, 0, 151 ,250);
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));

        return result;
    }

    /**
     * 50-'A' + 50-'C' where the middle 10-'A'+10-'C' is substituted with 10-'G' (forward strand representation)
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forLongRangeSubstitution() throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();

        // long range substitution '+' strand representation
        final byte[] leftRefFlank = makeDummySequence(50, (byte)'A');
        final byte[] rightRefFlank = makeDummySequence(50, (byte)'G');
        final byte[] substitution = makeDummySequence(10, (byte)'C');
        byte[] contigSeq = new byte[leftRefFlank.length+rightRefFlank.length-10];
        System.arraycopy(leftRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(substitution, 0, contigSeq, 40, substitution.length);
        System.arraycopy(rightRefFlank, 0, contigSeq, 50, 40);
        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100040), TextCigarCodec.decode("40M50S"), true, 60, 0, 1 ,40);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100061, 100100), TextCigarCodec.decode("50S40M"), true, 60, 0, 51 ,90);
        NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // long range substitution '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(substitution);
        System.arraycopy(rightRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(substitution, 0, contigSeq, 40, substitution.length);
        System.arraycopy(leftRefFlank, 0, contigSeq, 40 + substitution.length, 40);
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100061, 100100), TextCigarCodec.decode("40M50S"), false, 60, 0, 1 ,40);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100040), TextCigarCodec.decode("50S40M"), false, 60, 0, 51 ,90);
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));

        return result;
    }

    /**
     * 40-'C' + 'ATCG' + 34 bases of unique sequence + 'ATCG' + 40-'T' is shrunk to 40-'C' + 'ATCG' + 40-'T' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forDeletionWithHomology(final ByteArrayOutputStream outputStream) throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();

        // simple deletion with homology '+' strand representation
        final byte[] leftRefFlank = makeDummySequence(40, (byte)'C');
        final byte[] rightRefFlank = makeDummySequence(40, (byte)'T');
        final byte[] homology = new byte[]{'A', 'T', 'C', 'G'};
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(homology);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100044), TextCigarCodec.decode("44M40S"), true, 60, 0, 1 ,44);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100079, 100122), TextCigarCodec.decode("40S44M"), true, 60, 0, 41 ,84);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // simple deletion with homology '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(homology);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(homology);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100079, 100122), TextCigarCodec.decode("44M40S"), false, 60, 0, 1 ,44);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100044), TextCigarCodec.decode("40S44M"), false, 60, 0, 41 ,84);
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));

        return result;
    }

    /**
     * 40-'A' + 20-'C' + 40-'G' is shrunk to 40-'A' + 10-'C' + 40-'G' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forSimpleTandemDuplicationContraction() throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();

        // simple tandem duplication contraction '+' strand representation
        final byte[] leftRefFlank = makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = makeDummySequence(40, (byte)'G');
        final byte[] doubleDup = makeDummySequence(20, (byte)'C');
        final byte[] contigSeq = new byte[90];
        System.arraycopy(leftRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(doubleDup, 0, contigSeq, 40, 10);
        System.arraycopy(rightRefFlank, 0, contigSeq, 50, 40);

        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100050), TextCigarCodec.decode("50M40S"), true, 60, 0, 1 ,50);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100051, 100100), TextCigarCodec.decode("40S50M"), true, 60, 0, 41 ,100);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // simple tandem duplication contraction '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(doubleDup);
        System.arraycopy(rightRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(doubleDup, 0, contigSeq, 40, 10);
        System.arraycopy(leftRefFlank, 0, contigSeq, 50, 40);
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100051, 100100), TextCigarCodec.decode("50M40S"), false, 60, 0, 1 ,50);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100050), TextCigarCodec.decode("40S50M"), false, 60, 0, 41 ,100);
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));

        return result;
    }

    /**
     * 40-'A' + 10-'C' + 40-'G' is expanded to 40-'A' + 20-'C' + 40-'G' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forSimpleTandemDuplicationExpansion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();

        // simple tandem duplication expansion '+' strand representation
        final byte[] leftRefFlank = makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = makeDummySequence(40, (byte)'G');
        final byte[] doubleDup = makeDummySequence(20, (byte)'C');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(doubleDup);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();

        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100050), TextCigarCodec.decode("50M50S"), true, 60, 0, 1 ,50);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100041, 100090), TextCigarCodec.decode("50S50M"), true, 60, 0, 51 ,100);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // simple tandem duplication expansion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(doubleDup);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(doubleDup);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100041, 100090), TextCigarCodec.decode("50M50S"), false, 60, 0, 1 ,50);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 100001, 100050), TextCigarCodec.decode("50S50M"), false, 60, 0, 51 ,100);
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));

        return result;
    }

    /**
     * System.out.println(new String(reference.getReferenceBases(dummyOptions, new SimpleInterval("21", 25297100, 25297300)).getBases()));
     * leftFlank:  chr21:25297101-25297163
     * repeat:     chr21:25297164-25297252
     * rightFlank: chr21:25297253-25297300
     * GTTAGTAGATATTCTAGCTGACTCAGTTCAGTGTTGCTATGATTAAACAAGAGTGAGTTCCCT
     * AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG
     * CATTATTGATATTTCATTATGTTCAACAGATGGAGTTAATGTGAATGT
     *
     * insertedSequenceForwardStrandRep: CTCTCTCTCT
     *
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forSimpleTandemDuplicationExpansionWithNovelInsertion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();
        // simple tandem duplication expansion with novel insertion '+' strand representation
        final byte[] leftRefFlank = "GTTAGTAGATATTCTAGCTGACTCAGTTCAGTGTTGCTATGATTAAACAAGAGTGAGTTCCCT".getBytes();                     //63
        final byte[] rightRefFlank = "CATTATTGATATTTCATTATGTTCAACAGATGGAGTTAATGTGAATGT".getBytes();                                   //48
        final byte[] insertedSeq = "CTCTCTCTCT".getBytes();                                                                           //10
        final byte[] dup = "AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG".getBytes();    //89
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(dup);outputStream.write(insertedSeq);outputStream.write(dup);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();

        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 25297101, 25297252), TextCigarCodec.decode("152M147S"), true, 60, 0, 1 ,152);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 25297164, 25297300), TextCigarCodec.decode("162S137M"), true, 60, 0, 163 ,299);
        final NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // simple tandem duplication expansion with novel insertion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(insertedSeq);
        SequenceUtil.reverseComplement(dup);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(dup);outputStream.write(insertedSeq);outputStream.write(dup);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();

        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 25297164, 25297300), TextCigarCodec.decode("137M162S"), false, 60, 0, 1 ,137);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("21", 25297101, 25297252), TextCigarCodec.decode("147S152M"), false, 60, 0, 148 ,299);
        final NovelAdjacencyReferenceLocations breakpointsDetectedFromReverseStrand = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeq, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpointsDetectedFromReverseStrand));

        return result;
    }

    /**
     * These test data was based on a real observation on a locally-assembled contig
     * "TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC"
     * with two alignment records chr18:312579-312718 140M135S
     *                            chr18:312610-312757 127S148M
     * for a tandem repeat expansion event from 1 copy to 2 copies with also a pseudo-homologyForwardStrandRep

     * Return a list of eight entries for positive and reverse strand representations for:
     * 1. expansion from 1 unit to 2 units with pseudo-homology
     * 2. contraction from 2 units to 1 unit with pseudo-homology
     * 3. contraction from 3 units to 2 units without pseudo-homology
     * 4. expansion from 2 units to 3 units without pseudo-homology
     */
    private static List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>>
    forComplexTandemDuplication() throws IOException {

        final List<Tuple3<AlignmentRegion, AlignmentRegion, NovelAdjacencyReferenceLocations>> result = new ArrayList<>();
        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";                                                                    // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";                                                            // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA";   // 96
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA";   // 96
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                                                                      // 13


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        final byte[] fakeRefSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, pseudoHomology, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, pseudoHomology, rightRefFlank).getBytes();
        AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312718), TextCigarCodec.decode("140M135S"), true, 60, 0, 1 ,140);
        AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312610, 312757), TextCigarCodec.decode("127S148M"), true, 60, 0, 128 ,275);
        NovelAdjacencyReferenceLocations breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexExpansionWithPseudoHomology, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        final byte[] fakeRefSeqForComplexExpansionWithPseudoHomology_reverseStrand = Arrays.copyOf(fakeRefSeqForComplexExpansionWithPseudoHomology, fakeRefSeqForComplexExpansionWithPseudoHomology.length);
        final byte[] contigSeqForComplexExpansionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionWithPseudoHomology, contigSeqForComplexExpansionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(fakeRefSeqForComplexExpansionWithPseudoHomology_reverseStrand);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionWithPseudoHomology_reverseStrand);

        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312610, 312757), TextCigarCodec.decode("148M127S"), false, 60, 0, 1 ,148);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312718), TextCigarCodec.decode("135S140M"), false, 60, 0, 136 ,275);
        breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexExpansionWithPseudoHomology_reverseStrand, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        final byte[] fakeRefSeqForComplexContractionWithPseudoHomology = contigSeqForComplexExpansionWithPseudoHomology;
        final byte[] contigSeqForComplexContractionWithPseudoHomology = fakeRefSeqForComplexExpansionWithPseudoHomology;
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312718), TextCigarCodec.decode("140M39S"), true, 60, 0, 1, 140);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312706, 312853), TextCigarCodec.decode("31S148M"), true, 60, 0, 32, 179);
        breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexContractionWithPseudoHomology, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        final byte[] fakeRefSeqForComplexContractionWithPseudoHomology_reverseStrand = Arrays.copyOf(fakeRefSeqForComplexContractionWithPseudoHomology, fakeRefSeqForComplexContractionWithPseudoHomology.length);
        final byte[] contigSeqForComplexContractionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionWithPseudoHomology, contigSeqForComplexContractionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(fakeRefSeqForComplexContractionWithPseudoHomology_reverseStrand);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionWithPseudoHomology_reverseStrand);
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312706, 312853), TextCigarCodec.decode("148M31S"), false, 60, 0, 1, 148);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312718), TextCigarCodec.decode("39S140M"), false, 60, 0, 40, 179);
        breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexContractionWithPseudoHomology_reverseStrand, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // third test: contraction from 3 units to 2 units without pseudo-homology
        final byte[] fakeRefSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, firstRepeat, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, rightRefFlank).getBytes();

        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312801), TextCigarCodec.decode("223M39S"), true, 60, 0, 1, 223);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312706, 312936), TextCigarCodec.decode("31S231M"), true, 60, 0, 32, 262);
        breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexContractionNoPseudoHomology, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        final byte[] fakeRefSeqForComplexContractionNoPseudoHomology_reverseStrand = Arrays.copyOf(fakeRefSeqForComplexContractionNoPseudoHomology, fakeRefSeqForComplexContractionNoPseudoHomology.length);
        final byte[] contigSeqForComplexContractionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionNoPseudoHomology, contigSeqForComplexContractionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(fakeRefSeqForComplexContractionNoPseudoHomology_reverseStrand);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionNoPseudoHomology_reverseStrand);
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312706, 312936), TextCigarCodec.decode("231M31S"), false, 60, 0, 1, 231);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312801), TextCigarCodec.decode("39S223M"), false, 60, 0, 40, 262);
        breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexContractionNoPseudoHomology_reverseStrand, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        final byte[] fakeRefSeqForComplexExpansionNoPseudoHomology = contigSeqForComplexContractionNoPseudoHomology;
        final byte[] contigSeqForComplexExpansionNoPseudoHomology = fakeRefSeqForComplexContractionNoPseudoHomology;
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312801), TextCigarCodec.decode("223M135S"), true, 60, 0, 1, 223);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312610, 312840), TextCigarCodec.decode("127S231M"), true, 60, 0, 128, 358);
        breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexExpansionNoPseudoHomology, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        final byte[] fakeRefSeqForComplexExpansionNoPseudoHomology_reverseStrand = Arrays.copyOf(fakeRefSeqForComplexExpansionNoPseudoHomology, fakeRefSeqForComplexExpansionNoPseudoHomology.length);
        final byte[] contigSeqForComplexExpansionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionNoPseudoHomology, contigSeqForComplexExpansionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(fakeRefSeqForComplexExpansionNoPseudoHomology_reverseStrand);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionNoPseudoHomology_reverseStrand);
        region1 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312610, 312840), TextCigarCodec.decode("231M127S"), false, 60, 0, 1, 231);
        region2 = new AlignmentRegion("1", "contig-1", new SimpleInterval("20", 312579, 312801), TextCigarCodec.decode("135S223M"), false, 60, 0, 136, 358);
        breakpoints = new NovelAdjacencyReferenceLocations(new ChimericAlignment(region1, region2, contigSeqForComplexExpansionNoPseudoHomology_reverseStrand, Collections.emptyList()));
        result.add(new Tuple3<>(region1, region2, breakpoints));

        return result;
    }

    static byte[] makeDummySequence(final int length, byte base) {
        final byte[] result = new byte[length];
        Arrays.fill(result, base);
        return result;
    }
}
