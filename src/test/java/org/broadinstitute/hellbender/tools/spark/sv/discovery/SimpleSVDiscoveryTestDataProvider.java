package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryUtils.getCanonicalChromosomes;

/**
 * Provides test data for testing several methods involved in the SV variant caller,
 * but specifically on simple types.
 * NO TESTS ARE RUN IN THIS PARTICULAR CLASS.
 */
public final class SimpleSVDiscoveryTestDataProvider {

    public static final class TestDataForSimpleSVs {
        public final AlignmentInterval firstAlignment;
        public final AlignmentInterval secondAlignment;
        public final NovelAdjacencyAndAltHaplotype biPathBubble;
        public final String evidenceAssemblyContigName;

        TestDataForSimpleSVs(final AlignmentInterval firstAlignment, final AlignmentInterval secondAlignment,
                             final NovelAdjacencyAndAltHaplotype biPathBubble, final String evidenceAssemblyContigName) {
            this.firstAlignment = firstAlignment;
            this.secondAlignment = secondAlignment;
            this.biPathBubble = biPathBubble;
            this.evidenceAssemblyContigName = evidenceAssemblyContigName;
        }
    }

    public static final ReferenceMultiSource b37_reference = new ReferenceMultiSource(
            GATKBaseTest.b37_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b37_seqDict = b37_reference.getReferenceSequenceDictionary(null);
    public static final Set<String> canonicalChromosomes = getCanonicalChromosomes(null, b37_seqDict);

    public static final ReferenceMultiSource b38_reference = new ReferenceMultiSource(
            GATKBaseTest.b38_reference_20_21, ReferenceWindowFunctions.IDENTITY_FUNCTION);
    public static final SAMSequenceDictionary b38_seqDict = b38_reference.getReferenceSequenceDictionary(null);
    public static final Set<String> b38_canonicalChromosomes = getCanonicalChromosomes(null, b38_seqDict);

    // the chromosome that the long contig1 is supposed to be mapped to is actually chr19, but to make tests runnable, we could only use "20" or "21"
    // todo: this should be fixed, but since the exact mapped to chromosome is not important now, we push it to later
    public static final String chrForLongContig1 = "20";
    public static final String homologyForLongContig1 = "AAAAAAGCTAAATAAATCTCACTCCAGACTATCTAGTACTTCAGATTTATTGCTTTGTTTTTGGTTCTGAAAAGTTTAGTTTTTTCATGTTTAGTCTTTTAAAGTAGACACAGGTTTGTTTAGATAAAACCTCATTTTAAGAGCACACAGAAGCTTAACACACAGATAAATTTAGCAATACAAAATGATAAAGACTAAAAGATACTAAATTCCTTTGATAGAAATCTGATGATCCAAGGTAATTATCCAGTATTTGCAGGCAGAAATACTTACACTGCAAAAGCAAGAACAGTTCAGTGCATAAACTGAACAGTGGAGTCTGTAGTTGTGCCTTGCTTTCTATTTATTACTTCAGAACAATTAGCATAGTTATGTGTAGTGTTTGCAGAAAACCTGCATTCATAAAAATTAAACAGTTTTTTTTTTTTTTGAGACAGTCTCACTCTCTTGCCCAGGCTGGGGTTCAGTGGCATGATCATGGCTCACTGCAACCATTTCCTCCCTTGTTCAAGCAATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTACAAGTGTGCACAACCACCCTGGGTAAGTTTTGTATTTTTAGTAAAGATGGGATCAGGAACTCCTGACCTCA";
    public static final String LONG_CONTIG1 =
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

    public static final boolean testDataInitialized;
    public static final TestDataForSimpleSVs forSimpleInversionWithNovelInsertion;
    public static final TestDataForSimpleSVs forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
    public static final TestDataForSimpleSVs forSimpleInversionWithHom_leftPlus;
    public static final TestDataForSimpleSVs forSimpleInversionWithHom_leftMinus;
    public static final TestDataForSimpleSVs forSimpleInversionWithHom_rightPlus;
    public static final TestDataForSimpleSVs forSimpleInversionWithHom_rightMinus;
    public static final TestDataForSimpleSVs forSimpleDeletion_plus;
    public static final TestDataForSimpleSVs forSimpleDeletion_minus;
    public static final TestDataForSimpleSVs forSimpleInsertion_plus;
    public static final TestDataForSimpleSVs forSimpleInsertion_minus;
    public static final TestDataForSimpleSVs forLongRangeSubstitution_plus;
    public static final TestDataForSimpleSVs forLongRangeSubstitution_minus;
    public static final TestDataForSimpleSVs forDeletionWithHomology_plus;
    public static final TestDataForSimpleSVs forDeletionWithHomology_minus;
    public static final TestDataForSimpleSVs forSimpleTanDupContraction_plus;
    public static final TestDataForSimpleSVs forSimpleTanDupContraction_minus;
    public static final TestDataForSimpleSVs forSimpleTanDupExpansion_plus;
    public static final TestDataForSimpleSVs forSimpleTanDupExpansion_minus;
    public static final TestDataForSimpleSVs forSimpleTanDupExpansionWithNovelIns_plus;
    public static final TestDataForSimpleSVs forSimpleTanDupExpansionWithNovelIns_minus;
    public static final TestDataForSimpleSVs forComplexTanDup_1to2_pseudoHom_plus;
    public static final TestDataForSimpleSVs forComplexTanDup_1to2_pseudoHom_minus;
    public static final TestDataForSimpleSVs forComplexTanDup_2to1_pseudoHom_plus;
    public static final TestDataForSimpleSVs forComplexTanDup_2to1_pseudoHom_minus;
    public static final TestDataForSimpleSVs forComplexTanDup_3to2_noPseudoHom_plus;
    public static final TestDataForSimpleSVs forComplexTanDup_3to2_noPseudoHom_minus;
    public static final TestDataForSimpleSVs forComplexTanDup_2to3_noPseudoHom_plus;
    public static final TestDataForSimpleSVs forComplexTanDup_2to3_noPseudoHom_minus;

    static {
        try ( final ByteArrayOutputStream outputStream = new ByteArrayOutputStream() ){

            forSimpleInversionWithNovelInsertion = forSimpleInversionWithNovelInsertion_leftFlankingForwardStrandOnly();
            forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint = forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint();
            final List<TestDataForSimpleSVs> inversion4 = forSimpleInversionWithHomology(outputStream);
            forSimpleInversionWithHom_leftPlus = inversion4.get(0);
            forSimpleInversionWithHom_leftMinus = inversion4.get(1);
            forSimpleInversionWithHom_rightPlus = inversion4.get(2);
            forSimpleInversionWithHom_rightMinus = inversion4.get(3);

            final List<TestDataForSimpleSVs> simpleDeletion = forSimpleDeletion(outputStream);
            forSimpleDeletion_plus = simpleDeletion.get(0);
            forSimpleDeletion_minus = simpleDeletion.get(1);

            final List<TestDataForSimpleSVs> simpleInsertion = forSimpleInsertion(outputStream);
            forSimpleInsertion_plus = simpleInsertion.get(0);
            forSimpleInsertion_minus = simpleInsertion.get(1);

            final List<TestDataForSimpleSVs> longRangeSubstitution = forLongRangeSubstitution();
            forLongRangeSubstitution_plus = longRangeSubstitution.get(0);
            forLongRangeSubstitution_minus = longRangeSubstitution.get(1);

            final List<TestDataForSimpleSVs> deletionWithHomology = forDeletionWithHomology(outputStream);
            forDeletionWithHomology_plus = deletionWithHomology.get(0);
            forDeletionWithHomology_minus = deletionWithHomology.get(1);

            final List<TestDataForSimpleSVs> simpleTandemDuplicationContraction = forSimpleTandemDuplicationContraction();
            forSimpleTanDupContraction_plus = simpleTandemDuplicationContraction.get(0);
            forSimpleTanDupContraction_minus = simpleTandemDuplicationContraction.get(1);

            final List<TestDataForSimpleSVs> simpleTandemDuplicationExpansion = forSimpleTandemDuplicationExpansion(outputStream);
            forSimpleTanDupExpansion_plus = simpleTandemDuplicationExpansion.get(0);
            forSimpleTanDupExpansion_minus = simpleTandemDuplicationExpansion.get(1);

            final List<TestDataForSimpleSVs> simpleTandemDuplicationExpansionWithNovelInsertion = forSimpleTandemDuplicationExpansionWithNovelInsertion(outputStream);
            forSimpleTanDupExpansionWithNovelIns_plus = simpleTandemDuplicationExpansionWithNovelInsertion.get(0);
            forSimpleTanDupExpansionWithNovelIns_minus = simpleTandemDuplicationExpansionWithNovelInsertion.get(1);

            final List<TestDataForSimpleSVs> complexTandemDuplication = forComplexTandemDuplication();
            forComplexTanDup_1to2_pseudoHom_plus = complexTandemDuplication.get(0);
            forComplexTanDup_1to2_pseudoHom_minus = complexTandemDuplication.get(1);
            forComplexTanDup_2to1_pseudoHom_plus = complexTandemDuplication.get(2);
            forComplexTanDup_2to1_pseudoHom_minus = complexTandemDuplication.get(3);
            forComplexTanDup_3to2_noPseudoHom_plus = complexTandemDuplication.get(4);
            forComplexTanDup_3to2_noPseudoHom_minus = complexTandemDuplication.get(5);
            forComplexTanDup_2to3_noPseudoHom_plus = complexTandemDuplication.get(6);
            forComplexTanDup_2to3_noPseudoHom_minus = complexTandemDuplication.get(7);

            testDataInitialized = true;
        } catch (final Exception ioex) {
            throw new GATKException("Failed to create test data ", ioex);
        }
    }

    public static List<TestDataForSimpleSVs> getAllTestData() {
        final List<TestDataForSimpleSVs> testDataForSimpleSVs = Arrays.asList(forSimpleInversionWithNovelInsertion,
                forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint,
                forSimpleInversionWithHom_leftPlus,
                forSimpleInversionWithHom_leftMinus,
                forSimpleInversionWithHom_rightPlus,
                forSimpleInversionWithHom_rightMinus,
                forSimpleDeletion_plus,
                forSimpleDeletion_minus,
                forSimpleInsertion_plus,
                forSimpleInsertion_minus,
                forLongRangeSubstitution_plus,
                forLongRangeSubstitution_minus,
                forDeletionWithHomology_plus,
                forDeletionWithHomology_minus,
                forSimpleTanDupContraction_plus,
                forSimpleTanDupContraction_minus,
                forSimpleTanDupExpansion_plus,
                forSimpleTanDupExpansion_minus,
                forSimpleTanDupExpansionWithNovelIns_plus,
                forSimpleTanDupExpansionWithNovelIns_minus,
                forComplexTanDup_1to2_pseudoHom_plus,
                forComplexTanDup_1to2_pseudoHom_minus,
                forComplexTanDup_2to1_pseudoHom_plus,
                forComplexTanDup_2to1_pseudoHom_minus,
                forComplexTanDup_3to2_noPseudoHom_plus,
                forComplexTanDup_3to2_noPseudoHom_minus,
                forComplexTanDup_2to3_noPseudoHom_plus,
                forComplexTanDup_2to3_noPseudoHom_minus);
        return Collections.unmodifiableList(testDataForSimpleSVs);
    }

    // same event, two representations
    public static List<Tuple2<TestDataForSimpleSVs, TestDataForSimpleSVs>> getAllTestDataPaired() {
        final List<Tuple2<TestDataForSimpleSVs, TestDataForSimpleSVs>> testDataForSimpleSVs =
                Arrays.asList(
                        new Tuple2<>(forSimpleInversionWithHom_leftPlus, forSimpleInversionWithHom_leftMinus),
                        new Tuple2<>(forSimpleInversionWithHom_rightPlus, forSimpleInversionWithHom_rightMinus),
                        new Tuple2<>(forSimpleDeletion_plus, forSimpleDeletion_minus),
                        new Tuple2<>(forSimpleInsertion_plus, forSimpleInsertion_minus),
                        new Tuple2<>(forLongRangeSubstitution_plus, forLongRangeSubstitution_minus),
                        new Tuple2<>(forDeletionWithHomology_plus, forDeletionWithHomology_minus),
                        new Tuple2<>(forSimpleTanDupContraction_plus, forSimpleTanDupContraction_minus),
                        new Tuple2<>(forSimpleTanDupExpansion_plus, forSimpleTanDupExpansion_minus),
                        new Tuple2<>(forSimpleTanDupExpansionWithNovelIns_plus, forSimpleTanDupExpansionWithNovelIns_minus),
                        new Tuple2<>(forComplexTanDup_1to2_pseudoHom_plus, forComplexTanDup_1to2_pseudoHom_minus),
                        new Tuple2<>(forComplexTanDup_2to1_pseudoHom_plus, forComplexTanDup_2to1_pseudoHom_minus),
                        new Tuple2<>(forComplexTanDup_3to2_noPseudoHom_plus, forComplexTanDup_3to2_noPseudoHom_minus),
                        new Tuple2<>(forComplexTanDup_2to3_noPseudoHom_plus, forComplexTanDup_2to3_noPseudoHom_minus));
        return Collections.unmodifiableList(testDataForSimpleSVs);
    }

    private static TestDataForSimpleSVs
    forSimpleInversionWithNovelInsertion_leftFlankingForwardStrandOnly() {
        // inversion with inserted sequence
        final byte[] leftFlank = SVTestUtils.makeDummySequence(146, (byte)'A');
        final byte[] rightFlankRC = SVTestUtils.makeDummySequence(50, (byte)'C');
        final byte[] contigSeq = new byte[leftFlank.length+1+rightFlankRC.length];
        System.arraycopy(leftFlank, 0, contigSeq, 0, leftFlank.length);
        contigSeq[leftFlank.length] = (byte) 'T';
        System.arraycopy(rightFlankRC, 0, contigSeq, leftFlank.length+1, rightFlankRC.length);

        // reference intervals are changed from real values, which were generated on a different chromosome, to a fake but valid value.
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 69149, 69294), 1, 146, TextCigarCodec.decode("146M51S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 69315, 69364), 148, 197, TextCigarCodec.decode("147S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignedContig alignedContig = new AlignedContig("asm000001:tig00001", contigSeq, Arrays.asList(region1, region2));
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), alignedContig.getContigName(), b37_seqDict), alignedContig.getContigSequence(), b37_seqDict);
        return new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001");
    }

    private static TestDataForSimpleSVs
    forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint() {
        // inversion with strange left breakpoint
        final byte[] contigSequence = LONG_CONTIG1.getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval(chrForLongContig1, 20138007, 20142231), 1, contigSequence.length - 1986, TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval(chrForLongContig1, 20152030, 20154634), 3604, contigSequence.length, TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final AlignedContig alignedContig = new AlignedContig("asm702700:tig00001", contigSequence, Arrays.asList(region1, region2));
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(DiscoverVariantsFromContigAlignmentsSAMSpark.parseOneContig(alignedContig, b37_seqDict, true, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, true).get(0), contigSequence, b37_seqDict);
        return new TestDataForSimpleSVs(region1, region2, breakpoints, "asm702700:tig00001");
    }

    /**
     * The following four tests are all going to be for the same inversion, testing if implementations are correct for
     * identifying the breakpoints by looking at different representations of evidence.
     * The inversion we are looking at is
     *
     * '+' strand representation: G....100....G|ACACA|C....100....C               A....100....A|TGTGT|T....100....T
     *
     * 100-bases of 'G' is the left flanking before the homology |ACACA| and the region starting with 100-bases of 'C' and
     * ending with 100-bases of 'A' and maybe (homologyForwardStrandRep uncertainty) the homology |TGTGT| is inverted.
     * 100-bases of 'T' is the right flanking region.
     *
     * Returns a list of four Tuple4's with left flanking evidence '+'/'-' strand representation and right flanking side.
     */
    private static List<TestDataForSimpleSVs>
    forSimpleInversionWithHomology(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();

        final byte[] leftLeftPlus = SVTestUtils.makeDummySequence(100, (byte)'G');
        final byte[] leftLeftMinus = SVTestUtils.makeDummySequence(100, (byte)'C');
        final byte[] leftRightPlus = SVTestUtils.makeDummySequence(100, (byte)'C');
        final byte[] leftRightMinus = SVTestUtils.makeDummySequence(100, (byte)'G');
        final byte[] rightLeftPlus = SVTestUtils.makeDummySequence(100, (byte)'A');
        final byte[] rightLeftMinus = SVTestUtils.makeDummySequence(100, (byte)'T');
        final byte[] rightRightPlus = SVTestUtils.makeDummySequence(100, (byte)'T');
        final byte[] rightRightMinus = SVTestUtils.makeDummySequence(100, (byte)'A');
        final byte[] leftHomology = "ACACA".getBytes();
        final byte[] rightHomology = "TGTGT".getBytes();
        {// left flanking evidence '+'/'-' strand representation
            outputStream.reset();
            outputStream.write(leftLeftPlus);outputStream.write(leftHomology);outputStream.write(rightLeftMinus);
            byte[] contigSeq = outputStream.toByteArray();

            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 101, 205), 1, 105, TextCigarCodec.decode("105M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 501, 605), 101, 205, TextCigarCodec.decode("100S105M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, new ArrayList<>(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
            result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

            outputStream.reset();
            outputStream.write(rightLeftPlus);outputStream.write(rightHomology);outputStream.write(leftLeftMinus);
            contigSeq = outputStream.toByteArray();
            region1 = new AlignmentInterval(new SimpleInterval("20", 501, 605), 1, 105, TextCigarCodec.decode("105M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("20", 101, 205), 101, 205, TextCigarCodec.decode("100S105M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, new ArrayList<>(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
            result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));
        }
        {// right flanking evidence '+'/'-' strand representation
            outputStream.reset();
            outputStream.write(leftRightMinus);outputStream.write(rightHomology);outputStream.write(rightRightPlus);
            byte[] contigSeq = outputStream.toByteArray();

            AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 201, 305), 1, 105, TextCigarCodec.decode("105M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 601, 705), 101, 205, TextCigarCodec.decode("100S105M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, new ArrayList<>(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
            result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

            outputStream.reset();
            outputStream.write(rightRightMinus);outputStream.write(leftHomology);outputStream.write(leftRightPlus);
            contigSeq = outputStream.toByteArray();

            region1 = new AlignmentInterval(new SimpleInterval("20", 601, 705), 1, 105, TextCigarCodec.decode("105M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            region2 = new AlignmentInterval(new SimpleInterval("20", 201, 305), 101, 205, TextCigarCodec.decode("100S105M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
            final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, new ArrayList<>(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
            result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));
        }
        return result;
    }

    /**
     * 40-'A' + 10-'C'+10-'T' + 40-'G' where the segment 10-'C'+10-'T' is deleted (forward strand representation description).
     *
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSVs>
    forSimpleDeletion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();
        // simple deletion '+' strand representation
        final byte[] leftRefFlank = SVTestUtils.makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = SVTestUtils.makeDummySequence(40, (byte)'G');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100040), 1 ,40, TextCigarCodec.decode("40M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100061, 100100), 41 ,80, TextCigarCodec.decode("40S40M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // simple deletion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentInterval(new SimpleInterval("21", 100061, 100100), 1 ,40, TextCigarCodec.decode("40M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100040), 41 ,80, TextCigarCodec.decode("40S40M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));

        return result;
    }

    /**
     * 100-'A' + 100-'T' and a 50 bases of 'C' is inserted at the A->T junction point (forward strand description)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSVs>
    forSimpleInsertion(final ByteArrayOutputStream outputStream) throws IOException {
        final List<TestDataForSimpleSVs> result = new ArrayList<>();

        // simple insertion '+' strand representation
        final byte[] leftRefFlank = SVTestUtils.makeDummySequence(100, (byte)'A');
        final byte[] insertedSeq  = SVTestUtils.makeDummySequence(50, (byte)'C');
        final byte[] rightRefFlank = SVTestUtils.makeDummySequence(100, (byte)'T');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(insertedSeq);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100100), 1 ,100, TextCigarCodec.decode("100M100S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100101, 100200), 151 ,250, TextCigarCodec.decode("100S100M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // simple insertion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(insertedSeq);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(insertedSeq);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentInterval(new SimpleInterval("21", 100101, 100200), 1 ,100, TextCigarCodec.decode("100M100S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100100), 151 ,250, TextCigarCodec.decode("100S100M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));

        return result;
    }

    /**
     * 50-'A' + 50-'C' where the middle 10-'A'+10-'C' is substituted with 10-'G' (forward strand representation)
     */
    private static List<TestDataForSimpleSVs>
    forLongRangeSubstitution() {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();

        // long range substitution '+' strand representation
        final byte[] leftRefFlank = SVTestUtils.makeDummySequence(50, (byte)'A');
        final byte[] rightRefFlank = SVTestUtils.makeDummySequence(50, (byte)'G');
        final byte[] substitution = SVTestUtils.makeDummySequence(10, (byte)'C');
        byte[] contigSeq = new byte[leftRefFlank.length+rightRefFlank.length-10];
        System.arraycopy(leftRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(substitution, 0, contigSeq, 40, substitution.length);
        System.arraycopy(rightRefFlank, 0, contigSeq, 50, 40);
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100040), 1 ,40, TextCigarCodec.decode("40M50S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100061, 100100), 51 ,90, TextCigarCodec.decode("50S40M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // long range substitution '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(substitution);
        System.arraycopy(rightRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(substitution, 0, contigSeq, 40, substitution.length);
        System.arraycopy(leftRefFlank, 0, contigSeq, 40 + substitution.length, 40);
        region1 = new AlignmentInterval(new SimpleInterval("21", 100061, 100100), 1 ,40, TextCigarCodec.decode("40M50S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100040), 51 ,90, TextCigarCodec.decode("50S40M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));

        return result;
    }

    /**
     * 40-'C' + 'ATCG' + 34 bases of unique sequence + 'ATCG' + 40-'T' is shrunk to 40-'C' + 'ATCG' + 40-'T' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSVs>
    forDeletionWithHomology(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();

        // simple deletion with homology '+' strand representation
        final byte[] leftRefFlank = SVTestUtils.makeDummySequence(40, (byte)'C');
        final byte[] rightRefFlank = SVTestUtils.makeDummySequence(40, (byte)'T');
        final byte[] homology = new byte[]{'A', 'T', 'C', 'G'};
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(homology);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100044), 1 ,44, TextCigarCodec.decode("44M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100079, 100122), 41 ,84, TextCigarCodec.decode("40S44M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // simple deletion with homology '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(homology);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(homology);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentInterval(new SimpleInterval("21", 100079, 100122), 1 ,44, TextCigarCodec.decode("44M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100044), 41 ,84, TextCigarCodec.decode("40S44M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));

        return result;
    }

    /**
     * 40-'A' + 20-'C' + 40-'G' is shrunk to 40-'A' + 10-'C' + 40-'G' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSVs>
    forSimpleTandemDuplicationContraction() {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();

        // simple tandem duplication contraction '+' strand representation
        final byte[] leftRefFlank = SVTestUtils.makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = SVTestUtils.makeDummySequence(40, (byte)'G');
        final byte[] doubleDup = SVTestUtils.makeDummySequence(20, (byte)'C');
        final byte[] contigSeq = new byte[90];
        System.arraycopy(leftRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(doubleDup, 0, contigSeq, 40, 10);
        System.arraycopy(rightRefFlank, 0, contigSeq, 50, 40);

        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100050), 1 ,50, TextCigarCodec.decode("50M40S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100051, 100100), 41 ,90, TextCigarCodec.decode("40S50M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // simple tandem duplication contraction '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(doubleDup);
        System.arraycopy(rightRefFlank, 0, contigSeq, 0, 40);
        System.arraycopy(doubleDup, 0, contigSeq, 40, 10);
        System.arraycopy(leftRefFlank, 0, contigSeq, 50, 40);
        region1 = new AlignmentInterval(new SimpleInterval("21", 100051, 100100), 1 ,50, TextCigarCodec.decode("50M40S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100050), 41 ,90, TextCigarCodec.decode("40S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));

        return result;
    }

    /**
     * 40-'A' + 10-'C' + 40-'G' is expanded to 40-'A' + 20-'C' + 40-'G' (forward strand representation)
     * Return a list of two entries for positive and reverse strand representations.
     */
    private static List<TestDataForSimpleSVs>
    forSimpleTandemDuplicationExpansion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();

        // simple tandem duplication expansion '+' strand representation
        final byte[] leftRefFlank = SVTestUtils.makeDummySequence(40, (byte)'A');
        final byte[] rightRefFlank = SVTestUtils.makeDummySequence(40, (byte)'G');
        final byte[] doubleDup = SVTestUtils.makeDummySequence(20, (byte)'C');
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(doubleDup);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();

        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 100001, 100050), 1 ,50, TextCigarCodec.decode("50M50S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 100041, 100090), 51 ,100, TextCigarCodec.decode("50S50M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // simple tandem duplication expansion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(doubleDup);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(doubleDup);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();
        region1 = new AlignmentInterval(new SimpleInterval("21", 100041, 100090), 1 ,50, TextCigarCodec.decode("50M50S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 100001, 100050), 51 ,100, TextCigarCodec.decode("50S50M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));

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
    private static List<TestDataForSimpleSVs>
    forSimpleTandemDuplicationExpansionWithNovelInsertion(final ByteArrayOutputStream outputStream) throws IOException {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();
        // simple tandem duplication expansion with novel insertion '+' strand representation
        final byte[] leftRefFlank = "GTTAGTAGATATTCTAGCTGACTCAGTTCAGTGTTGCTATGATTAAACAAGAGTGAGTTCCCT".getBytes();                     //63
        final byte[] rightRefFlank = "CATTATTGATATTTCATTATGTTCAACAGATGGAGTTAATGTGAATGT".getBytes();                                   //48
        final byte[] insertedSeq = "CTCTCTCTCT".getBytes();                                                                           //10
        final byte[] dup = "AAAAGTAAATGTTATAAGAAATCTTAAGTATTATTTTCTTATGTTTCTAGCCTAATAAAGTGCTTTTATTAAAGCACTTTATTTAAAGG".getBytes();    //89
        outputStream.reset();
        outputStream.write(leftRefFlank);outputStream.write(dup);outputStream.write(insertedSeq);outputStream.write(dup);outputStream.write(rightRefFlank);
        byte[] contigSeq = outputStream.toByteArray();

        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("21", 25297101, 25297252), 1 ,152, TextCigarCodec.decode("152M147S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("21", 25297164, 25297300), 163 ,299, TextCigarCodec.decode("162S137M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpoints = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // simple tandem duplication expansion with novel insertion '-' strand representation
        SequenceUtil.reverseComplement(leftRefFlank);
        SequenceUtil.reverseComplement(rightRefFlank);
        SequenceUtil.reverseComplement(insertedSeq);
        SequenceUtil.reverseComplement(dup);
        outputStream.reset();
        outputStream.write(rightRefFlank);outputStream.write(dup);outputStream.write(insertedSeq);outputStream.write(dup);outputStream.write(leftRefFlank);
        contigSeq = outputStream.toByteArray();

        region1 = new AlignmentInterval(new SimpleInterval("21", 25297164, 25297300), 1 ,137, TextCigarCodec.decode("137M162S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("21", 25297101, 25297252), 148 ,299, TextCigarCodec.decode("147S152M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final NovelAdjacencyAndAltHaplotype breakpointsDetectedFromReverseStrand = new NovelAdjacencyAndAltHaplotype(new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict), contigSeq, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpointsDetectedFromReverseStrand, "asm000001:tig00001"));

        return result;
    }

    /**
     * These test data was based on a real observation on a locally-assembled contig
     * "TGCCAGGTTACATGGCAAAGAGGGTAGATATGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGAGGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGAGGGCAGCTGTGGATGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC"
     * with two alignment records chr18:312579-312718 140M135S
     *                            chr18:312610-312757 127S148M
     * for a tandem repeat expansion event from 1 copy to 2 copies with also a pseudo-homology

     * Return a list of eight entries for positive and reverse strand representations for:
     * 1. expansion from 1 unit to 2 units with pseudo-homology
     * 2. contraction from 2 units to 1 unit with pseudo-homology
     * 3. contraction from 3 units to 2 units without pseudo-homology
     * 4. expansion from 2 units to 3 units without pseudo-homology
     */
    private static List<TestDataForSimpleSVs>
    forComplexTandemDuplication() {

        final List<TestDataForSimpleSVs> result = new ArrayList<>();
        final String leftRefFlank       = "TGCCAGGTTACATGGCAAAGAGGGTAGATAT";                                                                    // 31
        final String rightRefFlank      = "TGGTGCAAATGCCATTTATGCTCCTCTCCACCCATATCC";                                                            // 39
        final String firstRepeat        = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCATGA";   // 96
        final String secondRepeat       = "GGGGAGCTGTGAAGAATGGAGCCAGTAATTAAATTCACTGAAGTCTCCACAGGAGGGCAAGGTGGACAATCTGTCCCATAGGAGGGGGATTCAGGA";   // 96
        final String pseudoHomology     = "GGGCAGCTGTGGA";                                                                                      // 13


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        final byte[] fakeRefSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, pseudoHomology, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexExpansionWithPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, pseudoHomology, rightRefFlank).getBytes();
        AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 1 ,140, TextCigarCodec.decode("140M135S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312757), 128 ,275, TextCigarCodec.decode("127S148M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        NovelAdjacencyAndAltHaplotype breakpoints =
                new NovelAdjacencyAndAltHaplotype(
                        new ChimericAlignment(region1, region2, Collections.emptyList(), "asm000001:tig00001", b37_seqDict),
                        contigSeqForComplexExpansionWithPseudoHomology, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        final byte[] contigSeqForComplexExpansionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionWithPseudoHomology, contigSeqForComplexExpansionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionWithPseudoHomology_reverseStrand);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312757), 1 ,148, TextCigarCodec.decode("148M127S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 136 ,275, TextCigarCodec.decode("135S140M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        breakpoints = new NovelAdjacencyAndAltHaplotype(
                new ChimericAlignment(region1, region2, Collections.emptyList(),
                        "asm000001:tig00001", b37_seqDict),
                contigSeqForComplexExpansionWithPseudoHomology_reverseStrand, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // second test: contraction from 2 units to 1 unit with pseudo-homology
        final byte[] contigSeqForComplexContractionWithPseudoHomology = fakeRefSeqForComplexExpansionWithPseudoHomology;
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 1, 140, TextCigarCodec.decode("140M39S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312706, 312853), 32, 179, TextCigarCodec.decode("31S148M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        breakpoints =
                new NovelAdjacencyAndAltHaplotype(
                        new ChimericAlignment(region1, region2, Collections.emptyList(),
                                "asm000001:tig00001", b37_seqDict),
                        contigSeqForComplexContractionWithPseudoHomology, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        final byte[] contigSeqForComplexContractionWithPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionWithPseudoHomology, contigSeqForComplexContractionWithPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionWithPseudoHomology_reverseStrand);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312706, 312853), 1, 148, TextCigarCodec.decode("148M31S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312718), 40, 179, TextCigarCodec.decode("39S140M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        breakpoints = new NovelAdjacencyAndAltHaplotype(
                new ChimericAlignment(region1, region2, Collections.emptyList(),
                        "asm000001:tig00001", b37_seqDict),
                contigSeqForComplexContractionWithPseudoHomology_reverseStrand, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // third test: contraction from 3 units to 2 units without pseudo-homology
        final byte[] fakeRefSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, firstRepeat, rightRefFlank).getBytes();
        final byte[] contigSeqForComplexContractionNoPseudoHomology = String.format("%s%s%s%s", leftRefFlank, firstRepeat, secondRepeat, rightRefFlank).getBytes();
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 1, 223, TextCigarCodec.decode("223M39S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312706, 312936), 32, 262, TextCigarCodec.decode("31S231M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        breakpoints =
                new NovelAdjacencyAndAltHaplotype(
                        new ChimericAlignment(region1, region2, Collections.emptyList(),
                                "asm000001:tig00001", b37_seqDict),
                        contigSeqForComplexContractionNoPseudoHomology, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        final byte[] contigSeqForComplexContractionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexContractionNoPseudoHomology, contigSeqForComplexContractionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexContractionNoPseudoHomology_reverseStrand);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312706, 312936), 1, 231, TextCigarCodec.decode("231M31S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 40, 262, TextCigarCodec.decode("39S223M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        breakpoints = new NovelAdjacencyAndAltHaplotype(
                new ChimericAlignment(region1, region2, Collections.emptyList(),
                        "asm000001:tig00001", b37_seqDict),
                contigSeqForComplexContractionNoPseudoHomology_reverseStrand, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        final byte[] contigSeqForComplexExpansionNoPseudoHomology = fakeRefSeqForComplexContractionNoPseudoHomology;
        region1 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 1, 223, TextCigarCodec.decode("223M135S"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312610, 312840), 128, 358, TextCigarCodec.decode("127S231M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        breakpoints =
                new NovelAdjacencyAndAltHaplotype(
                        new ChimericAlignment(region1, region2, Collections.emptyList(),
                                "asm000001:tig00001", b37_seqDict),
                        contigSeqForComplexExpansionNoPseudoHomology, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        final byte[] contigSeqForComplexExpansionNoPseudoHomology_reverseStrand = Arrays.copyOf(contigSeqForComplexExpansionNoPseudoHomology, contigSeqForComplexExpansionNoPseudoHomology.length);
        SequenceUtil.reverseComplement(contigSeqForComplexExpansionNoPseudoHomology_reverseStrand);
        region1 = new AlignmentInterval(new SimpleInterval("20", 312610, 312840), 1, 231, TextCigarCodec.decode("231M127S"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        region2 = new AlignmentInterval(new SimpleInterval("20", 312579, 312801), 136, 358, TextCigarCodec.decode("135S223M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        breakpoints = new NovelAdjacencyAndAltHaplotype(
                new ChimericAlignment(region1, region2, Collections.emptyList(),
                        "asm000001:tig00001", b37_seqDict),
                contigSeqForComplexExpansionNoPseudoHomology_reverseStrand, b37_seqDict);
        result.add(new TestDataForSimpleSVs(region1, region2, breakpoints, "asm000001:tig00001"));

        return result;
    }
}
