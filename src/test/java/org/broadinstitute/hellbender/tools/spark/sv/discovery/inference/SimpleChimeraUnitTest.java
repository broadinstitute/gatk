package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple3;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.*;

public class SimpleChimeraUnitTest extends GATKBaseTest {

    static List<Tuple3<AlignmentInterval, AlignmentInterval, SAMSequenceDictionary>> alignmentPairsForSimpleChimeraAndRefSeqDict() {

        final List<Tuple3<AlignmentInterval, AlignmentInterval, SAMSequenceDictionary>> result = new ArrayList<>(20);

        // simple inversion
        TestDataForSimpleSVs testData = forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = forSimpleInversionWithHom_leftPlus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = forSimpleInversionWithHom_leftMinus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = forSimpleInversionWithHom_rightPlus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = forSimpleInversionWithHom_rightMinus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        // simple deletion
        testData = SimpleSVDiscoveryTestDataProvider.forSimpleDeletion_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forSimpleDeletion_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        // simple insertion
        testData = SimpleSVDiscoveryTestDataProvider.forSimpleInsertion_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forSimpleInsertion_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        // long range substitution
        testData = SimpleSVDiscoveryTestDataProvider.forLongRangeSubstitution_fudgedDel_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forLongRangeSubstitution_fudgedDel_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        // simple deletion with homology
        testData = SimpleSVDiscoveryTestDataProvider.forDeletionWithHomology_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forDeletionWithHomology_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        // tandem duplication simple contraction
        testData = SimpleSVDiscoveryTestDataProvider.forSimpleTanDupContraction_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        // tandem duplication simple expansion
        testData = SimpleSVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forSimpleTanDupExpansion_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        // tandem duplication simple expansion with novel insertion
        testData = SimpleSVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));


        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));


        // second test: contraction from 2 units to 1 unit with pseudo-homology
        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));


        // third test: contraction from 3 units to 2 units without pseudo-homology
        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));


        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));

        testData = SimpleSVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;
        result.add(new Tuple3<>(testData.firstAlignment, testData.secondAlignment, SimpleSVDiscoveryTestDataProvider.b37_seqDict));


        ////////// ABOVE ARE FOR SIMPLE VARIANTS: INS/DEL, DUP EXPANSION, DUP CONTRACTION, INVERSION, BELOW ARE FOR TRANSLOCATION SUSPECTS AND INV DUP

        // inverted duplication
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(21, 1, 46709983);
        SAMRecord one =
                ArtificialReadUtils.createArtificialRead(header, "asm017968:tig00020", 20, 25625477,
                        "CCTTTCTATTCTAACATTATTGACCCACTACAATAAAATTAAGTATACTTTAGGCTGGGCATGGTGGTTTACACCTGTAATCCCAACACTTTGGGAGGCCGAGGCGGGTGG".getBytes(),
                        ArtificialReadUtils.createRandomReadQuals(111), "212H111M").convertToSAMRecord(header);
        one.setSupplementaryAlignmentFlag(true);
        one.setMappingQuality(60);
        one.setReadNegativeStrandFlag(true);
        one.setAttribute("NM", 0);
        one.setAttribute("AS", 111);

        SAMRecord two = ArtificialReadUtils.createArtificialRead(header, "asm017968:tig00020", 20, 25625379,
                "CCACCCGCCTCGGCCTCCCAAAGTGTTGGGATTACAGGTGTAAACCACCATGCCCAGCCTAAAGTATACTTAATTTTATTGTAGTGGGTCAATAATGTTAGAATAGAAAGGAGAATCAGAAGAGAGAAAAGAAACAGATGATCTTATGAATCCCCAATTTATATCCCCCAATTAAGGCATCTCTTTCTTCTGTCTCCCTATATCCCTTTCTATTCTAACATTATTGACCCACTACAATAAAATTAAGTATACTTTAGGCTGGGCATGGTGGTTTACACCTGTAATCCCAACACTTTGGGAGGCCGAGGCGGGTGGATCACTTG".getBytes(),
                ArtificialReadUtils.createRandomReadQuals(323), "106S217M").convertToSAMRecord(header);
        two.setMappingQuality(60);
        two.setAttribute("NM", 0);
        two.setAttribute("AS", 217);

        AlignmentInterval intervalOne = new AlignmentInterval(one);
        AlignmentInterval intervalTwo = new AlignmentInterval(two);
        result.add(new Tuple3<>(intervalOne, intervalTwo, SimpleSVDiscoveryTestDataProvider.b38_seqDict));

        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 48513458, 48513545), 1, 88, TextCigarCodec.decode("88M227H"), true, 39, 1, 83, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 48513297, 48513578), 84, 365, TextCigarCodec.decode("83S282M"), false, 60, 0, 282, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new Tuple3<>(intervalOne, intervalTwo, SimpleSVDiscoveryTestDataProvider.b38_seqDict));


        // same-chr translocation suspect, forward and reverse representation
        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 61015129, 61015272), 1, 144, TextCigarCodec.decode("144M148H"), true, 60, 1, 139, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 60992732, 60992880), 144, 292, TextCigarCodec.decode("143S149M"), true, 60, 0, 149, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new Tuple3<>(intervalOne, intervalTwo, SimpleSVDiscoveryTestDataProvider.b38_seqDict));

        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 28861368, 28861775), 1, 409, TextCigarCodec.decode("387M1I21M623H"), false, 60, 22, 286, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28896473, 28897229), 276, 1032, TextCigarCodec.decode("275S757M"), false, 60, 1, 752, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new Tuple3<>(intervalOne, intervalTwo, SimpleSVDiscoveryTestDataProvider.b38_seqDict));

        // diff-chr translocation suspect without SS
        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 24923683, 24923715), 1, 33, TextCigarCodec.decode("33M130H"), true, 60, 0, 33, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 11590055, 11590197), 21, 163, TextCigarCodec.decode("20S143M"), true, 60, 3, 128, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new Tuple3<>(intervalOne, intervalTwo, SimpleSVDiscoveryTestDataProvider.b38_seqDict));

        // diff-chr translocation suspect with SS
        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 5374092, 5374747), 1, 656, TextCigarCodec.decode("656M322S"), true, 60, 14, 586, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28764673, 28765145), 506, 978, TextCigarCodec.decode("473M505H"), false, 60, 16, 393, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new Tuple3<>(intervalOne, intervalTwo, SimpleSVDiscoveryTestDataProvider.b38_seqDict));

        // same-chr reference order switch, but overlaps (hence incomplete picture)
        intervalOne = new AlignmentInterval(new SimpleInterval("20", 283, 651), 383, 751, TextCigarCodec.decode("382H369M274H"), true, 60, 23, 254, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("20", 1, 413), 613, 1025, TextCigarCodec.decode("612H413M"), true, 60, 0, 413, ContigAlignmentsModifier.AlnModType.NONE);
        result.add(new Tuple3<>(intervalOne, intervalTwo, SimpleSVDiscoveryTestDataProvider.b38_seqDict));

        return result;
    }

    @DataProvider(name = "testRepresentationAndSerialization")
    private Object[][] testRepresentationAndSerialization() {
        final List<Tuple3<AlignmentInterval, AlignmentInterval, SAMSequenceDictionary>> tuple3s = alignmentPairsForSimpleChimeraAndRefSeqDict();
        final List<Object[]> data = new ArrayList<>(tuple3s.size());

        int i=0;
        // simple inversion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, false}); ++i;

        // simple deletion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;

        // simple insertion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;

        // long range substitution
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;

        // simple deletion with homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;

        // tandem duplication simple contraction
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;

        // tandem duplication simple expansion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;

        // tandem duplication simple expansion with novel insertion
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;



        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;


        // second test: contraction from 2 units to 1 unit with pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;


        // third test: contraction from 3 units to 2 units without pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;


        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;


        data.add(new Object[]{tuple3s.get(i), StrandSwitch.REVERSE_TO_FORWARD, false}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, true}); ++i;


        // same-chr translocation suspect, forward and reverse representation
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, false}); ++i;

        // diff-chr translocation suspect without SS
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        // diff-chr translocation suspect with SS
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.FORWARD_TO_REVERSE, false}); ++i;

        // same-chr reference order switch, but overlaps (hence incomplete picture)
        data.add(new Object[]{tuple3s.get(i), StrandSwitch.NO_SWITCH, true}); ++i;

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "testRepresentationAndSerialization", groups = "sv")
    public void testRepresentationAndSerialization(Tuple3<AlignmentInterval, AlignmentInterval, SAMSequenceDictionary> chimericPairsAndRefSeqDict,
                                                   final StrandSwitch expectedStrandSwitch,
                                                   final boolean expectedIsForwardStrandRepresentation) {
        final AlignmentInterval region1 = chimericPairsAndRefSeqDict._1();
        final AlignmentInterval region2 = chimericPairsAndRefSeqDict._2();
        final SAMSequenceDictionary refDict = chimericPairsAndRefSeqDict._3();

        Assert.assertEquals(SimpleChimera.determineStrandSwitch(region1, region2), expectedStrandSwitch);
        Assert.assertEquals(SimpleChimera.isForwardStrandRepresentation(region1, region2, expectedStrandSwitch, refDict), expectedIsForwardStrandRepresentation);

        final SimpleChimera simpleChimera = new SimpleChimera(region1, region2, Collections.emptyList(), "dummyName", refDict);
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, simpleChimera);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final SimpleChimera roundTrip = (SimpleChimera) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, simpleChimera);
        Assert.assertEquals(roundTrip.hashCode(), simpleChimera.hashCode());
    }
}