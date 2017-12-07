package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SVDiscoveryTestDataProvider;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple4;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SVDiscoveryTestDataProvider.*;

public class ChimericAlignmentUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testFilterByRegionTooSmall() {
        final byte[] contigSequence = SVDiscoveryTestDataProvider.LONG_CONTIG1.getBytes();
        final AlignmentInterval region1 = new AlignmentInterval(new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20138007, 20142231), 1, contigSequence.length - 1986, TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval region2 = new AlignmentInterval(new SimpleInterval(SVDiscoveryTestDataProvider.chrForLongContig1, 20152030, 20154634), 3604, contigSequence.length, TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, 60, 36, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( ChimericAlignment.firstAlignmentIsTooShort(region1, region2, 3000) );
        Assert.assertTrue( ChimericAlignment.firstAlignmentIsTooShort(region2, region1, 3000) );
    }

    @Test(groups = "sv")
    public void testFilterByNextAlignmentMayBeInsertion() {
        final AlignmentInterval overlappingRegion1 = new AlignmentInterval(new SimpleInterval("19", 48699881, 48700035), 1, 154, TextCigarCodec.decode("47S154M"), false, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval overlappingRegion2 = new AlignmentInterval(new SimpleInterval("19", 48700584, 48700669), 117, 201, TextCigarCodec.decode("116H85M"), true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        Assert.assertTrue(ChimericAlignment.nextAlignmentMayBeInsertion(overlappingRegion1, overlappingRegion2,  CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD, 50,true));
    }

    @DataProvider(name = "forBooleanSeriesAndSerialization")
    private Object[][] createTestData() {
        final List<Object[]> data = new ArrayList<>(20);

        // simple inversion
        Tuple4<AlignmentInterval, AlignmentInterval, NovelAdjacencyReferenceLocations, String> testData = forSimpleInversionFromLongCtg1WithStrangeLeftBreakpoint;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.REVERSE_TO_FORWARD, true, false, false, false, false});

        testData = forSimpleInversionWithHom_leftPlus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.FORWARD_TO_REVERSE, true, false, false, false, false});

        testData = forSimpleInversionWithHom_leftMinus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.FORWARD_TO_REVERSE, false, false, false, false, false});

        testData = forSimpleInversionWithHom_rightPlus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.REVERSE_TO_FORWARD, true, false, false, false, false});

        testData = forSimpleInversionWithHom_rightMinus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.REVERSE_TO_FORWARD, false, false, false, false, false});

        // simple deletion
        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleDeletion_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});

        // simple insertion
        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleInsertion_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});

        // long range substitution
        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forLongRangeSubstitution_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});

        // simple deletion with homology
        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forDeletionWithHomology_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});

        // tandem duplication simple contraction
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleTanDupContraction_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});

        // tandem duplication simple expansion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansion_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});

        // tandem duplication simple expansion with novel insertion
        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forSimpleTanDupExpansionWithNovelIns_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});



        // first test (the original observed event, but assigned to a different chromosome): expansion from 1 unit to 2 units with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_1to2_pseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});


        // second test: contraction from 2 units to 1 unit with pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to1_pseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});


        // third test: contraction from 3 units to 2 units without pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_3to2_noPseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});


        // fourth test: expansion from 2 units to 3 units without pseudo-homology
        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_plus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, true, false, false, false, false});

        testData = SVDiscoveryTestDataProvider.forComplexTanDup_2to3_noPseudoHom_minus;
        data.add(new Object[]{testData._1(), testData._2(), StrandSwitch.NO_SWITCH, false, false, false, false, false});

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

        data.add(new Object[]{intervalOne, intervalTwo, StrandSwitch.REVERSE_TO_FORWARD, false, false, true, true, true});

        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 48513458, 48513545), 1, 88, TextCigarCodec.decode("88M227H"), true, 39, 1, 83, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 48513297, 48513579), 84, 365, TextCigarCodec.decode("83S282M"), false, 60, 0, 282, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{intervalOne, intervalTwo, StrandSwitch.FORWARD_TO_REVERSE, true, false, true, true, true});


        // same-chr translocation suspect, forward and reverse representation
        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 61015129, 61015272), 1, 144, TextCigarCodec.decode("144M148H"), true, 60, 1, 139, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 60992732, 60992880), 144, 292, TextCigarCodec.decode("143S149M"), true, 60, 0, 149, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{intervalOne, intervalTwo, StrandSwitch.NO_SWITCH, true, true, false, false, true});

        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 28861368, 28861775), 1, 409, TextCigarCodec.decode("387M1I21M623H"), false, 60, 22, 286, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28896473, 28897229), 276, 1032, TextCigarCodec.decode("275S757M"), false, 60, 1, 752, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{intervalOne, intervalTwo, StrandSwitch.NO_SWITCH, false, true, false, false, true});

        // diff-chr translocation suspect without SS
        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 24923683, 24923715), 1, 33, TextCigarCodec.decode("33M130H"), true, 60, 0, 33, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 11590055, 11590197), 21, 163, TextCigarCodec.decode("20S143M"), true, 60, 3, 128, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{intervalOne, intervalTwo, StrandSwitch.NO_SWITCH, true, true, false, false, true});

        // diff-chr translocation suspect with SS
        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 5374092, 5374748), 1, 656, TextCigarCodec.decode("656M322S"), true, 60, 14, 586, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr20", 28764673, 28765145), 506, 978, TextCigarCodec.decode("473M505H"), false, 60, 16, 393, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{intervalOne, intervalTwo, StrandSwitch.FORWARD_TO_REVERSE, false, true, false, false, true});

        // same-chr reference order switch, but overlaps (hence incomplete picture)
        intervalOne = new AlignmentInterval(new SimpleInterval("20", 283, 651), 383, 751, TextCigarCodec.decode("382H369M274H"), true, 60, 23, 254, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("20", 1, 413), 613, 1025, TextCigarCodec.decode("612H413M"), true, 60, 0, 413, ContigAlignmentsModifier.AlnModType.NONE);
        data.add(new Object[]{intervalOne, intervalTwo, StrandSwitch.NO_SWITCH, true, false, false, true, true});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forBooleanSeriesAndSerialization", groups = "sv")
    public void test(final AlignmentInterval region1, final AlignmentInterval region2,
                     final StrandSwitch expectedStrandSwitch,
                     final boolean expectedIsForwardStrandRepresentation,
                     final boolean expectedIsLikelySimpleTranslocation,
                     final boolean expectedIsLikelyInvDup,
                     final boolean expectedIsIncompletePicture,
                     final boolean shouldUseRefVersion_b38) {

        final SAMSequenceDictionary refDict = shouldUseRefVersion_b38 ? SVDiscoveryTestDataProvider.b38_seqDict : SVDiscoveryTestDataProvider.seqDict;

        Assert.assertEquals(ChimericAlignment.determineStrandSwitch(region1, region2), expectedStrandSwitch);
        Assert.assertEquals(ChimericAlignment.isForwardStrandRepresentation(region1, region2, expectedStrandSwitch, refDict), expectedIsForwardStrandRepresentation);
        Assert.assertEquals(ChimericAlignment.isLikelySimpleTranslocation(region1, region2, expectedStrandSwitch), expectedIsLikelySimpleTranslocation);
        Assert.assertEquals(ChimericAlignment.isLikelyInvertedDuplication(region1, region2), expectedIsLikelyInvDup);
        Assert.assertEquals(ChimericAlignment.hasIncompletePictureFromTwoAlignments(region1, region2), expectedIsIncompletePicture);

        final ChimericAlignment chimericAlignment = new ChimericAlignment(region1, region2, Collections.emptyList(), "dummyName", refDict);
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, chimericAlignment);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final ChimericAlignment roundTrip = (ChimericAlignment) kryo.readClassAndObject(in);
        Assert.assertEquals(roundTrip, chimericAlignment);
        Assert.assertEquals(roundTrip.hashCode(), chimericAlignment.hashCode());
    }
}