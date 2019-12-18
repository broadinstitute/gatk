package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.*;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.CodonTracker.NO_FRAME_SHIFT_CODON;
import static org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis.*;

public class AnalyzeSaturationMutagenesisUnitTest extends GATKBaseTest {
    final private static byte CALL_A = (byte)'A';
    final private static byte CALL_C = (byte)'C';
    final private static byte CALL_G = (byte)'G';
    final private static byte CALL_T = (byte)'T';
    final private static byte NO_CALL = (byte)'-';
    final private static byte QUAL_10 = (byte)10;
    final private static byte QUAL_20 = (byte)20;
    final private static byte QUAL_30 = (byte)30;
    final private static byte[] refSeq = "ACATGCGTCTAGTACGT".getBytes();
    final private static String orfCoords = "3-6,8-12";
    final private static CodonTracker localCodonTracker = new CodonTracker(orfCoords, refSeq, logger);

    @Test
    public void testInterval() {
        final Interval interval1 = new Interval(1, 3);
        Assert.assertEquals(interval1.size(), 2);
        final Interval interval2 = new Interval(1, 3);
        Assert.assertEquals(interval1, interval2);
        Assert.assertEquals(interval1.hashCode(), interval2.hashCode());
        final Interval interval3 = new Interval(1, 4);
        Assert.assertNotEquals(interval1, interval3);
        Assert.assertNotEquals(interval1.hashCode(), interval3.hashCode());
    }

    @Test
    public void testSNV() {
        final SNV snv1 = new SNV(0, CALL_A, CALL_C, QUAL_30);
        final SNV snv2 = new SNV(0, CALL_A, CALL_C, QUAL_20);
        Assert.assertEquals(snv1.hashCode(), snv2.hashCode());
        Assert.assertEquals(snv1, snv2);
        Assert.assertEquals(snv1.compareTo(snv2), 0);
        final SNV snv3 = new SNV(1, CALL_A, CALL_C, QUAL_30);
        Assert.assertNotEquals(snv1.hashCode(), snv3.hashCode());
        Assert.assertNotEquals(snv1, snv3);
        Assert.assertTrue(snv1.compareTo(snv3) < 0);
        Assert.assertTrue(snv3.compareTo(snv2) > 0);
        final SNV snv4 = new SNV(0, CALL_G, CALL_C, QUAL_30);
        Assert.assertNotEquals(snv1.hashCode(), snv4.hashCode());
        Assert.assertNotEquals(snv1, snv4);
        Assert.assertTrue(snv1.compareTo(snv4) < 0);
        Assert.assertTrue(snv4.compareTo(snv1) > 0);
        final SNV snv5 = new SNV(0, CALL_A, CALL_G, QUAL_30);
        Assert.assertNotEquals(snv1.hashCode(), snv5.hashCode());
        Assert.assertNotEquals(snv1, snv5);
        Assert.assertTrue(snv1.compareTo(snv5) < 0);
        Assert.assertTrue(snv5.compareTo(snv1) > 0);
        Assert.assertEquals(snv5.toString(), "1:A>G");
    }

    @Test
    public void testIntervalCounter() {
        final IntervalCounter intervalCounter = new IntervalCounter(10);
        intervalCounter.addCount(new Interval(0, 10));
        intervalCounter.addCount(new Interval(1, 9));
        intervalCounter.addCount(new Interval(2, 8));
        intervalCounter.addCount(new Interval(3, 7));
        intervalCounter.addCount(new Interval(4, 6));
        Assert.assertEquals(intervalCounter.countSpanners(0, 10), 1);
        Assert.assertEquals(intervalCounter.countSpanners(5, 5), 5);
        Assert.assertEquals(intervalCounter.countSpanners(2, 5), 3);
        Assert.assertEquals(intervalCounter.countSpanners(5, 8), 3);
        intervalCounter.addCount(new Interval(0, 10));
        Assert.assertEquals(intervalCounter.countSpanners(0, 10), 2);
    }

    @Test
    public void testSNVCollection() {
        final List<SNV> snvList1 = new ArrayList<>(Arrays.asList(
                new SNV(0, CALL_A, CALL_C, QUAL_30),
                new SNV(1, CALL_A, CALL_C, QUAL_20)));
        final SNVCollectionCount cc1 = new SNVCollectionCount(snvList1, 10);

        // equality, key, compare, and hash should be independent of count and coverage
        final List<SNV> snvList2 = Arrays.asList(
                new SNV(0, CALL_A, CALL_C, QUAL_30),
                new SNV(1, CALL_A, CALL_C, QUAL_20));
        final SNVCollectionCount cc2 = new SNVCollectionCount(snvList2, 20);
        Assert.assertEquals(cc1.hashCode(), cc2.hashCode());
        Assert.assertEquals(cc1, cc2);
        Assert.assertEquals(cc1.compareTo(cc2), 0);
        cc2.bumpCount(30);
        Assert.assertEquals(cc1.hashCode(), cc2.hashCode());
        Assert.assertEquals(cc1, cc2);
        Assert.assertEquals(cc1.compareTo(cc2), 0);

        Assert.assertEquals(cc2.getCount(), 2);
        Assert.assertEquals(cc2.getMeanRefCoverage(), 25., .0000001);

        // changing the list shouldn't change the hash or the key
        final int cc1Hash = cc1.hashCode();
        final List<SNV> key1 = cc1.getSNVs();
        snvList1.add(new SNV(2, CALL_A, CALL_C, QUAL_10));
        Assert.assertEquals(cc1.hashCode(), cc1Hash);
        Assert.assertEquals(cc1.getSNVs(), key1);

        // different lists should mean unequal to each other, unequal hashes, and non-zero compare
        final SNVCollectionCount cc3 = new SNVCollectionCount(snvList1, 20);
        Assert.assertNotEquals(cc1.hashCode(), cc3.hashCode());
        Assert.assertNotEquals(cc1, cc3);
        Assert.assertTrue(cc1.compareTo(cc3) < 0);
        Assert.assertTrue(cc3.compareTo(cc1) > 0);
    }

    @Test
    public void testCodonVariationBasics() {
        Assert.assertTrue(CodonVariation.createDeletion(0).isDeletion());
        Assert.assertTrue(CodonVariation.createFrameshift(0).isFrameshift());
        Assert.assertTrue(CodonVariation.createInsertion(0,0).isInsertion());
        Assert.assertTrue(CodonVariation.createModification(0,0).isModification());
        Assert.assertFalse(CodonVariation.createDeletion(0).isModification());

        Assert.assertEquals(CodonVariation.createDeletion(1), CodonVariation.createDeletion(1));
        Assert.assertNotEquals(CodonVariation.createDeletion(0), CodonVariation.createDeletion(1));
        Assert.assertEquals(CodonVariation.createFrameshift(1), CodonVariation.createFrameshift(1));
        Assert.assertNotEquals(CodonVariation.createFrameshift(0), CodonVariation.createFrameshift(1));
        Assert.assertNotEquals(CodonVariation.createFrameshift(0), CodonVariation.createDeletion(0));
        Assert.assertEquals(CodonVariation.createInsertion(1, 0),
                            CodonVariation.createInsertion(1, 0));
        Assert.assertNotEquals(CodonVariation.createInsertion(0, 0),
                               CodonVariation.createInsertion(1, 0));
        Assert.assertNotEquals(CodonVariation.createInsertion(0, 0),
                               CodonVariation.createInsertion(0, 1));
        Assert.assertEquals(CodonVariation.createModification(1, 0),
                            CodonVariation.createModification(1, 0));
        Assert.assertNotEquals(CodonVariation.createModification(0, 0),
                               CodonVariation.createModification(1, 0));
        Assert.assertNotEquals(CodonVariation.createModification(0, 0),
                               CodonVariation.createModification(0, 1));

        Assert.assertEquals(CodonVariation.createDeletion(1).hashCode(), CodonVariation.createDeletion(1).hashCode());
        Assert.assertNotEquals(CodonVariation.createDeletion(0).hashCode(), CodonVariation.createDeletion(1).hashCode());
        Assert.assertEquals(CodonVariation.createFrameshift(1).hashCode(), CodonVariation.createFrameshift(1).hashCode());
        Assert.assertNotEquals(CodonVariation.createFrameshift(0).hashCode(), CodonVariation.createFrameshift(1).hashCode());
        Assert.assertNotEquals(CodonVariation.createFrameshift(0).hashCode(), CodonVariation.createDeletion(0).hashCode());
        Assert.assertEquals(CodonVariation.createInsertion(1, 0).hashCode(),
                CodonVariation.createInsertion(1, 0).hashCode());
        Assert.assertNotEquals(CodonVariation.createInsertion(0, 0).hashCode(),
                CodonVariation.createInsertion(1, 0).hashCode());
        Assert.assertNotEquals(CodonVariation.createInsertion(0, 0).hashCode(),
                CodonVariation.createInsertion(0, 1).hashCode());
        Assert.assertEquals(CodonVariation.createModification(1, 0).hashCode(),
                CodonVariation.createModification(1, 0).hashCode());
        Assert.assertNotEquals(CodonVariation.createModification(0, 0).hashCode(),
                CodonVariation.createModification(1, 0).hashCode());
        Assert.assertNotEquals(CodonVariation.createModification(0, 0).hashCode(),
                CodonVariation.createModification(0, 1).hashCode());
    }

    @Test
    public void testEncodingModifications() {
        // no SNVs implies no codon variations
        Assert.assertTrue(localCodonTracker.encodeSNVsAsCodons(Collections.emptyList()).isEmpty());

        // changes outside the ORF shouldn't produce codon variations
        Assert.assertTrue(localCodonTracker
                .encodeSNVsAsCodons(Collections.singletonList(new SNV(1, CALL_C, CALL_G, QUAL_30)))
                .isEmpty());
        Assert.assertTrue(localCodonTracker
                .encodeSNVsAsCodons(Collections.singletonList(new SNV(6, CALL_G, CALL_A, QUAL_30)))
                .isEmpty());
        Assert.assertTrue(localCodonTracker
                .encodeSNVsAsCodons(Collections.singletonList(new SNV(12, CALL_T, CALL_C, QUAL_30)))
                .isEmpty());

        // changes to a single codon should produce a single-codon variations
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(2, CALL_A, CALL_C, QUAL_30))),
                Collections.singletonList(CodonVariation.createModification(0, 30)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(5, CALL_C, CALL_G, QUAL_30))),
                Collections.singletonList(CodonVariation.createModification(1, 45)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(7, CALL_T, CALL_G, QUAL_30))),
                Collections.singletonList(CodonVariation.createModification(1, 25)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(8, CALL_C, CALL_A, QUAL_30))),
                Collections.singletonList(CodonVariation.createModification(1, 28)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(11, CALL_G, CALL_A, QUAL_30))),
                Collections.singletonList(CodonVariation.createModification(2, 48)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(5, CALL_C, CALL_G, QUAL_30),
                        new SNV(7, CALL_T, CALL_G, QUAL_30),
                        new SNV(8, CALL_C, CALL_A, QUAL_30))),
                Collections.singletonList(CodonVariation.createModification(1, 40)));

        // even if the change produces a nonsense codon
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(5, CALL_C, CALL_T, QUAL_30),
                        new SNV(7, CALL_T, CALL_A, QUAL_30),
                        new SNV(8, CALL_C, CALL_A, QUAL_30))),
                Collections.singletonList(CodonVariation.createModification(1, 48)));
    }

    @Test
    public void testEncodingDeletions() {
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Collections.singletonList(new SNV(2, CALL_A, NO_CALL, QUAL_30))),
                Arrays.asList(
                        CodonVariation.createFrameshift(0),
                        CodonVariation.createModification(0, 57),
                        CodonVariation.createModification(1, 55),
                        CodonVariation.createModification(2, 11)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(10, CALL_A, NO_CALL, QUAL_30),
                        new SNV(12, CALL_T, NO_CALL, QUAL_30))),
                Arrays.asList(
                        CodonVariation.createFrameshift(2),
                        CodonVariation.createModification(2, 56)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(5, CALL_C, NO_CALL, QUAL_30),
                        new SNV(7, CALL_T, NO_CALL, QUAL_30),
                        new SNV(8, CALL_C, NO_CALL, QUAL_30))),
                Collections.singletonList(CodonVariation.createDeletion(1)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(3, CALL_T, NO_CALL, QUAL_30),
                        new SNV(5, CALL_C, NO_CALL, QUAL_30),
                        new SNV(8, CALL_C, NO_CALL, QUAL_30))),
                Arrays.asList(
                        CodonVariation.createModification(0, 11),
                        CodonVariation.createDeletion(1)));

    }

    @Test
    public void testEncodingInsertions() {
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Collections.singletonList(
                        new SNV(5, NO_CALL, CALL_T, QUAL_30))),
                Arrays.asList(
                        CodonVariation.createFrameshift(1),
                        CodonVariation.createModification(1, 55),
                        CodonVariation.createModification(2, 28)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(2, NO_CALL, CALL_T, QUAL_30),
                        new SNV(2, NO_CALL, CALL_T, QUAL_30),
                        new SNV(2, NO_CALL, CALL_T, QUAL_30))),
                Collections.emptyList());
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(5, NO_CALL, CALL_T, QUAL_30),
                        new SNV(5, NO_CALL, CALL_T, QUAL_30),
                        new SNV(5, NO_CALL, CALL_T, QUAL_30))),
                Collections.singletonList(
                        CodonVariation.createInsertion(1, 63)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(4, NO_CALL, CALL_G, QUAL_30),
                        new SNV(4, NO_CALL, CALL_C, QUAL_30),
                        new SNV(4, NO_CALL, CALL_T, QUAL_30),
                        new SNV(4, NO_CALL, CALL_G, QUAL_30),
                        new SNV(4, NO_CALL, CALL_C, QUAL_30),
                        new SNV(4, NO_CALL, CALL_T, QUAL_30))),
                Arrays.asList(
                        CodonVariation.createInsertion(1, 30),
                        CodonVariation.createInsertion(1, 30)));
    }

    @Test
    public void testFrameRecoveringIndels() {
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(3, CALL_T, NO_CALL, QUAL_30),
                        new SNV(7, NO_CALL, CALL_A, QUAL_30))),
                Arrays.asList(
                        CodonVariation.createModification(0, 9),
                        CodonVariation.createModification(1, 13)));
        Assert.assertEquals(
                localCodonTracker.encodeSNVsAsCodons(Arrays.asList(
                        new SNV(3, NO_CALL, CALL_T, QUAL_30),
                        new SNV(7, CALL_T, NO_CALL, QUAL_30))),
                Arrays.asList(
                        CodonVariation.createModification(0, 15),
                        CodonVariation.createModification(1, 37)));
    }

    @Test
    public void testWildTypeCodonCounts() {
        final CodonTracker codonTracker = new CodonTracker(orfCoords, refSeq, logger);
        codonTracker.reportWildCodonCounts(new Interval(0, 17));
        codonTracker.reportWildCodonCounts(new Interval(0, 17));
        final int[] wildTypeCodonValues = codonTracker.getRefCodonValues();
        final long[][] codonCounts = codonTracker.getCodonCounts();
        for ( int codonId = 0; codonId != wildTypeCodonValues.length; ++codonId ) {
            final long[] row = codonCounts[codonId];
            final int wtValue = wildTypeCodonValues[codonId];
            for ( int codonValue = 0; codonValue != 64; ++codonValue ) {
                Assert.assertEquals(row[codonValue], codonValue==wtValue ? 2 : 0);
            }
        }
        codonTracker.reportWildCodonCounts(new Interval(6, 12));
        Assert.assertEquals(codonCounts[1][wildTypeCodonValues[1]], 2);
        Assert.assertEquals(codonCounts[2][wildTypeCodonValues[2]], 3);
    }

    @Test
    public void testVariantCodonCounts() {
        final CodonTracker codonTracker = new CodonTracker(orfCoords, refSeq, null);
        codonTracker.reportVariantCodonCounts(new Interval(0, 17), Arrays.asList(
                CodonVariation.createModification(0, 0),
                CodonVariation.createModification(1, 1),
                CodonVariation.createDeletion(1),
                CodonVariation.createInsertion(2, 2),
                CodonVariation.createFrameshift(2),
                CodonVariation.createModification(2, 6)));
        final long[][] codonCounts = codonTracker.getCodonCounts();
        Assert.assertEquals(codonCounts[0][0], 1);
        Assert.assertEquals(codonCounts[1][1], 1);
        Assert.assertEquals(codonCounts[1][CodonTracker.FRAME_PRESERVING_INDEL_INDEX], 1);
        Assert.assertEquals(codonCounts[2][6], 1);
        Assert.assertEquals(codonCounts[2][CodonTracker.FRAME_PRESERVING_INDEL_INDEX], 1);
        Assert.assertEquals(codonCounts[2][CodonTracker.FRAME_SHIFTING_INDEL_INDEX], 1);
        Assert.assertEquals(Arrays.stream(codonCounts).mapToLong(row -> Arrays.stream(row).sum()).sum(), 6L);
    }

    @Test
    public void testWildTypeCodonValuesAndExons() {
        final List<Interval> exons = CodonTracker.getExons(orfCoords, refSeq.length);
        Assert.assertEquals(exons, Arrays.asList(new Interval(2, 6), new Interval(7, 12)));
        final int[] expectedCodonValues = {14, 29, 50};
        Assert.assertTrue(Arrays.equals(CodonTracker.parseReferenceIntoCodons(refSeq, exons, null), expectedCodonValues));
    }

    @Test
    public void testFindFrameShift() {
        Assert.assertEquals(
                localCodonTracker.findFrameShift(Collections.singletonList(new SNV(1, CALL_A, NO_CALL, QUAL_30))),
                NO_FRAME_SHIFT_CODON);
        Assert.assertEquals(
                localCodonTracker.findFrameShift(Collections.singletonList(new SNV(2, CALL_A, NO_CALL, QUAL_30))),
                0);
        Assert.assertEquals(
                localCodonTracker.findFrameShift(Collections.singletonList(new SNV(4, NO_CALL, CALL_C, QUAL_30))),
                0);
        Assert.assertEquals(
                localCodonTracker.findFrameShift(Collections.singletonList(new SNV(6, CALL_G, NO_CALL, QUAL_30))),
                NO_FRAME_SHIFT_CODON);
        Assert.assertEquals(
                localCodonTracker.findFrameShift(Collections.singletonList(new SNV(12, NO_CALL, CALL_C, QUAL_30))),
                NO_FRAME_SHIFT_CODON);
        Assert.assertEquals(
                localCodonTracker.findFrameShift(Arrays.asList(
                        new SNV(5, CALL_C, NO_CALL, QUAL_30),
                        new SNV(6, CALL_G, NO_CALL, QUAL_30),
                        new SNV(7, CALL_T, NO_CALL, QUAL_30),
                        new SNV(8, CALL_C, NO_CALL, QUAL_30))),
                NO_FRAME_SHIFT_CODON);
    }

    @Test
    public void testIsStop() {
        Assert.assertFalse(CodonTracker.isStop(-1));
        Assert.assertFalse(CodonTracker.isStop(47));
        Assert.assertTrue(CodonTracker.isStop(48));
        Assert.assertFalse(CodonTracker.isStop(49));
        Assert.assertTrue(CodonTracker.isStop(50));
        Assert.assertTrue(CodonTracker.isStop(56));
        Assert.assertFalse(CodonTracker.isStop(57));
        Assert.assertFalse(CodonTracker.isStop(99));
    }

    @Test
    public void testIsExonic() {
        Assert.assertFalse(localCodonTracker.isExonic(-1));
        Assert.assertFalse(localCodonTracker.isExonic(0));
        Assert.assertFalse(localCodonTracker.isExonic(1));
        Assert.assertTrue(localCodonTracker.isExonic(2));
        Assert.assertTrue(localCodonTracker.isExonic(3));
        Assert.assertTrue(localCodonTracker.isExonic(4));
        Assert.assertTrue(localCodonTracker.isExonic(5));
        Assert.assertFalse(localCodonTracker.isExonic(6));
        Assert.assertTrue(localCodonTracker.isExonic(7));
        Assert.assertTrue(localCodonTracker.isExonic(8));
        Assert.assertTrue(localCodonTracker.isExonic(9));
        Assert.assertTrue(localCodonTracker.isExonic(10));
        Assert.assertTrue(localCodonTracker.isExonic(11));
        Assert.assertFalse(localCodonTracker.isExonic(12));
        Assert.assertFalse(localCodonTracker.isExonic(99));
    }

    @Test
    public void testExonicBaseCount() {
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(-1), 0);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(0), 0);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(1), 0);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(2), 0);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(3), 1);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(4), 2);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(5), 3);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(6), 4);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(7), 4);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(8), 5);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(9), 6);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(10), 7);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(11), 8);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(12), 9);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(13), 9);
        Assert.assertEquals(localCodonTracker.exonicBaseIndex(99), 9);
    }

    @Test
    public void testReadReportFlanks() {
        final List<SNV> snvList = Collections.singletonList(new SNV(8, CALL_C, CALL_A, QUAL_30));
        final ReadReport readReport1 =
                new ReadReport(Collections.singletonList(new Interval(2,15)), snvList);
        Assert.assertTrue(readReport1.hasCleanLeftFlank(6));
        Assert.assertFalse(readReport1.hasCleanLeftFlank(7));
        Assert.assertTrue(readReport1.hasCleanRightFlank(6, refSeq.length));
        Assert.assertFalse(readReport1.hasCleanRightFlank(7, refSeq.length));

        // unpenalized for going off the end of the reference
        final ReadReport readReport2 =
                new ReadReport(Collections.singletonList(new Interval(0,17)), snvList);
        Assert.assertTrue(readReport2.hasCleanLeftFlank(9));
        Assert.assertTrue(readReport2.hasCleanRightFlank(9, refSeq.length));
    }

    @Test
    public void testCombiningCoverage() {
        final ReadReport report1 = new ReadReport(Arrays.asList(
                new Interval(0, 9), new Interval(20, 40), new Interval(60,70), new Interval(100, 120)),
                Collections.emptyList());
        final ReadReport report2 = new ReadReport(Arrays.asList(
                new Interval(17, 21), new Interval(22, 27), new Interval(38, 50), new Interval(58, 60), new Interval(70, 72)),
                Collections.emptyList());
        Assert.assertEquals(new ReadReport(report1, report2).getRefCoverage(),
                Arrays.asList(new Interval(0, 9), new Interval(17, 50), new Interval(58, 72), new Interval(100, 120)));
    }

    @Test
    public void testCombiningSNVLists() {
        final List<SNV> list1 = Arrays.asList(
                new SNV(10, CALL_A, CALL_C, QUAL_30),
                new SNV(20, CALL_C, CALL_A, QUAL_30),
                new SNV(30, CALL_C, CALL_G, QUAL_30));
        final List<SNV> list2 = Arrays.asList(
                new SNV(35, CALL_T, CALL_C, QUAL_30),
                new SNV(45, CALL_G, CALL_A, QUAL_30),
                new SNV(55, CALL_C, CALL_A, QUAL_30));
        final ReadReport report1 = new ReadReport(Collections.singletonList(new Interval(0, 32)), list1);
        final ReadReport report2 = new ReadReport(Collections.singletonList(new Interval(33, 60)), list2);
        final List<SNV> combinedList = new ArrayList<>(list1);
        combinedList.addAll(list2);
        Assert.assertEquals(new ReadReport(report1, report2).getVariations(), combinedList);

        // region of overlap has no conflicting SNVs
        final ReadReport report3 = new ReadReport(Collections.singletonList(new Interval(0, 35)), list1);
        final ReadReport report4 = new ReadReport(Collections.singletonList(new Interval(33, 60)), list2);
        Assert.assertEquals(new ReadReport(report3, report4).getVariations(), combinedList);

        // region of overlap has conflicting SNV: SNV at 35 is missing from report5
        final ReadReport report5 = new ReadReport(Collections.singletonList(new Interval(0, 36)), list1);
        final ReadReport report6 = new ReadReport(Collections.singletonList(new Interval(33, 60)), list2);
        Assert.assertNull(new ReadReport(report5, report6).getVariations());

        // repair flaw above
        final List<SNV> list1Plus = new ArrayList<>(list1);
        list1Plus.add(new SNV(35, CALL_T, CALL_C, QUAL_10));
        final ReadReport report7 = new ReadReport(Collections.singletonList(new Interval(0, 36)), list1Plus);
        Assert.assertEquals(new ReadReport(report7, report6).getVariations(), combinedList);

        // now mess it up again by making the SNVs unequal
        list1Plus.set(3, new SNV(35, CALL_T, CALL_A, QUAL_30));
        Assert.assertNull(new ReadReport(report7, report6).getVariations());
    }

    @Test
    public void testPairedApplication() {
        final byte[] longRefSeq = "TTTTTTTTTTAAAAAAAAAACCCCCCCCCCCCCCCTTTTTTTTTTGGGGGGGGGGCCCCCCCCCC".getBytes();
        reference = new Reference(longRefSeq);
        codonTracker = new CodonTracker("5-25", longRefSeq, logger);

        final List<SNV> list1 = Arrays.asList(
                new SNV(10, CALL_A, CALL_C, QUAL_30),
                new SNV(20, CALL_C, CALL_A, QUAL_30),
                new SNV(30, CALL_C, CALL_G, QUAL_30));
        final List<SNV> list2 = Arrays.asList(
                new SNV(35, CALL_T, CALL_C, QUAL_30),
                new SNV(45, CALL_G, CALL_A, QUAL_30),
                new SNV(55, CALL_C, CALL_A, QUAL_30));

        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final List<GATKRead> reads =
            ArtificialReadUtils.createPair(header, "blahBlah", 150, 1, 151, true, false);
        final GATKRead read1 = reads.get(0);
        final GATKRead read2 = reads.get(1);
        final ReadReport report1 = new ReadReport(Collections.singletonList(new Interval(0, 32)), list1);
        final ReadReport report2 = new ReadReport(Collections.singletonList(new Interval(33, 60)), list2);
        updateCountsForPair(read1, report1, read2, report2);

        // shouldn't be applied -- fails flanking bases test
        Assert.assertEquals(reference.countSpanners(0, 32), 0);
        Assert.assertEquals(reference.countSpanners(32, 33), 0);
        Assert.assertEquals(reference.countSpanners(33, 60), 0);

        // only read1 should be applied
        minFlankingLength = 1;
        updateCountsForPair(read1, report1, read2, report2);
        Assert.assertEquals(reference.countSpanners(0, 32), 1);
        Assert.assertEquals(reference.countSpanners(32, 33), 0);
        Assert.assertEquals(reference.countSpanners(33, 60), 0);

        // should be applied as one
        final ReadReport report3 = new ReadReport(Collections.singletonList(new Interval(0, 35)), list1);
        final ReadReport report4 = new ReadReport(Collections.singletonList(new Interval(33, 60)), list2);
        updateCountsForPair(read1, report3, read2, report4);
        Assert.assertEquals(reference.countSpanners(0, 32), 2);
        Assert.assertEquals(reference.countSpanners(32, 33), 1);
        Assert.assertEquals(reference.countSpanners(33, 60), 1);
    }

    @Test
    public void testCalculateTrim() {
        final byte BAD_Q = (byte)(minQ - 1);
        final byte GOOD_Q = (byte)minQ;
        final byte[] quals = new byte[100];
        Arrays.fill(quals, BAD_Q);
        Assert.assertEquals(calculateQualityTrim(quals).size(), 0);
        Arrays.fill(quals, GOOD_Q);
        Assert.assertEquals(calculateQualityTrim(quals), new Interval(0, 100));
        quals[minLength] = BAD_Q;
        Assert.assertEquals(calculateQualityTrim(quals), new Interval(0, 100));
        quals[minLength - 1] = BAD_Q;
        Assert.assertEquals(calculateQualityTrim(quals), new Interval(minLength + 1, 100));
        quals[quals.length - minLength] = BAD_Q;
        Assert.assertEquals(calculateQualityTrim(quals), new Interval(minLength + 1, quals.length - minLength));
        Arrays.fill(quals, GOOD_Q);
        for ( int idx = minLength - 1; idx < quals.length; idx += minLength ) {
            quals[idx] = BAD_Q;
        }
        Assert.assertEquals(calculateQualityTrim(quals).size(), 0);
    }

    @Test
    public void testGetReadReport() {
        final SAMFileHeader header =
                ArtificialReadUtils.createArtificialSamHeader(1, 0, 17);
        final byte[] bases = Arrays.copyOf(refSeq, refSeq.length);
        final byte[] quals = new byte[17];
        Arrays.fill(quals, QUAL_30);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1,
                                                                          bases, quals, "17M");
        read.setMappingQuality(0);
        reference = new Reference(refSeq);
        Assert.assertEquals(getReadReport(read).getRefCoverage(), Collections.emptyList());

        read.setMappingQuality(60);
        final ReadReport readReport = getReadReport(read);
        Assert.assertEquals(readReport.getRefCoverage(), Collections.singletonList(new Interval(0, 17)));
        Assert.assertEquals(readReport.getVariations().size(), 0);
        read.getBasesNoCopy()[3] = CALL_A;
        final ReadReport readReport2 = getReadReport(read);
        Assert.assertEquals(readReport2.getVariations(),
                Collections.singletonList(new SNV(3, CALL_T, CALL_A, QUAL_30)));
    }

    @Test
    public void testFindLargeDeletions() {
        AnalyzeSaturationMutagenesis.minAltLength = 2;
        AnalyzeSaturationMutagenesis.findLargeDels = true;
        final SAMFileHeader header =
                ArtificialReadUtils.createArtificialSamHeader(1, 0, 17);
        final byte[] bases = new byte[10];
        System.arraycopy(refSeq, 0, bases, 0, 5);
        System.arraycopy(refSeq, 12, bases, 5, 5);
        final byte[] quals = new byte[10];
        Arrays.fill(quals, QUAL_30);

        // OK: primary preceeds supplementary
        final GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1,
                bases, quals, "5M5S");
        read1.setAttribute("SA", "0,13,+,5S5M,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read1).toString(), "5M7D5M");

        // OK: supplementary preceeds primary
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 13,
                bases, quals, "5S5M");
        read2.setAttribute("SA", "0,1,+,5M5S,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read2).toString(), "5M7D5M");

        // too much overlap on read
        final GATKRead read3 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1,
                bases, quals, "6M4S");
        read3.setAttribute("SA", "0,11,+,3S7M,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read3).toString(), "6M4S");

        // too far apart on read
        final GATKRead read4 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1,
                bases, quals, "3M7S");
        read4.setAttribute("SA", "0,14,+,6S4M,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read4).toString(), "3M7S");

        // wrong strand
        final GATKRead read5 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1,
                bases, quals, "5M5S");
        read5.setAttribute("SA", "0,13,-,5S5M,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read5).toString(), "5M5S");

        // wrong order
        final GATKRead read6 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 13,
                bases, quals, "5M5S");
        read6.setAttribute("SA", "0,1,+,5S5M,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read6).toString(), "5M5S");

        // OK: 1-base overlap on read
        final GATKRead read7 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1,
                bases, quals, "6M4S");
        read7.setAttribute("SA", "0,13,+,5S5M,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read7).toString(), "6M7D4M");

        // OK: 1-base separation on read
        final GATKRead read8 = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 1,
                bases, quals, "4M6S");
        read8.setAttribute("SA", "0,13,+,5S5M,60,0;");
        Assert.assertEquals(AnalyzeSaturationMutagenesis.ReadReport.findLargeDeletions(read8).toString(), "4M7D6M");
    }
}
