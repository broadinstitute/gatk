package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.tools.dragstr.BufferedReferenceBases;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CollectSTRSlippageCountsUnitTest extends GATKBaseTest {
    private static final int MIN_LENGTH = CollectSTRSlippageCounts.DEFAULT_MIN_REPEAT_LENGTH;
    private static final int MIN_UNITS = CollectSTRSlippageCounts.DEFAULT_MIN_REPEAT_UNITS;
    private static final int MAX_UNIT_LENGTH = CollectSTRSlippageCounts.DEFAULT_MAX_REPEAT_UNIT_LENGTH;
    private static final int FLANK = CollectSTRSlippageCounts.DEFAULT_FLANKING_BASES;

    @DataProvider
    public Object[][] repeatAtData() {
        // bases, start, max unit length, min units, min length, expected unit (null if none), expected units
        return new Object[][] {
                {"AAAAAAAT", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "A", 7},
                {"AAAAAAT", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, null, 0},      // 6 bases is too short
                {"GAAAAAAAT", 2, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "A", 7},    // starting partway into the contig
                {"GAAAAAAA", 2, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "A", 7},     // ending the contig
                {"aaaaaaaT", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "A", 7},     // lower case, as in soft-masked references
                {"AAAAAAAAAAT", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "A", 10}, // AA and AAA are not units
                {"ACACACACG", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "AC", 4},
                {"ACACACACAG", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "AC", 4},  // a partial unit at the end is left out
                {"ACACACACACACT", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "AC", 6},   // ACAC is not a unit
                {"ACACACG", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, null, 0},     // 6 bases is too short
                {"ACACACG", 1, MAX_UNIT_LENGTH, MIN_UNITS, 6, "AC", 3},
                {"ACGACGACGT", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "ACG", 3},
                {"ACGTAACGTAC", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, "ACGTA", 2},
                {"ACGTAACGTAC", 1, 4, MIN_UNITS, MIN_LENGTH, null, 0},               // the unit is longer than the maximum
                {"ACGTACGTAC", 1, MAX_UNIT_LENGTH, 3, MIN_LENGTH, null, 0},          // too few units
                {"NNNNNNNNNN", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, null, 0},  // N is not a base
                {"ACNACNACNACN", 1, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH, null, 0},    // so units can't contain it
        };
    }

    @Test(dataProvider = "repeatAtData")
    public void testRepeatAt(final String bases, final int start, final int maxUnitLength, final int minUnits,
                             final int minLength, final String expectedUnit, final int expectedUnits) {
        final CollectSTRSlippageCounts.Repeat repeat = makeScanner(bases, 1000, maxUnitLength, minUnits, minLength).repeatAt(start);
        if (expectedUnit == null) {
            Assert.assertNull(repeat);
        } else {
            Assert.assertEquals(repeat.getStart(), start);
            Assert.assertEquals(repeat.getUnit(), expectedUnit);
            Assert.assertEquals(repeat.getUnits(), expectedUnits);
        }
    }

    @Test
    public void testScannerFindsNonOverlappingRepeatsAcrossBuffers() {
        // A x 8 at 3-10, AC x 4 at 13-20 (the A at 21 is a partial unit), A x 8 at 21-28
        final String contig = "GT" + "AAAAAAAA" + "CT" + "ACACACAC" + "aaaaaaaa" + "GCTTGA";
        final List<String> expected = List.of("3 A 8", "13 AC 4", "21 A 8");

        // a buffer of 4 bases makes runs of copies cross buffer boundaries, and scanning in steps must give the same result
        for (final int[] steps : new int[][] { {100}, {5, 15, 100}, {1, 2, 3, 4, 20, 21, 34} }) {
            final CollectSTRSlippageCounts.RepeatScanner scanner = makeScanner(contig, 4, MAX_UNIT_LENGTH, MIN_UNITS, MIN_LENGTH);
            final List<String> found = new ArrayList<>();
            for (final int step : steps) {
                scanner.scanThrough(step, repeat -> found.add(repeat.getStart() + " " + repeat.getUnit() + " " + repeat.getUnits()));
            }
            Assert.assertEquals(found, expected);
        }
    }

    @DataProvider
    public Object[][] addReadData() {
        // a homopolymer A x 10 at 101-110 and a dinucleotide AC x 5 at 101-110, so the flanks end at 96 and 115
        // unit, read start, cigar, whether the fragment counts, expected bin
        return new Object[][] {
                {"A", 91, "50M", true, 4},
                {"A", 96, "50M", true, 4},          // just reaches the left flank
                {"A", 97, "50M", false, -1},        // starts too late to flank the STR
                {"A", 66, "50M", true, 4},          // ends exactly at the right flank (66-115)
                {"A", 91, "12M2D38M", true, 2},     // -2 units
                {"A", 91, "12M1I37M", true, 5},     // +1 unit
                {"A", 91, "12M6D38M", true, 0},     // at most -4 units
                {"A", 91, "12M5I33M", true, 8},     // at least +4 units
                {"A", 91, "8M1D41M", true, 3},      // a deletion in the left flank counts too
                {"A", 91, "5M3D45M", false, -1},    // the left flank base is deleted
                {"A", 91, "5S45M", true, 4},        // soft clips elsewhere don't matter
                {"A", 97, "6S44M", false, -1},      // the left flank base is soft-clipped
                {"AC", 91, "12M2D38M", true, 3},    // -1 unit
                {"AC", 91, "12M3D38M", true, 9},    // not whole units
        };
    }

    @Test(dataProvider = "addReadData")
    public void testAddRead(final String unit, final int readStart, final String cigar, final boolean counts, final int expectedBin) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        final CollectSTRSlippageCounts.Repeat repeat = new CollectSTRSlippageCounts.Repeat("1", 101, unit, 10 / unit.length());
        final GATKRead read = makeRead(header, "read", readStart, cigar);

        Assert.assertEquals(repeat.addRead(read, FLANK), counts);
        final int[] expectedCounts = new int[CollectSTRSlippageCounts.COUNT_COLUMNS.size()];
        if (counts) {
            expectedCounts[expectedBin] = 1;
        }
        Assert.assertEquals(repeat.getCounts(), expectedCounts);
    }

    @Test
    public void testFragmentCountsOnce() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        final CollectSTRSlippageCounts.Repeat repeat = new CollectSTRSlippageCounts.Repeat("1", 101, "A", 10);
        Assert.assertTrue(repeat.addRead(makeRead(header, "fragment", 91, "50M"), FLANK));
        Assert.assertFalse(repeat.addRead(makeRead(header, "fragment", 95, "12M1I37M"), FLANK));    // its mate
        Assert.assertTrue(repeat.addRead(makeRead(header, "other", 91, "12M1I37M"), FLANK));
        Assert.assertEquals(repeat.getFragmentCount(), 2);
        Assert.assertEquals(repeat.getCounts()[4], 1);
        Assert.assertEquals(repeat.getCounts()[5], 1);
    }

    @Test
    public void testBin() {
        Assert.assertEquals(CollectSTRSlippageCounts.bin(0, 1), 4);
        Assert.assertEquals(CollectSTRSlippageCounts.bin(-3, 1), 1);
        Assert.assertEquals(CollectSTRSlippageCounts.bin(-9, 1), 0);
        Assert.assertEquals(CollectSTRSlippageCounts.bin(9, 1), 8);
        Assert.assertEquals(CollectSTRSlippageCounts.bin(-4, 2), 2);
        Assert.assertEquals(CollectSTRSlippageCounts.bin(6, 3), 6);
        Assert.assertEquals(CollectSTRSlippageCounts.bin(-3, 2), CollectSTRSlippageCounts.NOT_WHOLE_UNITS_BIN);
        Assert.assertEquals(CollectSTRSlippageCounts.bin(5, 4), CollectSTRSlippageCounts.NOT_WHOLE_UNITS_BIN);
    }

    private static CollectSTRSlippageCounts.RepeatScanner makeScanner(final String bases, final int bufferSize, final int maxUnitLength,
                                                                       final int minUnits, final int minLength) {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(new SAMSequenceRecord("1", bases.length())));
        final ReferenceMemorySource reference = new ReferenceMemorySource(
                new ReferenceBases(bases.getBytes(), new SimpleInterval("1", 1, bases.length())), dictionary);
        return new CollectSTRSlippageCounts.RepeatScanner(BufferedReferenceBases.of(reference, "1", bufferSize),
                maxUnitLength, minUnits, minLength);
    }

    private static GATKRead makeRead(final SAMFileHeader header, final String name, final int start, final String cigar) {
        final int length = htsjdk.samtools.TextCigarCodec.decode(cigar).getReadLength();
        final byte[] bases = new byte[length];
        final byte[] quals = new byte[length];
        Arrays.fill(bases, (byte) 'A');
        Arrays.fill(quals, (byte) 30);
        return ArtificialReadUtils.createArtificialRead(header, name, 0, start, bases, quals, cigar);
    }
}
