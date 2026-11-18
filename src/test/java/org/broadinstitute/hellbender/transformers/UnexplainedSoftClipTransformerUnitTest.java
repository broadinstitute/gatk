package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class UnexplainedSoftClipTransformerUnitTest {
    private static final int CONTIG_LENGTH = 2000;
    private static final int ALIGNED_LENGTH = 60;
    private static final byte GOOD_QUALITY = 30;
    private static final int MIN_TAIL_QUALITY = 9;

    // (CA)10 at positions 1201-1220, flanked by 50 bases of G and T that a clip of C and A can't match
    private static final int REPEAT_START = 1201;
    private static final String REPEAT = "CACACACACACACACACACA";
    private static final int REPEAT_END = REPEAT_START + REPEAT.length();    // exclusive
    // mostly A, but too interrupted to be a tandem repeat
    private static final int MOSTLY_A_START = 1601;
    private static final String MOSTLY_A = "AAAACAAAAGAAAATAAAAC";

    // novel sequence inserted after INSERTION_POINT, for the tests of stage 2
    private static final int INSERTION_POINT = 700;
    private static final int INSERTION_LENGTH = 30;
    // an aligned deletion, for the tests of stage 2
    private static final int DELETION_START = 1001;
    private static final int DELETION_LENGTH = 8;

    private final RandomDNA random = new RandomDNA(13);
    // two contigs, but reference bases only for the first: reads on the second have no soft clips
    private final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, CONTIG_LENGTH);
    private final byte[] reference = makeReference();
    private final byte[] inserted = new RandomDNA(17).nextBases(INSERTION_LENGTH);
    private final ReferenceDataSource referenceSource = ReferenceDataSource.of(
            new ReferenceBases(reference, new SimpleInterval(header.getSequence(0).getSequenceName(), 1, CONTIG_LENGTH)),
            header.getSequenceDictionary());

    private byte[] makeReference() {
        final byte[] result = random.nextBases(CONTIG_LENGTH);
        put(result, REPEAT_START - 50, randomGT(50));
        put(result, REPEAT_START, REPEAT.getBytes());
        put(result, REPEAT_START + REPEAT.length(), randomGT(50));
        put(result, MOSTLY_A_START, MOSTLY_A.getBytes());
        return result;
    }

    @DataProvider(name = "clips")
    public Object[][] clips() {
        final List<Object[]> rows = new ArrayList<>();
        for (final boolean right : new boolean[] {true, false}) {
            // the clip point is between the last aligned base and the junction, the reference position of the first clipped base
            final int junction = 501;
            rows.add(new Object[] {"continues the reference", right, junction, referenceClip(junction, right, 0, 10), GOOD_QUALITY, true});
            rows.add(new Object[] {"continues the reference with a mismatch", right, junction, withMismatches(referenceClip(junction, right, 0, 10), 4), GOOD_QUALITY, true});
            rows.add(new Object[] {"short clip continues the reference with a mismatch", right, junction, withMismatches(referenceClip(junction, right, 0, 3), 1), GOOD_QUALITY, true});
            rows.add(new Object[] {"too many mismatches", right, junction, withMismatches(referenceClip(junction, right, 0, 10), 2, 6), GOOD_QUALITY, false});
            rows.add(new Object[] {"two mismatches in twenty bases", right, junction, withMismatches(referenceClip(junction, right, 0, 20), 2, 12), GOOD_QUALITY, true});
            rows.add(new Object[] {"7-base deletion", right, junction, referenceClip(junction, right, 7, 19), GOOD_QUALITY, true});
            rows.add(new Object[] {"5-base duplication", right, junction, referenceClip(junction, right, -5, 7), GOOD_QUALITY, true});
            rows.add(new Object[] {"3-base insertion", right, junction, Utils.concat("GTA".getBytes(), referenceClip(junction, right, 0, 10)), GOOD_QUALITY, true});
            rows.add(new Object[] {"deletion too long", right, junction, referenceClip(junction, right, 60, 72), GOOD_QUALITY, false});
            rows.add(new Object[] {"unexplained", right, junction, random.nextBases(20), GOOD_QUALITY, false});
            rows.add(new Object[] {"low quality", right, junction, referenceClip(junction, right, 0, 10), (byte) 15, false});

            final int repeatJunction = REPEAT_START + 10;
            rows.add(new Object[] {"stays in a repeat", right, repeatJunction, referenceClip(repeatJunction, right, 0, 10), GOOD_QUALITY, false});
            rows.add(new Object[] {"spans a repeat", right, repeatJunction, referenceClip(repeatJunction, right, 0, 20), GOOD_QUALITY, true});
            rows.add(new Object[] {"short clip in a repeat", right, repeatJunction, referenceClip(repeatJunction, right, 0, 6), GOOD_QUALITY, true});

            // a short clip of the repeat's flank, after the read skips the rest of the repeat, or of an added unit and the flank
            final int pastRepeat = right ? REPEAT_END - repeatJunction : repeatJunction - (REPEAT_START - 1);
            rows.add(new Object[] {"short clip past deleted repeat units", right, repeatJunction, referenceClip(repeatJunction, right, pastRepeat, pastRepeat + 4), GOOD_QUALITY, true});
            rows.add(new Object[] {"short clip past deleted repeat units with a mismatch", right, repeatJunction, withMismatches(referenceClip(repeatJunction, right, pastRepeat, pastRepeat + 4), 3), GOOD_QUALITY, false});
            final int repeatEndJunction = right ? REPEAT_END : REPEAT_START - 1;
            rows.add(new Object[] {"short clip of an added repeat unit", right, repeatEndJunction, referenceClip(repeatEndJunction, right, -2, 3), GOOD_QUALITY, true});
            rows.add(new Object[] {"short clip of an added repeat unit with a mismatch", right, repeatEndJunction, withMismatches(referenceClip(repeatEndJunction, right, -2, 3), 4), GOOD_QUALITY, false});
            rows.add(new Object[] {"short clip in a repeat out of phase", right, repeatJunction, referenceClip(repeatJunction, right, -1, 5), GOOD_QUALITY, false});

            final int mostlyAJunction = right ? MOSTLY_A_START : MOSTLY_A_START + MOSTLY_A.length() - 1;
            rows.add(new Object[] {"mostly one base", right, mostlyAJunction, referenceClip(mostlyAJunction, right, 0, 20), GOOD_QUALITY, false});
        }
        return rows.toArray(new Object[][] {});
    }

    @Test(dataProvider = "clips")
    public void testOneClip(final String description, final boolean right, final int junction, final byte[] orientedClip,
                            final byte clipQuality, final boolean kept) {
        final byte[] orientedQuals = new byte[orientedClip.length];
        Arrays.fill(orientedQuals, clipQuality);
        checkClip(description, right, junction, orientedClip, orientedQuals, kept);
    }

    // the transformer judges a clip without its tail of bases with quality at most MIN_TAIL_QUALITY
    @DataProvider(name = "tails")
    public Object[][] tails() {
        final List<Object[]> rows = new ArrayList<>();
        for (final boolean right : new boolean[] {true, false}) {
            final int junction = 501;
            final byte[] tail = random.nextBases(10);
            rows.add(new Object[] {"explained, then a low-quality tail", right, junction, referenceClip(junction, right, 0, 4), tail, true});
            rows.add(new Object[] {"unexplained, then a low-quality tail", right, junction, withMismatches(referenceClip(junction, right, 0, 4), 0, 2), tail, false});
            rows.add(new Object[] {"low-quality tail only", right, junction, new byte[0], tail, false});
        }
        return rows.toArray(new Object[][] {});
    }

    @Test(dataProvider = "tails")
    public void testLowQualityTail(final String description, final boolean right, final int junction, final byte[] orientedGoodBases,
                                   final byte[] orientedTail, final boolean kept) {
        final byte[] orientedQuals = new byte[orientedGoodBases.length + orientedTail.length];
        Arrays.fill(orientedQuals, (byte) MIN_TAIL_QUALITY);
        Arrays.fill(orientedQuals, 0, orientedGoodBases.length, GOOD_QUALITY);
        checkClip(description, right, junction, Utils.concat(orientedGoodBases, orientedTail), orientedQuals, kept);
    }

    // make a read with ALIGNED_LENGTH aligned bases ending (or, for a left clip, starting) next to the junction and a soft
    // clip of the given bases and qualities, oriented away from the aligned bases; check whether the transformer keeps it
    private void checkClip(final String description, final boolean right, final int junction, final byte[] orientedClip,
                           final byte[] orientedQuals, final boolean kept) {
        final int alignmentStart = right ? junction - ALIGNED_LENGTH : junction + 1;
        final byte[] aligned = referenceBases(alignmentStart, ALIGNED_LENGTH);
        final byte[] clip = orientedClip.clone();
        final byte[] clipQuals = orientedQuals.clone();
        if (!right) {
            ArrayUtils.reverse(clip);
            ArrayUtils.reverse(clipQuals);
        }
        final byte[] bases = right ? Utils.concat(aligned, clip) : Utils.concat(clip, aligned);
        final byte[] alignedQuals = new byte[ALIGNED_LENGTH];
        Arrays.fill(alignedQuals, GOOD_QUALITY);
        final byte[] quals = right ? Utils.concat(alignedQuals, clipQuals) : Utils.concat(clipQuals, alignedQuals);
        final String cigar = right ? ALIGNED_LENGTH + "M" + clip.length + "S" : clip.length + "S" + ALIGNED_LENGTH + "M";

        final GATKRead transformed = transformOne(ArtificialReadUtils.createArtificialRead(header, "read", 0, alignmentStart, bases, quals, cigar));
        Assert.assertEquals(transformed.getCigar().toString(), kept ? cigar : cigar.replace('S', 'H'), description);
        Assert.assertEquals(transformed.getStart(), alignmentStart, description);
        Assert.assertEquals(transformed.getBases(), kept ? bases : aligned, description);
    }

    @Test
    public void testNoSoftClips() {
        final GATKRead read = makeRead(801, new byte[0], "60M", new byte[0]);
        Assert.assertSame(transformOne(read), read);
    }

    @Test
    public void testBothEnds() {
        final int start = 801;
        final byte[] aligned = referenceBases(start, ALIGNED_LENGTH);
        final byte[] unexplained = random.nextBases(10);
        final byte[] continuation = referenceBases(start + ALIGNED_LENGTH, 10);

        final GATKRead oneExplained = transformOne(makeRead(start, unexplained, "10S60M10S", continuation));
        Assert.assertEquals(oneExplained.getCigar().toString(), "10H60M10S");
        Assert.assertEquals(oneExplained.getBases(), Utils.concat(aligned, continuation));

        final GATKRead neitherExplained = transformOne(makeRead(start, unexplained, "10S60M10S", random.nextBases(10)));
        Assert.assertEquals(neitherExplained.getCigar().toString(), "10H60M10H");
        Assert.assertEquals(neitherExplained.getBases(), aligned);
    }

    @Test
    public void testSoftClipInsideHardClip() {
        final GATKRead transformed = transformOne(makeRead(801, random.nextBases(10), "5H10S60M", new byte[0]));
        Assert.assertEquals(transformed.getCigar().toString(), "15H60M");
        Assert.assertEquals(transformed.getStart(), 801);
    }

    // the reference ends just past the clip point, so it can't explain the clip
    @Test
    public void testContigEnds() {
        final GATKRead atStart = transformOne(makeRead(1, random.nextBases(10), "10S60M", new byte[0]));
        Assert.assertEquals(atStart.getCigar().toString(), "10H60M");

        final GATKRead atEnd = transformOne(makeRead(CONTIG_LENGTH - ALIGNED_LENGTH + 1, new byte[0], "60M10S", random.nextBases(10)));
        Assert.assertEquals(atEnd.getCigar().toString(), "60M10H");
    }

    // ------------------------------------------------------------------------------------------------ stage 2

    @Test
    public void testClipsFromBothSidesOfAnInsertion() {
        // the reference can't explain either clip of novel inserted sequence alone
        assertHardClipped(transformOne(rightClipped("right", INSERTION_POINT, insertedBases(0, 20))));
        assertHardClipped(transformOne(leftClipped("left", INSERTION_POINT + 1, insertedBases(8, 30))));

        // together they share 12 clip bases
        transform(rightClipped("right", INSERTION_POINT, insertedBases(0, 20)), leftClipped("left", INSERTION_POINT + 1, insertedBases(8, 30)))
                .forEach(UnexplainedSoftClipTransformerUnitTest::assertSoftClipped);

        // 9 shared clip bases are too few
        transform(rightClipped("right", INSERTION_POINT, insertedBases(0, 20)), leftClipped("left", INSERTION_POINT + 1, insertedBases(11, 30)))
                .forEach(UnexplainedSoftClipTransformerUnitTest::assertHardClipped);

        // the aligner may align a few more bases, which happen to match the reference, before clipping
        final byte[] shiftedClip = Utils.concat(insertedBases(8, 30), referenceBases(INSERTION_POINT + 1, 3));
        transform(rightClipped("right", INSERTION_POINT, insertedBases(0, 20)), leftClipped("left", INSERTION_POINT + 4, shiftedClip))
                .forEach(UnexplainedSoftClipTransformerUnitTest::assertSoftClipped);

        // a read's mate doesn't count
        transform(rightClipped("pair", INSERTION_POINT, insertedBases(0, 20)), leftClipped("pair", INSERTION_POINT + 1, insertedBases(8, 30)))
                .forEach(UnexplainedSoftClipTransformerUnitTest::assertHardClipped);
    }

    // chimeric junk and reads from paralogs agree on one side too
    @Test
    public void testClipsFromOneSideDontCount() {
        transform(rightClipped("a", INSERTION_POINT, insertedBases(0, 20)), rightClipped("b", INSERTION_POINT, insertedBases(0, 25)))
                .forEach(UnexplainedSoftClipTransformerUnitTest::assertHardClipped);
    }

    @Test
    public void testClipsMatchingAnAlignedInsertion() {
        assertHardClipped(transformOne(rightClipped("right", INSERTION_POINT, insertedBases(0, 4))));
        assertHardClipped(transformOne(leftClipped("left", INSERTION_POINT + 1, insertedBases(6, 10))));

        final List<GATKRead> reads = transform(rightClipped("right", INSERTION_POINT, insertedBases(0, 4)),
                withInsertion("aligned", insertedBases(0, 10)), leftClipped("left", INSERTION_POINT + 1, insertedBases(6, 10)));
        Assert.assertEquals(reads.stream().map(GATKRead::getName).toList(), List.of("right", "aligned", "left"));
        assertSoftClipped(reads.get(0));
        assertSoftClipped(reads.get(2));

        // a clip with a mismatch, a clip of 2 bases, and a clip of the read's own mate don't count
        final byte[] mismatched = insertedBases(0, 4);
        mismatched[3] = (byte) (mismatched[3] == 'A' ? 'C' : 'A');
        assertHardClipped(transform(rightClipped("right", INSERTION_POINT, mismatched), withInsertion("aligned", insertedBases(0, 10))).get(0));
        assertHardClipped(transform(rightClipped("right", INSERTION_POINT, insertedBases(0, 2)), withInsertion("aligned", insertedBases(0, 10))).get(0));
        assertHardClipped(transform(rightClipped("aligned", INSERTION_POINT, insertedBases(0, 4)), withInsertion("aligned", insertedBases(0, 10))).get(0));
    }

    @Test
    public void testClipsMatchingAnAlignedDeletion() {
        final int lastBefore = DELETION_START - 1;
        final int firstAfter = DELETION_START + DELETION_LENGTH;
        assertHardClipped(transformOne(rightClipped("right", lastBefore, referenceBases(firstAfter, 4))));
        assertHardClipped(transformOne(leftClipped("left", firstAfter, referenceBases(lastBefore - 3, 4))));

        final List<GATKRead> reads = transform(rightClipped("right", lastBefore, referenceBases(firstAfter, 4)),
                withDeletion("aligned"), leftClipped("left", firstAfter, referenceBases(lastBefore - 3, 4)));
        assertSoftClipped(reads.get(0));
        assertSoftClipped(reads.get(2));
    }

    // a read with an unexplained clip waits for reads that could support it, and reads leave in order
    @Test
    public void testReadsLeaveInOrder() {
        final UnexplainedSoftClipTransformer transformer = newTransformer();
        final GATKRead first = makeRead(801, new byte[0], "60M", new byte[0]);
        final GATKRead waiting = makeRead(805, new byte[0], "60M10S", random.nextBases(10));
        final GATKRead blocked = makeRead(810, new byte[0], "60M", new byte[0]);
        final GATKRead otherContig = ArtificialReadUtils.createArtificialRead(header, "other", 1, 5, random.nextBases(ALIGNED_LENGTH),
                Utils.dupBytes(GOOD_QUALITY, ALIGNED_LENGTH), "60M");

        transformer.submit(first);
        Assert.assertEquals(transformer.consumeFinalizedItems(), List.of(first));
        transformer.submit(waiting);
        transformer.submit(blocked);
        Assert.assertFalse(transformer.hasFinalizedItems());
        // nothing on another contig can support a clip, so moving to one releases everything
        transformer.submit(otherContig);
        final List<GATKRead> released = transformer.consumeFinalizedItems();
        Assert.assertEquals(released.size(), 3);
        Assert.assertEquals(released.get(0).getCigar().toString(), "60M10H");
        Assert.assertSame(released.get(1), blocked);
        Assert.assertSame(released.get(2), otherContig);
        transformer.signalEndOfInput();
        Assert.assertFalse(transformer.hasFinalizedItems());
    }

    private UnexplainedSoftClipTransformer newTransformer() {
        return new UnexplainedSoftClipTransformer(referenceSource, MIN_TAIL_QUALITY);
    }

    // push coordinate-sorted reads through a new transformer
    private List<GATKRead> transform(final GATKRead... reads) {
        final UnexplainedSoftClipTransformer transformer = newTransformer();
        final List<GATKRead> result = new ArrayList<>();
        for (final GATKRead read : reads) {
            transformer.submit(read);
            result.addAll(transformer.consumeFinalizedItems());
        }
        transformer.signalEndOfInput();
        result.addAll(transformer.consumeFinalizedItems());
        Assert.assertEquals(result.size(), reads.length);
        return result;
    }

    private GATKRead transformOne(final GATKRead read) {
        return transform(read).get(0);
    }

    private static void assertSoftClipped(final GATKRead read) {
        Assert.assertTrue(read.getCigar().toString().contains("S") && !read.getCigar().toString().contains("H"), read.getName() + " " + read.getCigar());
    }

    private static void assertHardClipped(final GATKRead read) {
        Assert.assertTrue(read.getCigar().toString().contains("H") && !read.getCigar().toString().contains("S"), read.getName() + " " + read.getCigar());
    }

    private byte[] insertedBases(final int from, final int to) {
        return Arrays.copyOfRange(inserted, from, to);
    }

    // ALIGNED_LENGTH aligned bases ending at lastAligned, then a soft clip
    private GATKRead rightClipped(final String name, final int lastAligned, final byte[] clip) {
        return makeNamedRead(name, lastAligned - ALIGNED_LENGTH + 1, Utils.concat(referenceBases(lastAligned - ALIGNED_LENGTH + 1, ALIGNED_LENGTH), clip),
                ALIGNED_LENGTH + "M" + clip.length + "S");
    }

    // a soft clip, then ALIGNED_LENGTH aligned bases starting at firstAligned
    private GATKRead leftClipped(final String name, final int firstAligned, final byte[] clip) {
        return makeNamedRead(name, firstAligned, Utils.concat(clip, referenceBases(firstAligned, ALIGNED_LENGTH)), clip.length + "S" + ALIGNED_LENGTH + "M");
    }

    // 40 aligned bases on each side of the given bases inserted after INSERTION_POINT
    private GATKRead withInsertion(final String name, final byte[] insertion) {
        return makeNamedRead(name, INSERTION_POINT - 39, Utils.concat(referenceBases(INSERTION_POINT - 39, 40), insertion, referenceBases(INSERTION_POINT + 1, 40)),
                "40M" + insertion.length + "I40M");
    }

    // 40 aligned bases on each side of the deletion of DELETION_LENGTH bases at DELETION_START
    private GATKRead withDeletion(final String name) {
        return makeNamedRead(name, DELETION_START - 40, Utils.concat(referenceBases(DELETION_START - 40, 40), referenceBases(DELETION_START + DELETION_LENGTH, 40)),
                "40M" + DELETION_LENGTH + "D40M");
    }

    private GATKRead makeNamedRead(final String name, final int start, final byte[] bases, final String cigar) {
        return ArtificialReadUtils.createArtificialRead(header, name, 0, start, bases, Utils.dupBytes(GOOD_QUALITY, bases.length), cigar);
    }

    // a read with ALIGNED_LENGTH reference bases from start and the given soft clips
    private GATKRead makeRead(final int start, final byte[] leftClip, final String cigar, final byte[] rightClip) {
        final byte[] bases = Utils.concat(leftClip, referenceBases(start, ALIGNED_LENGTH), rightClip);
        final byte[] quals = new byte[bases.length];
        Arrays.fill(quals, GOOD_QUALITY);
        return ArtificialReadUtils.createArtificialRead(header, "read", 0, start, bases, quals, cigar);
    }

    private byte[] referenceBases(final int start, final int length) {
        return Arrays.copyOfRange(reference, start - 1, start - 1 + length);
    }

    // reference offsets [from, to) from the junction, in the direction away from the aligned bases
    private byte[] referenceClip(final int junction, final boolean right, final int from, final int to) {
        final byte[] result = new byte[to - from];
        for (int t = from; t < to; t++) {
            result[t - from] = reference[junction + (right ? t : -t) - 1];
        }
        return result;
    }

    // random bases, each G or T
    private byte[] randomGT(final int length) {
        final byte[] result = random.nextBases(length);
        for (int i = 0; i < length; i++) {
            result[i] = (byte) (result[i] == 'A' || result[i] == 'G' ? 'G' : 'T');
        }
        return result;
    }

    private static byte[] withMismatches(final byte[] bases, final int... indices) {
        final byte[] result = Arrays.copyOf(bases, bases.length);
        for (final int i : indices) {
            result[i] = (byte) (result[i] == 'A' ? 'C' : 'A');
        }
        return result;
    }

    private static void put(final byte[] reference, final int position, final byte[] bases) {
        System.arraycopy(bases, 0, reference, position - 1, bases.length);
    }
}
