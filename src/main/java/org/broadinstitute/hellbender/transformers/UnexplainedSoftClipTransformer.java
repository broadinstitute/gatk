package org.broadinstitute.hellbender.transformers;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Tail;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.iterators.PushPullTransformer;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * Hard clips soft-clipped bases that neither the nearby reference nor other reads explain.
 *
 * Mutect2 restores soft-clipped bases for assembly, and its active region detection counts a soft clip as indel
 * evidence.  That helps when the aligner clips an indel or a cluster of mismatches near the end of a read rather than
 * aligning through it, but the clipped bases of chimeric reads, adapter read-through and reads from paralogs come from
 * elsewhere, so they trigger spurious active regions and, threading through the assembly graph, many spurious
 * haplotypes.
 *
 * Each soft clip is judged without its tail of bases with quality at most the caller's minimum tail quality, which
 * assembly removes anyway.  Stage 1 judges the clip alone and keeps it if what remains:
 * <ul>
 *     <li>is not empty, and has mean base quality at least {@link #MIN_MEAN_BASE_QUALITY}</li>
 *     <li>isn't mostly one base (eg a poly-A tail or two-color poly-G), if it has at least
 *     {@link #MIN_LENGTH_TO_CHECK_COMPOSITION} bases</li>
 *     <li>is explained by the nearby reference, up to one mismatch per {@link #BASES_PER_MISMATCH} bases compared:
 *     either it continues the reference past the clip point (mismatches near the end of the read) or, if it has at
 *     least {@link #MIN_LENGTH_TO_PLACE_INDEL} bases, it matches the reference after a deletion or duplication of
 *     up to {@link #MAX_INDEL_LENGTH} bases, or after an insertion of novel bases.  Where the clip point is in a short
 *     tandem repeat of the reference, these explanations must place at least {@link #MIN_BASES_BEYOND_STR} clipped
 *     bases past the repeat, since bases inside a repeat match it at many offsets.</li>
 *     <li>or, where the clip point is in a short tandem repeat, matches the reference exactly from a position in the
 *     repeat (or just past it) to at least one base past it, within {@link #MAX_INDEL_LENGTH} bases of the clip point:
 *     a change in the length of the repeat that the aligner clipped rather than aligned, however short the clip</li>
 * </ul>
 * Stage 2 keeps a clip that passes stage 1's quality checks but that the reference can't explain, such as one of
 * novel inserted sequence or one too short to place, if other reads place the same sequence at its clip point:
 * <ul>
 *     <li>a soft clip on the other side, with its clip point within {@link #MAX_CLIP_POINT_SHIFT} bases, overlaps it
 *     so that both fit one inserted sequence, sharing at least {@link #MIN_SHARED_CLIP_BASES} clip bases, up to one
 *     mismatch per {@link #BASES_PER_MISMATCH} bases of overlap</li>
 *     <li>or another read aligns an insertion or deletion within {@link #MAX_CLIP_POINT_SHIFT} bases of its clip point,
 *     and the clip has at least {@link #MIN_LENGTH_TO_MATCH_ALIGNED_INDEL} bases and matches the haplotype with that
 *     indel where the clipped read's alignment places it, exactly if it has fewer than {@link #BASES_PER_MISMATCH}
 *     bases and otherwise up to one mismatch per {@link #BASES_PER_MISMATCH} bases</li>
 * </ul>
 * Several soft clips that agree on the same side don't count, since chimeric junk and reads from paralogs do that too.
 * A read's mate doesn't count either.
 *
 * Since stage 2 needs other reads, the transformer buffers the coordinate-sorted stream of reads and releases them in
 * order, each once no read still to come could support its clips.  It should see the reads of all samples together,
 * and, running after downsampling, it holds roughly one read length of reads.  Running before active region
 * detection, it saves the cost of spurious regions as well as of spurious haplotypes.
 */
public final class UnexplainedSoftClipTransformer implements PushPullTransformer<GATKRead> {
    public static final double MIN_MEAN_BASE_QUALITY = 20;
    public static final int MIN_LENGTH_TO_CHECK_COMPOSITION = 5;
    public static final double MAX_FRACTION_OF_ONE_BASE = 0.8;
    public static final int BASES_PER_MISMATCH = 10;
    public static final int MIN_LENGTH_TO_PLACE_INDEL = 8;
    public static final int MAX_INDEL_LENGTH = 50;
    public static final int MAX_STR_PERIOD = 6;
    public static final int MIN_STR_LENGTH = 8;
    public static final int MIN_STR_UNITS = 3;
    public static final int MIN_BASES_BEYOND_STR = 5;
    public static final int MAX_CLIP_POINT_SHIFT = 5;
    public static final int MIN_SHARED_CLIP_BASES = 10;
    public static final int MIN_LENGTH_TO_MATCH_ALIGNED_INDEL = 3;

    // reference fetched past the farthest base a stage-1 explanation can reach, so that a repeat there is seen whole
    private static final int REFERENCE_PADDING = 60;
    // reference bases next to a clip point that must agree with an aligned indel's haplotype: enough to reach past
    // the indel from any clip point within MAX_CLIP_POINT_SHIFT
    private static final int ANCHOR_LENGTH = MAX_CLIP_POINT_SHIFT + 1;
    // a read spanning more reference than this (eg an RNA read with a long N) has its right clip judged without waiting
    // for every read that might support it, which would hold the whole span in the buffer
    private static final int MAX_SPAN_TO_WAIT_FOR = 1000;

    private final ReferenceDataSource referenceDataSource;
    private final int minTailQuality;

    private final Deque<PendingRead> pending = new ArrayDeque<>();
    private List<GATKRead> finalized = new ArrayList<>();
    private String currentContig = null;
    private int currentStart = 0;   // the start of the newest read on the current contig

    // On the current contig, by clip point: the soft clips that pass stage 1's quality checks, and the indels that
    // reads align, keyed where a clip that they explain would be.  A clip point is the gap after a reference position:
    // the last aligned base, for a right clip, or the position before the first aligned base, for a left clip.
    private final TreeMap<Integer, List<Clip>> rightClips = new TreeMap<>();
    private final TreeMap<Integer, List<Clip>> leftClips = new TreeMap<>();
    private final TreeMap<Integer, Map<AlignedIndel, List<String>>> indelsForRightClips = new TreeMap<>();
    private final TreeMap<Integer, Map<AlignedIndel, List<String>>> indelsForLeftClips = new TreeMap<>();

    /**
     * @param minTailQuality the caller later hard clips read tails of bases with at most this quality, so the transformer
     *                       ignores them in judging a soft clip
     */
    public UnexplainedSoftClipTransformer(final ReferenceDataSource referenceDataSource, final int minTailQuality) {
        this.referenceDataSource = Utils.nonNull(referenceDataSource);
        this.minTailQuality = minTailQuality;
    }

    @Override
    public void submit(final GATKRead read) {
        Utils.nonNull(read);
        if (read.isUnmapped()) {
            pending.add(new PendingRead(read, ClipJudgment.NONE, ClipJudgment.NONE, Integer.MIN_VALUE, currentStart));
        } else {
            if (!read.getContig().equals(currentContig)) {
                release(true);
                clearIndex();
                currentContig = read.getContig();
            }
            currentStart = read.getStart();
            pending.add(index(read));
        }
        release(false);
    }

    @Override
    public boolean hasFinalizedItems() {
        return !finalized.isEmpty();
    }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        final List<GATKRead> result = finalized;
        finalized = new ArrayList<>();
        return result;
    }

    @Override
    public void signalEndOfInput() {
        release(true);
        clearIndex();
        currentContig = null;
    }

    // release the decidable reads at the front of the queue, or all of them, in order
    private void release(final boolean all) {
        while (!pending.isEmpty() && (all || pending.peekFirst().decideAfter() < currentStart)) {
            finalized.add(decide(pending.pollFirst()));
        }
        // forget clips and indels that no read still to be decided, or still to come, could use
        final int oldestNeeded = (pending.isEmpty() ? currentStart : pending.peekFirst().start()) - 1 - MAX_CLIP_POINT_SHIFT;
        rightClips.headMap(oldestNeeded).clear();
        leftClips.headMap(oldestNeeded).clear();
        indelsForRightClips.headMap(oldestNeeded).clear();
        indelsForLeftClips.headMap(oldestNeeded).clear();
    }

    private void clearIndex() {
        rightClips.clear();
        leftClips.clear();
        indelsForRightClips.clear();
        indelsForLeftClips.clear();
    }

    // judge a read's soft clips alone (stage 1), record them and its aligned indels, and work out when stage 2 can decide
    private PendingRead index(final GATKRead read) {
        recordAlignedIndels(read);
        // a read with no reference span has no clip point (and countClippedBases rejects a fully clipped cigar)
        if (read.getCigar().getReferenceLength() == 0) {
            return new PendingRead(read, ClipJudgment.NONE, ClipJudgment.NONE, Integer.MIN_VALUE, read.getStart());
        }
        final ClipJudgment left = judge(read, false);
        final ClipJudgment right = judge(read, true);
        if (left.passesQuality()) {
            leftClips.computeIfAbsent(read.getStart() - 1, k -> new ArrayList<>()).add(new Clip(read.getName(), left.forwardBases()));
        }
        if (right.passesQuality()) {
            rightClips.computeIfAbsent(read.getEnd(), k -> new ArrayList<>()).add(new Clip(read.getName(), right.forwardBases()));
        }
        // Support for a right clip at point p comes from reads that start by p + 1 + MAX_CLIP_POINT_SHIFT: left clips
        // at nearby points, and reads that align nearby indels.  Support for a left clip comes from reads that start by
        // its point + MAX_CLIP_POINT_SHIFT.
        final int decideAfter = right.verdict() == Verdict.UNEXPLAINED ?
                Math.min(read.getEnd(), read.getStart() + MAX_SPAN_TO_WAIT_FOR) + 1 + MAX_CLIP_POINT_SHIFT :
                (left.verdict() == Verdict.UNEXPLAINED ? read.getStart() - 1 + MAX_CLIP_POINT_SHIFT : Integer.MIN_VALUE);
        return new PendingRead(read, left, right, decideAfter, read.getStart());
    }

    private void recordAlignedIndels(final GATKRead read) {
        int position = read.getStart();
        int readOffset = 0;
        for (final CigarElement element : read.getCigarElements()) {
            final CigarOperator operator = element.getOperator();
            final int length = element.getLength();
            if (operator == CigarOperator.I) {
                final AlignedIndel insertion = new AlignedIndel(position - 1, position - 1,
                        new String(read.getBasesNoCopy(), readOffset, length, StandardCharsets.US_ASCII));
                recordIndel(indelsForRightClips, position - 1, insertion, read.getName());
                recordIndel(indelsForLeftClips, position - 1, insertion, read.getName());
            } else if (operator == CigarOperator.D) {
                final AlignedIndel deletion = new AlignedIndel(position, position + length - 1, "");
                recordIndel(indelsForRightClips, position - 1, deletion, read.getName());
                recordIndel(indelsForLeftClips, position + length - 1, deletion, read.getName());
            }
            position += operator.consumesReferenceBases() ? length : 0;
            readOffset += operator.consumesReadBases() ? length : 0;
        }
    }

    // keep up to two read names per indel, enough to tell whether a read other than a clip's own aligns it
    private static void recordIndel(final TreeMap<Integer, Map<AlignedIndel, List<String>>> indels, final int point,
                                    final AlignedIndel indel, final String readName) {
        final List<String> names = indels.computeIfAbsent(point, k -> new HashMap<>()).computeIfAbsent(indel, k -> new ArrayList<>(2));
        if (names.size() < 2 && !names.contains(readName)) {
            names.add(readName);
        }
    }

    private GATKRead decide(final PendingRead pendingRead) {
        final GATKRead read = pendingRead.read();
        final boolean keepLeft = keep(read, false, pendingRead.left());
        final boolean keepRight = keep(read, true, pendingRead.right());
        if (keepLeft && keepRight) {
            return read;
        }
        // the clipper applies its ops in turn to the shrinking read, so clip the right end first
        final ReadClipper clipper = new ReadClipper(read);
        if (!keepRight) {
            final int rightClip = CigarUtils.countClippedBases(read.getCigar(), Tail.RIGHT, CigarOperator.SOFT_CLIP);
            clipper.addOp(new ClippingOp(read.getLength() - rightClip, read.getLength() - 1));
        }
        if (!keepLeft) {
            final int leftClip = CigarUtils.countClippedBases(read.getCigar(), Tail.LEFT, CigarOperator.SOFT_CLIP);
            clipper.addOp(new ClippingOp(0, leftClip - 1));
        }
        return clipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);
    }

    private boolean keep(final GATKRead read, final boolean right, final ClipJudgment judgment) {
        return switch (judgment.verdict()) {
            case NO_CLIP, EXPLAINED -> true;
            case POOR_QUALITY -> false;
            case UNEXPLAINED -> isPlacedByOtherReads(read.getName(), right, right ? read.getEnd() : read.getStart() - 1, judgment.forwardBases());
        };
    }

    // ------------------------------------------------------------------------------------------------ stage 1

    private ClipJudgment judge(final GATKRead read, final boolean right) {
        final int length = CigarUtils.countClippedBases(read.getCigar(), right ? Tail.RIGHT : Tail.LEFT, CigarOperator.SOFT_CLIP);
        if (length == 0) {
            return ClipJudgment.NONE;
        }
        final byte[] bases = read.getBasesNoCopy();
        final byte[] quals = read.getBaseQualitiesNoCopy();

        // judge the clip without its low-quality tail, at the outer end of the read
        int judged = length;
        while (judged > 0 && quals[right ? read.getLength() - length + judged - 1 : length - judged] <= minTailQuality) {
            judged--;
        }
        if (judged == 0) {
            return ClipJudgment.POOR;
        }
        final int readStart = right ? read.getLength() - length : length - judged;

        long qualSum = 0;
        final int[] baseCounts = new int[256];
        for (int i = readStart; i < readStart + judged; i++) {
            qualSum += quals[i];
            baseCounts[bases[i] & 0xFF]++;
        }
        if (qualSum < MIN_MEAN_BASE_QUALITY * judged) {
            return ClipJudgment.POOR;
        } else if (judged >= MIN_LENGTH_TO_CHECK_COMPOSITION && MathUtils.arrayMax(baseCounts) >= MAX_FRACTION_OF_ONE_BASE * judged) {
            return ClipJudgment.POOR;
        }
        final byte[] forwardBases = Arrays.copyOfRange(bases, readStart, readStart + judged);
        final boolean explained = new OrientedClip(read, right, forwardBases).isExplained();
        return new ClipJudgment(explained ? Verdict.EXPLAINED : Verdict.UNEXPLAINED, forwardBases);
    }

    /**
     * A soft clip and the reference around it, oriented away from the clip point: base j of the clip is the j-th base
     * past the clip point, and offset t is the t-th reference base past it, with negative offsets back over the
     * aligned part of the read.
     */
    private final class OrientedClip {
        private final byte[] clipBases;
        private final int direction;
        private final int firstPosition;    // the reference position of offset 0
        private final int windowStart;
        private final byte[] window;
        private final int repeatStart;      // the reference repeat at the clip point, as half-open positions
        private final int repeatEnd;        // repeatStart == repeatEnd if there is none

        // like the reference source, which upper-cases its bases, we assume upper-case read bases
        OrientedClip(final GATKRead read, final boolean right, final byte[] forwardBases) {
            clipBases = forwardBases.clone();
            if (!right) {
                ArrayUtils.reverse(clipBases);
            }
            direction = right ? 1 : -1;
            firstPosition = right ? read.getEnd() + 1 : read.getStart() - 1;

            final int reach = MAX_INDEL_LENGTH + clipBases.length + REFERENCE_PADDING;
            windowStart = Math.max(1, firstPosition - reach);
            window = referenceBases(windowStart, firstPosition + reach);

            // the clip point lies between reference positions boundary - 1 and boundary
            final int boundary = right ? read.getEnd() + 1 : read.getStart();
            final int[] repeat = findRepeat(boundary);
            repeatStart = repeat[0];
            repeatEnd = repeat[1];
        }

        boolean isExplained() {
            final int n = clipBases.length;
            if (matches(0, 0, n) && (n < MIN_LENGTH_TO_PLACE_INDEL || reachesPastRepeat(0, n))) {
                return true;    // continues the reference
            } else if (changesRepeatLength()) {
                return true;
            } else if (n < MIN_LENGTH_TO_PLACE_INDEL) {
                return false;
            }
            for (int offset = 1; offset <= MAX_INDEL_LENGTH; offset++) {
                if ((matches(0, offset, n) && reachesPastRepeat(offset, n)) ||        // a deletion
                        (matches(0, -offset, n) && reachesPastRepeat(-offset, n))) {  // a duplication
                    return true;
                }
            }
            for (int inserted = 1; inserted <= n - MIN_LENGTH_TO_PLACE_INDEL; inserted++) {
                if (matches(inserted, 0, n - inserted) && reachesPastRepeat(0, n - inserted)) {
                    return true;    // novel inserted bases, then the reference
                }
            }
            return false;
        }

        // whether the clip matches the reference exactly from a position in the repeat at the clip point, or the first
        // position past its far end, to at least one base past it: units of the repeat deleted or added, then its flank
        private boolean changesRepeatLength() {
            if (repeatStart == repeatEnd) {
                return false;
            }
            final int nearOffset = offset(direction > 0 ? repeatStart : repeatEnd - 1);
            final int beyondOffset = offset(direction > 0 ? repeatEnd : repeatStart - 1);
            final int lastOffset = Math.min(beyondOffset, MAX_INDEL_LENGTH);
            for (int offset = Math.max(nearOffset, Math.max(beyondOffset - clipBases.length + 1, -MAX_INDEL_LENGTH)); offset <= lastOffset; offset++) {
                if (matches(0, offset, clipBases.length, 0)) {
                    return true;
                }
            }
            return false;
        }

        private int position(final int offset) {
            return firstPosition + direction * offset;
        }

        private int offset(final int position) {
            return direction * (position - firstPosition);
        }

        // whether clip bases [clipStart, clipStart + count) match the reference from offset on, up to the mismatch limit
        private boolean matches(final int clipStart, final int offset, final int count) {
            return matches(clipStart, offset, count, Math.max(1, count / BASES_PER_MISMATCH));
        }

        private boolean matches(final int clipStart, final int offset, final int count, final int maxMismatches) {
            int mismatches = 0;
            for (int j = 0; j < count; j++) {
                final int index = position(offset + j) - windowStart;
                if (index < 0 || index >= window.length) {
                    return false;
                } else if (window[index] != clipBases[clipStart + j] && ++mismatches > maxMismatches) {
                    return false;
                }
            }
            return true;
        }

        // whether at least MIN_BASES_BEYOND_STR of count bases placed from offset on lie outside the repeat
        private boolean reachesPastRepeat(final int offset, final int count) {
            if (repeatStart == repeatEnd) {
                return true;
            }
            int beyond = 0;
            for (int j = 0; j < count; j++) {
                final int position = position(offset + j);
                if (position < repeatStart || position >= repeatEnd) {
                    beyond++;
                }
            }
            return beyond >= MIN_BASES_BEYOND_STR;
        }

        // the longest tandem repeat of the reference touching the boundary, as half-open positions {start, end}, or {0, 0}
        private int[] findRepeat(final int boundary) {
            final int boundaryIndex = boundary - windowStart;
            int bestStart = 0;  // indices into the window
            int bestEnd = 0;
            for (int period = 1; period <= MAX_STR_PERIOD; period++) {
                for (int anchor = boundaryIndex - 1; anchor <= boundaryIndex; anchor++) {
                    // the repeat unit may begin up to period - 1 bases before the anchor
                    for (int unitStart = anchor - period + 1; unitStart <= anchor; unitStart++) {
                        if (!isRegularUnit(unitStart, period)) {
                            continue;
                        }
                        // whole units before and after, then partial units, which can't run into an N since the unit has none
                        int start = unitStart - period * GATKVariantContextUtils.findNumberOfRepetitions(
                                window, unitStart, period, window, 0, unitStart, false);
                        int end = unitStart + period + period * GATKVariantContextUtils.findNumberOfRepetitions(
                                window, unitStart, period, window, unitStart + period, window.length - unitStart - period, true);
                        while (start > 0 && window[start - 1] == window[start - 1 + period]) {
                            start--;
                        }
                        while (end < window.length && window[end] == window[end - period]) {
                            end++;
                        }
                        if (end - start >= Math.max(MIN_STR_LENGTH, MIN_STR_UNITS * period) && start <= boundaryIndex &&
                                boundaryIndex <= end && end - start > bestEnd - bestStart) {
                            bestStart = start;
                            bestEnd = end;
                        }
                    }
                }
            }
            return bestStart == bestEnd ? new int[] {0, 0} : new int[] {bestStart + windowStart, bestEnd + windowStart};
        }

        // whether the window holds the period bases from unitStart, none of them N
        private boolean isRegularUnit(final int unitStart, final int period) {
            if (unitStart < 0 || unitStart + period > window.length) {
                return false;
            }
            for (int i = unitStart; i < unitStart + period; i++) {
                if (!BaseUtils.isRegularBase(window[i])) {
                    return false;
                }
            }
            return true;
        }
    }

    // ------------------------------------------------------------------------------------------------ stage 2

    private boolean isPlacedByOtherReads(final String readName, final boolean right, final int point, final byte[] forwardBases) {
        final NavigableMap<Integer, List<Clip>> otherSide = (right ? leftClips : rightClips)
                .subMap(point - MAX_CLIP_POINT_SHIFT, true, point + MAX_CLIP_POINT_SHIFT, true);
        for (final Map.Entry<Integer, List<Clip>> entry : otherSide.entrySet()) {
            for (final Clip other : entry.getValue()) {
                if (!other.readName().equals(readName) && (right ? fitOneInsertion(point, forwardBases, entry.getKey(), other.forwardBases()) :
                        fitOneInsertion(entry.getKey(), other.forwardBases(), point, forwardBases))) {
                    return true;
                }
            }
        }
        if (forwardBases.length < MIN_LENGTH_TO_MATCH_ALIGNED_INDEL) {
            return false;
        }
        final NavigableMap<Integer, Map<AlignedIndel, List<String>>> indels = (right ? indelsForRightClips : indelsForLeftClips)
                .subMap(point - MAX_CLIP_POINT_SHIFT, true, point + MAX_CLIP_POINT_SHIFT, true);
        for (final Map<AlignedIndel, List<String>> indelsAtPoint : indels.values()) {
            for (final Map.Entry<AlignedIndel, List<String>> entry : indelsAtPoint.entrySet()) {
                if (entry.getValue().stream().anyMatch(name -> !name.equals(readName)) &&
                        matchesAlignedIndel(right, point, forwardBases, entry.getKey())) {
                    return true;
                }
            }
        }
        return false;
    }

    // Whether a right clip and a left clip fit one inserted sequence: with A the reference up to the right clip's point,
    // then its bases, and B the left clip's bases, then the reference after its point, some suffix of A matches a
    // prefix of B in which at least MIN_SHARED_CLIP_BASES bases of the two clips face each other.  Flanking reference
    // bases in the overlap must match too, which allows the two clip points to differ, but they don't count as shared.
    private boolean fitOneInsertion(final int rightPoint, final byte[] rightBases, final int leftPoint, final byte[] leftBases) {
        final byte[] a = Utils.concat(referenceBases(rightPoint - leftBases.length - 2 * MAX_CLIP_POINT_SHIFT + 1, rightPoint), rightBases);
        final byte[] b = Utils.concat(leftBases, referenceBases(leftPoint + 1, leftPoint + rightBases.length + 2 * MAX_CLIP_POINT_SHIFT));
        for (int overlap = MIN_SHARED_CLIP_BASES; overlap <= Math.min(a.length, b.length); overlap++) {
            final int shared = Math.min(rightBases.length, overlap) + Math.min(leftBases.length, overlap) - overlap;
            if (shared >= MIN_SHARED_CLIP_BASES &&
                    !tooManyMismatches(a, a.length - overlap, b, 0, overlap, Math.max(1, overlap / BASES_PER_MISMATCH))) {
                return true;
            }
        }
        return false;
    }

    // Whether a clip matches the haplotype with an aligned indel, placed where the clipped read's alignment puts it:
    // the reference bases on the aligned side of the clip point, which must match the haplotype exactly, anchor it.
    private boolean matchesAlignedIndel(final boolean right, final int point, final byte[] forwardBases, final AlignedIndel indel) {
        final byte[] inserted = indel.insertedBases().getBytes(StandardCharsets.US_ASCII);
        final int lastBefore = indel.isInsertion() ? indel.start() : indel.start() - 1;   // the haplotype is the reference
        final int firstAfter = indel.end() + 1;                                           // through lastBefore, the inserted
        final int shift = inserted.length - (firstAfter - lastBefore - 1);                // bases, then the reference from firstAfter
        final int anchorStart = right ? point - ANCHOR_LENGTH + 1 : point + 1;
        final int anchorEnd = right ? point : point + ANCHOR_LENGTH;
        if (anchorStart < 1 || anchorEnd > contigLength()) {
            return false;
        }
        final int windowStart = Math.max(1, Math.min(point, lastBefore) - ANCHOR_LENGTH - forwardBases.length - 1);
        final int windowEnd = Math.max(point, firstAfter) + ANCHOR_LENGTH + forwardBases.length + 1;
        final byte[] haplotype = Utils.concat(referenceBases(windowStart, lastBefore), inserted, referenceBases(firstAfter, windowEnd));

        final byte[] anchor = referenceBases(anchorStart, anchorEnd);
        final byte[] placed = right ? Utils.concat(anchor, forwardBases) : Utils.concat(forwardBases, anchor);
        // the haplotype index of the placed bases: a right clip's anchor begins left of the indel, a left clip's ends right of it
        final int start = right ? anchorStart - windowStart : anchorEnd - windowStart + shift - placed.length + 1;
        if (start < 0 || start + placed.length > haplotype.length) {
            return false;
        }
        final int anchorOffset = right ? 0 : forwardBases.length;
        final int clipOffset = right ? ANCHOR_LENGTH : 0;
        final int maxMismatches = forwardBases.length < BASES_PER_MISMATCH ? 0 : forwardBases.length / BASES_PER_MISMATCH;
        return !tooManyMismatches(placed, anchorOffset, haplotype, start + anchorOffset, ANCHOR_LENGTH, 0) &&
                !tooManyMismatches(placed, clipOffset, haplotype, start + clipOffset, forwardBases.length, maxMismatches);
    }

    // whether more than maxMismatches of count bases from xStart and yStart differ, stopping early
    private static boolean tooManyMismatches(final byte[] x, final int xStart, final byte[] y, final int yStart, final int count,
                                             final int maxMismatches) {
        int mismatches = 0;
        for (int i = 0; i < count; i++) {
            if (x[xStart + i] != y[yStart + i] && ++mismatches > maxMismatches) {
                return true;
            }
        }
        return false;
    }

    private int contigLength() {
        return referenceDataSource.getSequenceDictionary().getSequence(currentContig).getSequenceLength();
    }

    // the reference bases of the current contig from start to end, inclusive, clipped to the contig
    private byte[] referenceBases(final int start, final int end) {
        final int from = Math.max(1, start);
        final int to = Math.min(contigLength(), end);
        return from > to ? new byte[0] : referenceDataSource.queryAndPrefetch(currentContig, from, to).getBases();
    }

    private enum Verdict { NO_CLIP, POOR_QUALITY, EXPLAINED, UNEXPLAINED }

    // stage 1's verdict on a read end, and the bases it judged, in reference order
    private record ClipJudgment(Verdict verdict, byte[] forwardBases) {
        static final ClipJudgment NONE = new ClipJudgment(Verdict.NO_CLIP, null);
        static final ClipJudgment POOR = new ClipJudgment(Verdict.POOR_QUALITY, null);

        boolean passesQuality() {
            return verdict == Verdict.EXPLAINED || verdict == Verdict.UNEXPLAINED;
        }
    }

    // decidable once the newest read starts after decideAfter
    private record PendingRead(GATKRead read, ClipJudgment left, ClipJudgment right, int decideAfter, int start) { }

    private record Clip(String readName, byte[] forwardBases) { }

    // an insertion after reference position start == end, or a deletion of reference positions start through end
    private record AlignedIndel(int start, int end, String insertedBases) {
        boolean isInsertion() {
            return !insertedBases.isEmpty();
        }
    }
}
