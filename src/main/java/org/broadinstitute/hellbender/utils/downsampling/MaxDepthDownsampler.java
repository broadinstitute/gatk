package org.broadinstitute.hellbender.utils.downsampling;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Max-depth downsampler: caps the per-sample depth at every reference position at {@code maxDepth}, while
 * guaranteeing that no position ends up with fewer reads than {@code min(original depth, maxDepth)}.
 *
 * Reads are buffered in windows of {@code windowSize} bases of alignment start. When a window is finalized,
 * its reads are visited in random order (or arrival order in non-random mode) and a read is discarded only
 * if every position it spans already holds {@code maxDepth} kept reads; otherwise it is kept and its span is
 * counted. Kept reads are emitted in their original coordinate order.
 *
 * Because a read can only be discarded once {@code maxDepth} reads already cover its whole span, a window in
 * which no position can be covered by more than {@code maxDepth} reads (bounded by the still-open kept reads
 * plus the most window reads starting within one read span of each other) cannot discard anything; such
 * windows are passed through without any per-base work. In ordinary sequence that is almost every window, so
 * the downsampler only does real work in collapsed repeats and other ultra-deep regions. Any region whose
 * depth never exceeds {@code maxDepth} is passed through unchanged.
 *
 * Depth is counted per sample over the aligned span of each read (soft clips excluded, deletions included).
 * Depth is only tracked for one window's width past the end of the current window; a read whose span reaches
 * beyond that is always kept, which keeps memory and work proportional to the window rather than to the
 * longest read span. Reads without an assigned position are passed through untouched, after any buffered
 * positioned reads.
 *
 * The cap is a floor-preserving target rather than a hard ceiling: a read is judged before the later reads that
 * will cover its right end have arrived, so kept depth overshoots the cap by up to about one read's worth of
 * reads at the end of each window. Windows should therefore be several read lengths wide; the overshoot grows
 * as the window shrinks below a read length.
 *
 * The downsampler can be reused across inputs: {@link #signalEndOfInput()} flushes the current window and
 * resets all positional state, so a subsequent input may start at any position.
 */
public final class MaxDepthDownsampler extends ReadsDownsampler {

    private final int maxDepth;
    private final int windowSize;
    private final SAMFileHeader header;
    /** Source of the per-window visiting order; null means arrival order (deterministic, for tests). */
    private final Random random;

    private String windowContig;
    private int windowStart;
    private final List<GATKRead> pendingWindow;
    /** Kept reads, per sample, that may still overlap positions at or beyond the current window start. */
    private final Map<String, List<GATKRead>> openKeptBySample;
    private List<GATKRead> finalizedReads;
    private GATKRead previousRead;

    /**
     * @param maxDepth maximum kept depth per sample at any position; must be > 0
     * @param windowSize width, in bases of alignment start, of the windows within which reads are randomly
     *                   ordered before capping; must be > 0
     * @param header header used to resolve sample names and to check sort order
     * @param random source of randomness for the per-window visiting order, or null to visit reads in
     *               arrival order
     */
    public MaxDepthDownsampler(final int maxDepth, final int windowSize, final SAMFileHeader header, final Random random) {
        Utils.validateArg(maxDepth > 0, "maxDepth must be > 0");
        Utils.validateArg(windowSize > 0, "windowSize must be > 0");
        this.maxDepth = maxDepth;
        this.windowSize = windowSize;
        this.header = Utils.nonNull(header);
        this.random = random;
        this.pendingWindow = new ArrayList<>();
        this.openKeptBySample = new HashMap<>();
        this.finalizedReads = new ArrayList<>();
        clearItems();
        resetStats();
    }

    @Override
    public void submit(final GATKRead newRead) {
        Utils.nonNull(newRead, "newRead");
        if (ReadUtils.readHasNoAssignedPosition(newRead)) {
            // Unplaced reads sort after all positioned reads, so flush the window first to keep the output ordered.
            finalizeWindow();
            finalizedReads.add(newRead);
            return;
        }
        checkSortOrder(newRead);
        if (!isInCurrentWindow(newRead)) {
            finalizeWindow();
            if (windowContig == null || !windowContig.equals(newRead.getAssignedContig())) {
                openKeptBySample.clear();
            }
            windowContig = newRead.getAssignedContig();
            windowStart = newRead.getAssignedStart();
        }
        pendingWindow.add(newRead);
        previousRead = newRead;
    }

    private boolean isInCurrentWindow(final GATKRead read) {
        return windowContig != null && windowContig.equals(read.getAssignedContig())
                && read.getAssignedStart() < (long) windowStart + windowSize;
    }

    private void checkSortOrder(final GATKRead newRead) {
        if (previousRead != null && ReadCoordinateComparator.compareCoordinates(previousRead, newRead, header) > 0) {
            throw new IllegalStateException(
                    String.format("Reads must be coordinate sorted (earlier %s later %s)", previousRead, newRead));
        }
    }

    private void finalizeWindow() {
        if (pendingWindow.isEmpty()) {
            return;
        }
        // Group the window's reads by sample, remembering each read's position in the window so that the
        // survivors can be emitted in their original order.
        final Map<String, List<Integer>> indicesBySample = new HashMap<>();
        for (int i = 0; i < pendingWindow.size(); i++) {
            indicesBySample.computeIfAbsent(sampleKey(pendingWindow.get(i)), k -> new ArrayList<>()).add(i);
        }
        final boolean[] keep = new boolean[pendingWindow.size()];
        for (final Map.Entry<String, List<Integer>> entry : indicesBySample.entrySet()) {
            capSample(entry.getKey(), entry.getValue(), keep);
        }
        for (int i = 0; i < pendingWindow.size(); i++) {
            if (keep[i]) {
                finalizedReads.add(pendingWindow.get(i));
            } else {
                incrementNumberOfDiscardedItems(1);
            }
        }
        pendingWindow.clear();
    }

    /** Decides, for one sample, which of the window's reads (given by index into pendingWindow) survive. */
    private void capSample(final String sample, final List<Integer> indices, final boolean[] keep) {
        final List<GATKRead> open = openKeptBySample.computeIfAbsent(sample, k -> new ArrayList<>());
        // Kept reads ending before this window starts cannot overlap anything still to come.
        open.removeIf(r -> spanEnd(r) < windowStart);

        if (open.size() + maxReadsStartingWithinOneSpan(indices) <= maxDepth) {
            // No position can reach maxDepth before the last read covering it is placed, so nothing can be discarded.
            for (final int i : indices) {
                keep[i] = true;
                open.add(pendingWindow.get(i));
            }
            return;
        }

        // Depth is tracked from the window start to one window's width past its end; reads reaching beyond
        // that are always kept, so the array never grows with a single read's span.
        final int trackedEnd = (int) Math.min((long) windowStart + 2L * windowSize - 1, Integer.MAX_VALUE);
        final int[] depth = new int[trackedEnd - windowStart + 1];
        for (final GATKRead r : open) {
            addSpan(depth, Math.max(r.getAssignedStart(), windowStart), Math.min(spanEnd(r), trackedEnd));
        }

        final List<Integer> order = new ArrayList<>(indices);
        if (random != null) {
            Collections.shuffle(order, random);
        }
        for (final int i : order) {
            final GATKRead read = pendingWindow.get(i);
            final int start = read.getAssignedStart();
            final int end = spanEnd(read);
            if (end > trackedEnd || anyPositionBelowCap(depth, start, end)) {
                keep[i] = true;
                addSpan(depth, start, Math.min(end, trackedEnd));
                open.add(read);
            }
        }
    }

    /**
     * Upper bound on the number of the window's reads covering any single position: the most reads whose starts
     * fall within one maximal read span of each other. Reads are in start order, so this is a sliding count.
     */
    private int maxReadsStartingWithinOneSpan(final List<Integer> indices) {
        long maxSpan = 1;
        for (final int i : indices) {
            final GATKRead r = pendingWindow.get(i);
            maxSpan = Math.max(maxSpan, (long) spanEnd(r) - r.getAssignedStart() + 1);
        }
        int best = 0;
        int left = 0;
        for (int right = 0; right < indices.size(); right++) {
            final long rightStart = pendingWindow.get(indices.get(right)).getAssignedStart();
            while (pendingWindow.get(indices.get(left)).getAssignedStart() <= rightStart - maxSpan) {
                left++;
            }
            best = Math.max(best, right - left + 1);
        }
        return best;
    }

    /** Aligned end of the read's span, falling back to its assigned start for reads that report no span. */
    private static int spanEnd(final GATKRead read) {
        return Math.max(read.getEnd(), read.getAssignedStart());
    }

    private boolean anyPositionBelowCap(final int[] depth, final int start, final int end) {
        for (int p = start; p <= end; p++) {
            if (depth[p - windowStart] < maxDepth) {
                return true;
            }
        }
        return false;
    }

    private void addSpan(final int[] depth, final int start, final int end) {
        for (int p = start; p <= end; p++) {
            depth[p - windowStart]++;
        }
    }

    private String sampleKey(final GATKRead read) {
        final String sample = ReadUtils.getSampleName(read, header);
        return sample == null ? "" : sample;
    }

    @Override
    public boolean hasFinalizedItems() {
        return !finalizedReads.isEmpty();
    }

    @Override
    public List<GATKRead> consumeFinalizedItems() {
        final List<GATKRead> toReturn = finalizedReads;
        finalizedReads = new ArrayList<>();
        return toReturn;
    }

    @Override
    public boolean hasPendingItems() {
        return !pendingWindow.isEmpty();
    }

    @Override
    public GATKRead peekFinalized() {
        return finalizedReads.isEmpty() ? null : finalizedReads.get(0);
    }

    @Override
    public GATKRead peekPending() {
        return pendingWindow.isEmpty() ? null : pendingWindow.get(0);
    }

    @Override
    public int size() {
        return finalizedReads.size() + pendingWindow.size();
    }

    @Override
    public void signalEndOfInput() {
        finalizeWindow();
        // The next input, if any, may start anywhere, and kept reads from this input must not count against it.
        openKeptBySample.clear();
        windowContig = null;
        previousRead = null;
    }

    @Override
    public void clearItems() {
        pendingWindow.clear();
        openKeptBySample.clear();
        finalizedReads.clear();
        windowContig = null;
        previousRead = null;
    }

    @Override
    public boolean requiresCoordinateSortOrder() {
        return true;
    }

    @Override
    public void signalNoMoreReadsBefore(final GATKRead read) {
        Utils.nonNull(read, "read");
        if (windowContig != null && !isInCurrentWindow(read)) {
            finalizeWindow();
        }
    }
}
