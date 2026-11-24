package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import org.broadinstitute.hellbender.utils.BaseUtils;

import java.util.Arrays;

/**
 * A kmer graph whose nodes carry a rank, roughly a reference-window offset.  The reference is threaded first, one node
 * per window offset, so repeated reference kmers are distinct nodes.  A read kmer may only join a node ranked after
 * the node of the read's previous kmer (before, threading backwards), so every edge points to a higher rank and the
 * graph is acyclic for any kmer size; among several copies of a kmer the read's alignment and the bases that follow
 * choose one.  A read kmer that joins no node gets a new node ranked just after (before) its neighbor.
 */
final class PositionalGraph {
    static final double EPSILON = 1e-4;
    /** how many read bases {@link #agreement} compares at most */
    private static final int AGREEMENT_CAP = 150;
    /** how many read bases {@link #deletionBeatsAlternatives} compares at most against the reference */
    private static final int DELETION_AGREEMENT_CAP = 200;

    private final byte[] reference;
    private final int k;
    private final double delta;
    private final int maxJump;
    private final int lookahead;

    private int size = 0;
    private long[] kmers = new long[0];
    private double[] ranks = new double[0];
    private int[] referenceOffsets = new int[0];
    private byte[] firstBases = new byte[0];
    private byte[] lastBases = new byte[0];
    // insertion-ordered adjacency lists with edge weights (the number of read kmer pairs, 0 for reference-only edges)
    private IntArrayList[] successors = new IntArrayList[0];
    private IntArrayList[] successorWeights = new IntArrayList[0];
    private IntArrayList[] predecessors = new IntArrayList[0];
    private IntArrayList[] predecessorWeights = new IntArrayList[0];

    /** kmer -> nodes with that kmer, in order of creation */
    private final Long2ObjectOpenHashMap<IntArrayList> index = new Long2ObjectOpenHashMap<>();
    private final LongOpenHashSet repeatedInReference = new LongOpenHashSet();

    /** where a walk from a dangling end rejoins the reference: the nodes walked, ending with a reference node */
    record Walk(int[] nodes, int referenceNode) { }

    private static final Walk NO_WALK = new Walk(new int[0], -1);
    private final Long2ObjectOpenHashMap<Walk> walks = new Long2ObjectOpenHashMap<>();

    /**
     * @param reference the reference window, upper case
     * @param delta how far (in bases) a repeated kmer may join a node from where the read expects it
     * @param maxJump how far a kmer unique in the graph may join a node from where the read expects it ...
     * @param lookahead ... and, beyond delta, only if the read then matches the reference for this many more bases
     */
    PositionalGraph(final byte[] reference, final int k, final double delta, final int maxJump, final int lookahead) {
        this.reference = reference;
        this.k = k;
        this.delta = delta;
        this.maxJump = maxJump;
        this.lookahead = lookahead;
        final int n = Math.max(0, reference.length - k + 1);
        final Long2IntOpenHashMap counts = new Long2IntOpenHashMap();
        int previous = -1;
        for (int i = 0; i < n; i++) {
            final long kmer = encode(reference, i, k);
            final int node = newNode(kmer, i, i, reference[i], reference[i + k - 1], kmer >= 0);
            if (kmer >= 0) {
                counts.addTo(kmer, 1);
            }
            if (previous >= 0) {
                addEdge(previous, node, 0);
            }
            previous = node;
        }
        counts.long2IntEntrySet().forEach(e -> {
            if (e.getIntValue() > 1) {
                repeatedInReference.add(e.getLongKey());
            }
        });
    }

    /** 2-bit encoding of bases[start, start + k), or -1 if a base is not A, C, G or T */
    static long encode(final byte[] bases, final int start, final int k) {
        long code = 0;
        for (int i = start; i < start + k; i++) {
            final int base = BaseUtils.simpleBaseToBaseIndex(bases[i]);
            if (base < 0) {
                return -1;
            }
            code = (code << 2) | base;
        }
        return code;
    }

    int size() {
        return size;
    }

    int referenceOffset(final int node) {
        return referenceOffsets[node];
    }

    boolean isReference(final int node) {
        return referenceOffsets[node] >= 0;
    }

    byte firstBase(final int node) {
        return firstBases[node];
    }

    byte lastBase(final int node) {
        return lastBases[node];
    }

    private int newNode(final long kmer, final double rank, final int referenceOffset, final byte first, final byte last, final boolean indexed) {
        if (size == kmers.length) {
            final int capacity = Math.max(1024, 2 * size);
            kmers = Arrays.copyOf(kmers, capacity);
            ranks = Arrays.copyOf(ranks, capacity);
            referenceOffsets = Arrays.copyOf(referenceOffsets, capacity);
            firstBases = Arrays.copyOf(firstBases, capacity);
            lastBases = Arrays.copyOf(lastBases, capacity);
            successors = Arrays.copyOf(successors, capacity);
            successorWeights = Arrays.copyOf(successorWeights, capacity);
            predecessors = Arrays.copyOf(predecessors, capacity);
            predecessorWeights = Arrays.copyOf(predecessorWeights, capacity);
        }
        final int node = size++;
        kmers[node] = kmer;
        ranks[node] = rank;
        referenceOffsets[node] = referenceOffset;
        firstBases[node] = first;
        lastBases[node] = last;
        successors[node] = new IntArrayList(2);
        successorWeights[node] = new IntArrayList(2);
        predecessors[node] = new IntArrayList(2);
        predecessorWeights[node] = new IntArrayList(2);
        if (indexed) {
            index.computeIfAbsent(kmer, x -> new IntArrayList(1)).add(node);
        }
        return node;
    }

    private void addEdge(final int from, final int to, final int weight) {
        addWeight(successors[from], successorWeights[from], to, weight);
        addWeight(predecessors[to], predecessorWeights[to], from, weight);
    }

    private static void addWeight(final IntArrayList targets, final IntArrayList weights, final int target, final int weight) {
        final int i = targets.indexOf(target);
        if (i < 0) {
            targets.add(target);
            weights.add(weight);
        } else {
            weights.set(i, weights.getInt(i) + weight);
        }
    }

    private boolean isUnique(final long kmer) {
        final IntArrayList nodes = index.get(kmer);
        return nodes != null && nodes.size() == 1 && !repeatedInReference.contains(kmer);
    }

    /**
     * May a kmer unique in the graph join node c although it is dist (more than delta) from where the read expects
     * it?  Only within maxJump and, for a reference node, only if the read keeps matching the reference for lookahead
     * more bases.
     */
    private boolean farJoinAllowed(final int c, final double dist, final byte[] seq, final int pos, final boolean forward) {
        if (dist > maxJump) {
            return false;
        }
        final int r = referenceOffsets[c];
        if (r < 0) {
            return true;
        } else if (forward) {
            return pos + k + lookahead <= seq.length && r + k + lookahead <= reference.length &&
                    Arrays.equals(seq, pos, pos + k + lookahead, reference, r, r + k + lookahead);
        } else {
            return pos - lookahead >= 0 && r - lookahead >= 0 &&
                    Arrays.equals(seq, pos - lookahead, pos + k, reference, r - lookahead, r + k);
        }
    }

    /** the node, if any, that a read's first placed kmer joins */
    private int pickInitial(final long kmer, final double expected, final byte[] seq, final int pos) {
        final IntArrayList candidates = index.get(kmer);
        if (candidates == null || candidates.isEmpty()) {
            return -1;
        } else if (isUnique(kmer)) {
            final int c = candidates.getInt(0);
            final double dist = Math.abs(ranks[c] - expected);
            return dist <= delta || farJoinAllowed(c, dist, seq, pos, true) ? c : -1;
        }
        final IntArrayList within = new IntArrayList(candidates.size());
        for (int i = 0; i < candidates.size(); i++) {
            if (Math.abs(ranks[candidates.getInt(i)] - expected) <= delta) {
                within.add(candidates.getInt(i));
            }
        }
        if (within.size() <= 1) {
            return within.isEmpty() ? -1 : within.getInt(0);
        }
        // several copies of a repeated kmer: prefer the one the read agrees with longest on both sides
        int best = -1, bestAgreement = 0;
        double bestDist = 0;
        for (int i = 0; i < within.size(); i++) {
            final int c = within.getInt(i);
            final double dist = Math.abs(ranks[c] - expected);
            final int agreement = agreement(c, seq, pos, true, AGREEMENT_CAP) + agreement(c, seq, pos, false, AGREEMENT_CAP);
            if (best < 0 || agreement > bestAgreement || (agreement == bestAgreement && (dist < bestDist || (dist == bestDist && c < best)))) {
                best = c;
                bestAgreement = agreement;
                bestDist = dist;
            }
        }
        return best;
    }

    /**
     * How many read bases beyond the kmer at pos agree with what follows node c in the threading direction: the
     * reference for a reference node, otherwise a graph path that matches the read (and then the reference).
     */
    private int agreement(final int c, final byte[] seq, final int pos, final boolean forward, final int cap) {
        int n = 0;
        final int r = referenceOffsets[c];
        if (r >= 0) {
            int i = forward ? pos + k : pos - 1, j = forward ? r + k : r - 1;
            final int step = forward ? 1 : -1;
            while (i >= 0 && i < seq.length && j >= 0 && j < reference.length && seq[i] == reference[j] && n < cap) {
                i += step;
                j += step;
                n++;
            }
            return n;
        }
        int current = c, i = forward ? pos + k : pos - 1;
        while (i >= 0 && i < seq.length && n < cap) {
            final IntArrayList adjacent = forward ? successors[current] : predecessors[current];
            int next = -1;
            for (int a = 0; a < adjacent.size(); a++) {
                final int u = adjacent.getInt(a);
                if ((forward ? lastBases[u] : firstBases[u]) == seq[i]) {
                    next = u;
                    break;
                }
            }
            if (next < 0) {
                break;
            }
            current = next;
            n++;
            if (referenceOffsets[current] >= 0) {
                return n + agreement(current, seq, forward ? i - k + 1 : i, forward, cap - n);
            }
            i += forward ? 1 : -1;
        }
        return n;
    }

    /** the node, if any, that a read kmer joins after (before, threading backwards) node neighbor */
    private int pick(final long kmer, final int neighbor, final double q, final double expected, final boolean forward,
                     final byte[] seq, final int pos) {
        final int c = pickCandidate(kmer, neighbor, q, expected, forward, seq, pos);
        if (c < 0 || referenceOffsets[c] < 0) {
            return c;
        }
        final double impliedDeletion = forward ? ranks[c] - expected : expected - ranks[c];
        if (impliedDeletion > 1 && !isUnique(kmer)) {
            // Joining c would make the read skip reference bases.  In a repeat a later copy of the kmer is often
            // available although the read really carries an insertion, so demand that the read keep agreeing with the
            // reference after the jump and, where the read leaves the reference, longer than under a substitution or
            // an insertion.  Otherwise the kmer stays novel: the read rejoins a few kmers later, or its dangling end is
            // recovered.
            if (agreement(c, seq, pos, forward, AGREEMENT_CAP) < k) {
                return -1;
            }
            if (referenceOffsets[neighbor] >= 0 && !deletionBeatsAlternatives(neighbor, c, seq, pos, forward)) {
                return -1;
            }
        }
        return c;
    }

    private int pickCandidate(final long kmer, final int neighbor, final double q, final double expected, final boolean forward,
                              final byte[] seq, final int pos) {
        final IntArrayList candidates = index.get(kmer);
        if (candidates == null || candidates.isEmpty()) {
            return -1;
        }
        final double neighborRank = ranks[neighbor];
        final IntArrayList ok = new IntArrayList(candidates.size());
        for (int i = 0; i < candidates.size(); i++) {
            final int c = candidates.getInt(i);
            if (forward ? ranks[c] > neighborRank : ranks[c] < neighborRank) {
                ok.add(c);
            }
        }
        if (ok.isEmpty()) {
            return -1;
        } else if (isUnique(kmer)) {
            final int c = ok.getInt(0);
            final double dist = distance(c, q, expected);
            return dist <= delta || farJoinAllowed(c, dist, seq, pos, forward) ? c : -1;
        }
        // a copy that continues an existing edge from the neighbor
        final IntArrayList adjacent = forward ? successors[neighbor] : predecessors[neighbor];
        int best = -1;
        double bestDist = 0;
        for (int i = 0; i < ok.size(); i++) {
            final int c = ok.getInt(i);
            if (adjacent.contains(c)) {
                final double dist = Math.abs(ranks[c] - expected);
                if (best < 0 || dist < bestDist || (dist == bestDist && c < best)) {
                    best = c;
                    bestDist = dist;
                }
            }
        }
        if (best >= 0) {
            return best;
        }
        // several copies of a repeated kmer: prefer the one the rest of the read agrees with longest
        int bestAgreement = 0;
        for (int i = 0; i < ok.size(); i++) {
            final int c = ok.getInt(i);
            final double dist = distance(c, q, expected);
            if (dist > delta) {
                continue;
            }
            final int agreement = agreement(c, seq, pos, forward, AGREEMENT_CAP);
            if (best < 0 || agreement > bestAgreement || (agreement == bestAgreement && (dist < bestDist || (dist == bestDist && c < best)))) {
                best = c;
                bestAgreement = agreement;
                bestDist = dist;
            }
        }
        return best;
    }

    private double distance(final int c, final double q, final double expected) {
        return Math.min(Math.abs(ranks[c] - expected), Math.abs(ranks[c] - q));
    }

    /**
     * Does the read agree with the reference after jumping from reference node neighbor to reference node c for
     * longer than it would after a substitution or an insertion of up to delta bases at the point where it leaves
     * the reference?
     */
    private boolean deletionBeatsAlternatives(final int neighbor, final int c, final byte[] seq, final int pos, final boolean forward) {
        final int r = referenceOffsets[neighbor];
        final int jump = agreement(c, seq, pos, forward, AGREEMENT_CAP);
        // the first read base that disagrees with the reference continuation of the neighbor, and its reference offset
        final int i0 = forward ? pos + k - 1 : pos, j0 = forward ? r + k : r - 1, step = forward ? 1 : -1;
        int alternative = referenceAgreement(seq, i0 + step, j0 + step, step);
        for (int s = 1; s <= (int) delta; s++) {
            alternative = Math.max(alternative, referenceAgreement(seq, i0 + step * s, j0, step));
        }
        return jump >= k && jump > alternative;
    }

    private int referenceAgreement(final byte[] seq, final int start, final int referenceStart, final int step) {
        int n = 0;
        for (int i = start, j = referenceStart; i >= 0 && i < seq.length && j >= 0 && j < reference.length && seq[i] == reference[j] && n < DELETION_AGREEMENT_CAP; i += step, j += step) {
            n++;
        }
        return n;
    }

    /**
     * Thread seq[a, b), all of whose bases are A, C, G or T, into the graph.
     * @param projection the expected window offset of each read base
     * @return the node of each kmer of seq[a, b)
     */
    int[] thread(final byte[] seq, final double[] projection, final int a, final int b) {
        final int n = b - a - k + 1;
        final long[] codes = new long[n];
        final long mask = k == 32 ? -1L : (1L << (2 * k)) - 1;
        long code = encode(seq, a, k);
        codes[0] = code;
        for (int j = 1; j < n; j++) {
            code = ((code << 2) | BaseUtils.simpleBaseToBaseIndex(seq[a + j + k - 1])) & mask;
            codes[j] = code;
        }

        int start = -1, startNode = -1;
        for (int j = 0; j < n; j++) {
            final int v = pickInitial(codes[j], projection[a + j], seq, a + j);
            if (v >= 0) {
                start = j;
                startNode = v;
                break;
            }
        }
        if (start < 0) {
            start = 0;
            startNode = newNode(codes[0], projection[a], -1, seq[a], seq[a + k - 1], true);
        }
        final int[] path = new int[n];
        path[start] = startNode;

        // the last kmer on a reference node, as (kmer index, reference offset), from which to extrapolate positions
        int anchorIndex = referenceOffsets[startNode] >= 0 ? start : -1, anchorOffset = referenceOffsets[startNode];
        int previous = startNode;
        for (int j = start + 1; j < n; j++) {
            final double q = projection[a + j];
            final double expected = anchorIndex >= 0 ? anchorOffset + (j - anchorIndex) : q;
            int v = pick(codes[j], previous, q, expected, true, seq, a + j);
            if (v < 0) {
                v = newNode(codes[j], ranks[previous] + EPSILON, -1, seq[a + j], seq[a + j + k - 1], true);
            }
            addEdge(previous, v, 1);
            path[j] = v;
            previous = v;
            if (referenceOffsets[v] >= 0) {
                anchorIndex = j;
                anchorOffset = referenceOffsets[v];
            }
        }
        anchorIndex = referenceOffsets[startNode] >= 0 ? start : -1;
        anchorOffset = referenceOffsets[startNode];
        int next = startNode;
        for (int j = start - 1; j >= 0; j--) {
            final double q = projection[a + j];
            final double expected = anchorIndex >= 0 ? anchorOffset - (anchorIndex - j) : q;
            int v = pick(codes[j], next, q, expected, false, seq, a + j);
            if (v < 0) {
                v = newNode(codes[j], ranks[next] - EPSILON, -1, seq[a + j], seq[a + j + k - 1], true);
            }
            addEdge(v, next, 1);
            path[j] = v;
            next = v;
            if (referenceOffsets[v] >= 0) {
                anchorIndex = j;
                anchorOffset = referenceOffsets[v];
            }
        }
        return path;
    }

    /**
     * Follow the heaviest edges with at least minWeight reads from node v (forward or backward) until a reference node.
     * @return the walk, or {@code null} if it reaches a dead end first
     */
    Walk walk(final int v, final boolean forward, final int minWeight) {
        final long key = 2L * v + (forward ? 1 : 0);
        final Walk memo = walks.get(key);
        if (memo != null) {
            return memo == NO_WALK ? null : memo;
        }
        final IntArrayList walked = new IntArrayList();
        Walk result = NO_WALK;
        int current = v;
        while (true) {
            final IntArrayList adjacent = forward ? successors[current] : predecessors[current];
            final IntArrayList weights = forward ? successorWeights[current] : predecessorWeights[current];
            int next = -1, nextWeight = 0;
            for (int i = 0; i < adjacent.size(); i++) {
                final int u = adjacent.getInt(i), w = weights.getInt(i);
                if (w >= minWeight && (next < 0 || w > nextWeight || (w == nextWeight && u < next))) {
                    next = u;
                    nextWeight = w;
                }
            }
            if (next < 0) {
                break;
            }
            walked.add(next);
            if (referenceOffsets[next] >= 0) {
                result = new Walk(walked.toIntArray(), next);
                break;
            }
            current = next;
        }
        walks.put(key, result);
        return result == NO_WALK ? null : result;
    }
}
