package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import com.google.common.annotations.VisibleForTesting;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.AdaptiveChainPruner;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Finds the haplotypes of one assembly region from reads threaded through a {@link PositionalGraph}:
 * <ol>
 *     <li>Each read's path becomes events relative to the reference between consecutive reference nodes.  A dangling
 *     end (a read end that leaves the reference and doesn't come back) is completed by comparison with the reference or
 *     by walking the graph back to the reference; one shorter than --min-dangling-branch-length kmers only supports
 *     events that other reads make.  An event is kept if enough reads carry it, counting every position the same indel
 *     takes in different reads: at least --min-pruning and, with adaptive pruning, enough that sequencing errors are an
 *     unlikely explanation at the depth of the kmer before or after it.</li>
 *     <li>Kept events closer than a kmer, or in the same tandem repeat, form a site.  A site's alleles are the
 *     combinations of its events that reads spanning it carry, plus each event in the context of the site's events
 *     that most of its reads carry.</li>
 *     <li>Each fragment gets, per site, the allele its reads spanning the site carry or else the set of alleles that
 *     best match its bases.</li>
 *     <li>Haplotypes are chosen so that every recurring fragment pattern is consistent with one and every allele is on
 *     one, most frequent patterns first.  Sites no pattern decides take the site's most frequent allele.</li>
 * </ol>
 * Allele sets are bit masks, so a site keeps at most {@link #MAX_ALLELES_PER_SITE} alleles.
 */
final class PositionalHaplotypeFinder {
    static final int MAX_ALLELES_PER_SITE = 63;

    /** a haplotype as its events in window coordinates, sorted, with the number of fragments that support it */
    record AssembledHaplotype(List<PositionalEvent> events, int support) { }

    private record Site(int start, int end, List<PositionalEvent> events) { }

    /** a read end left dangling: the read, its segment and path, and the reference anchor it leaves from */
    private record DanglingEnd(int read, int segmentStart, int segmentEnd, int[] path, boolean tail, int anchorOffset, int anchorReference) { }

    /** reference anchors of a read: read offsets and the reference offsets of their kmers' nodes, in read order */
    private record Anchors(int[] readOffsets, int[] referenceOffsets) {
        int size() {
            return readOffsets.length;
        }
    }

    private final PositionalAssemblyParameters params;
    private final int k;
    private final byte[] windowBases;
    private final String window;
    private final List<PositionalRead> reads;

    // statistics of the last assembly
    int graphSize, supportedEvents, siteCount, patternCount, unresolvedDanglingEnds, truncatedSites;
    /** extra positions of supported events, (read, site) pairs given the allele the read spans, and pairs rescued by alignment */
    int extraEventPositions, spanningReadAlleles, alignedReadAlleles;

    /**
     * @param reads the reads, in any order: they are threaded in {@link PositionalRead#THREADING_ORDER}
     */
    PositionalHaplotypeFinder(final PositionalAssemblyParameters params, final byte[] windowBases, final List<PositionalRead> reads) {
        this.params = params;
        this.k = params.k();
        this.windowBases = windowBases;
        this.window = new String(windowBases, StandardCharsets.US_ASCII);
        this.reads = reads.stream().sorted(PositionalRead.THREADING_ORDER).toList();
    }

    private List<PositionalEvent> normalize(final int position, final String ref, final String alt, final int floor) {
        return PositionalEvent.normalize(window, position, ref, alt, floor, params.eventScores());
    }

    /** @return the haplotypes, most supported first, starting with the reference (no events) if it is the most supported */
    List<AssembledHaplotype> findHaplotypes() {
        final PositionalGraph graph = new PositionalGraph(windowBases, k, params.delta(), params.maxJump(), params.lookahead());
        for (final PositionalRead read : reads) {
            read.paths = new int[read.segments.size()][];
            for (int s = 0; s < read.segments.size(); s++) {
                read.paths[s] = graph.thread(read.bases, read.projection, read.segments.get(s)[0], read.segments.get(s)[1]);
            }
        }
        graphSize = graph.size();

        // ---- events from read paths ---------------------------------------------------------------------------------
        final int nReads = reads.size();
        final List<Set<PositionalEvent>> supportedByRead = new ArrayList<>(nReads);   // every event of each read
        final List<Set<PositionalEvent>> onPath = new ArrayList<>(nReads);            // events between reference anchors
        final List<Set<PositionalEvent>> shortEndsByRead = new ArrayList<>(nReads);   // events of dangling ends too short to make events
        final List<Anchors> anchors = new ArrayList<>(nReads);
        final List<DanglingEnd> danglingEnds = new ArrayList<>();
        for (int r = 0; r < nReads; r++) {
            final PositionalRead read = reads.get(r);
            final Set<PositionalEvent> events = new HashSet<>();
            final Set<PositionalEvent> shortEndEvents = new HashSet<>();
            final Set<PositionalEvent> pathEvents = new HashSet<>();
            final IntArrayList flatOffsets = new IntArrayList(), flatReferences = new IntArrayList();
            int floor = 0;
            for (int s = 0; s < read.segments.size(); s++) {
                final int a = read.segments.get(s)[0], b = read.segments.get(s)[1];
                final int[] path = read.paths[s];
                final IntArrayList offsets = new IntArrayList(), references = new IntArrayList();
                for (int j = 0; j < path.length; j++) {
                    if (graph.isReference(path[j])) {
                        offsets.add(a + j);
                        references.add(graph.referenceOffset(path[j]));
                    }
                }
                for (int x = 0; x + 1 < offsets.size(); x++) {
                    final int i1 = offsets.getInt(x), r1 = references.getInt(x), i2 = offsets.getInt(x + 1), r2 = references.getInt(x + 1);
                    if (i2 - i1 == 1 && r2 - r1 == 1) {
                        continue;
                    }
                    final List<PositionalEvent> segmentEvents = normalize(r1, window.substring(r1, r2 + k), read.basesString(i1, i2 + k), floor);
                    floor = eventsEnd(segmentEvents, floor);
                    pathEvents.addAll(segmentEvents);
                    events.addAll(segmentEvents);
                }
                flatOffsets.addAll(offsets);
                flatReferences.addAll(references);
                if (offsets.isEmpty()) {
                    continue;
                }

                // dangling tail: complete against the reference if the read is aligned there without clipping or indels.
                // As in a de Bruijn graph, a dangling end of fewer than minDanglingBranchLength kmers can't make an event
                // by itself, but it supports an event that other reads make.
                final int iL = offsets.getInt(offsets.size() - 1), rL = references.getInt(references.size() - 1);
                final int tailKmers = a + path.length - 1 - iL;
                if (tailKmers > 0) {
                    final boolean shortTail = tailKmers < params.minDanglingBranchLength();
                    boolean completed = false;
                    if (!anyClipped(read, iL, b) && unitSteps(read, iL, b - 1) && read.projection[iL] == rL && rL + (b - iL) <= windowBases.length) {
                        final List<PositionalEvent> tailEvents = normalize(rL, window.substring(rL, rL + b - iL), read.basesString(iL, b), floor);
                        if (tailEvents.size() <= params.maxTailMismatches()) {
                            (shortTail ? shortEndEvents : events).addAll(tailEvents);
                            completed = true;
                        }
                    }
                    if (!completed && !shortTail) {
                        danglingEnds.add(new DanglingEnd(r, a, b, path, true, iL, rL));
                    }
                }
                // dangling head
                final int iF = offsets.getInt(0), rF = references.getInt(0);
                final int headKmers = iF - a;
                if (headKmers > 0) {
                    final boolean shortHead = headKmers < params.minDanglingBranchLength();
                    final int headStart = rF - headKmers;
                    boolean completed = false;
                    if (!anyClipped(read, a, iF) && unitSteps(read, a, iF) && read.projection[iF] == rF && headStart >= 0) {
                        final List<PositionalEvent> headEvents = normalize(headStart, window.substring(headStart, rF + k), read.basesString(a, iF + k), 0);
                        if (headEvents.size() <= params.maxTailMismatches()) {
                            (shortHead ? shortEndEvents : events).addAll(headEvents);
                            completed = true;
                        }
                    }
                    if (!completed && !shortHead) {
                        danglingEnds.add(new DanglingEnd(r, a, b, path, false, iF, rF));
                    }
                }
            }
            supportedByRead.add(events);
            shortEndsByRead.add(shortEndEvents);
            onPath.add(pathEvents);
            anchors.add(new Anchors(flatOffsets.toIntArray(), flatReferences.toIntArray()));
        }

        // ---- dangling ends that rejoin the reference through the graph ---------------------------------------------
        unresolvedDanglingEnds = 0;
        for (final DanglingEnd end : danglingEnds) {
            final PositionalRead read = reads.get(end.read);
            final PositionalGraph.Walk walk = graph.walk(end.path[end.tail ? end.path.length - 1 : 0], end.tail, 1);
            if (walk == null) {
                unresolvedDanglingEnds++;
                continue;
            }
            final int rejoin = graph.referenceOffset(walk.referenceNode());
            final StringBuilder alt = new StringBuilder();
            if (end.tail) {
                if (rejoin + k < end.anchorReference) {
                    continue;
                }
                alt.append(read.basesString(end.anchorOffset, end.segmentEnd));
                for (final int node : walk.nodes()) {
                    alt.append((char) graph.lastBase(node));
                }
                supportedByRead.get(end.read).addAll(normalize(end.anchorReference, window.substring(end.anchorReference, rejoin + k), alt.toString(), 0));
            } else {
                if (rejoin > end.anchorReference) {
                    continue;
                }
                for (int i = walk.nodes().length - 1; i >= 0; i--) {
                    alt.append((char) graph.firstBase(walk.nodes()[i]));
                }
                alt.append(read.basesString(end.segmentStart, end.anchorOffset + k));
                supportedByRead.get(end.read).addAll(normalize(rejoin, window.substring(rejoin, end.anchorReference + k), alt.toString(), 0));
            }
        }

        // ---- events with enough support ---------------------------------------------------------------------------
        // A read's indel can't be left-aligned past the read's previous event, so the same indel can take different
        // positions in different reads.  Support is counted per event left-aligned as far as the reference allows.
        final Map<PositionalEvent, PositionalEvent> leftAligned = new HashMap<>();
        final Map<PositionalEvent, Integer> support = new HashMap<>();
        for (final Set<PositionalEvent> events : supportedByRead) {
            final Set<PositionalEvent> keys = new HashSet<>();
            for (final PositionalEvent event : events) {
                keys.add(leftAligned.computeIfAbsent(event, e -> e.leftAligned(window)));
            }
            for (final PositionalEvent key : keys) {
                support.merge(key, 1, Integer::sum);
            }
        }
        // short dangling ends add support only to events that other reads make
        for (int r = 0; r < nReads; r++) {
            final Set<PositionalEvent> keys = new HashSet<>();
            for (final PositionalEvent event : supportedByRead.get(r)) {
                keys.add(leftAligned.get(event));
            }
            final Set<PositionalEvent> shortEndKeys = new HashSet<>();
            for (final PositionalEvent event : shortEndsByRead.get(r)) {
                final PositionalEvent key = event.leftAligned(window);
                if (support.containsKey(key) && !keys.contains(key)) {
                    leftAligned.putIfAbsent(event, key);
                    shortEndKeys.add(key);
                }
            }
            for (final PositionalEvent key : shortEndKeys) {
                support.merge(key, 1, Integer::sum);
            }
        }
        final Set<PositionalEvent> supportedKeys = new HashSet<>(supportedEvents(support));
        supportedEvents = supportedKeys.size();
        // every position of a supported event that reads carry: they share a footprint, hence a site
        final Set<PositionalEvent> positions = new TreeSet<>();
        leftAligned.forEach((event, key) -> {
            if (supportedKeys.contains(key)) {
                positions.add(event);
            }
        });
        final List<PositionalEvent> supported = new ArrayList<>(positions);
        extraEventPositions = supported.size() - supportedKeys.size();

        // ---- sites ----------------------------------------------------------------------------------------------------
        final List<Site> sites = findSites(supported);
        siteCount = sites.size();
        final int nSites = sites.size();

        // ---- alleles ---------------------------------------------------------------------------------------------------
        truncatedSites = 0;
        final List<List<String>> siteAlleles = new ArrayList<>(nSites);
        final List<byte[][]> siteAlleleBases = new ArrayList<>(nSites);
        final List<List<List<PositionalEvent>>> siteAlleleEvents = new ArrayList<>(nSites);
        // per site, the allele each read carries across it, if the read spans the site and that sequence is an allele
        final List<int[]> carriedAlleles = new ArrayList<>(nSites);
        for (final Site site : sites) {
            final String[] readKeys = new String[nReads];
            List<String> alleles = siteAlleles(site, onPath, anchors, readKeys);
            if (alleles.size() > MAX_ALLELES_PER_SITE) {
                alleles = new ArrayList<>(alleles.subList(0, MAX_ALLELES_PER_SITE));
                truncatedSites++;
            }
            siteAlleles.add(alleles);
            final Map<String, Integer> alleleIndex = new HashMap<>();
            for (int a = 0; a < alleles.size(); a++) {
                alleleIndex.putIfAbsent(alleles.get(a), a);
            }
            final int[] carried = new int[nReads];
            for (int r = 0; r < nReads; r++) {
                carried[r] = readKeys[r] == null ? -1 : alleleIndex.getOrDefault(readKeys[r], -1);
            }
            carriedAlleles.add(carried);
            final byte[][] bases = new byte[alleles.size()][];
            final List<List<PositionalEvent>> events = new ArrayList<>(alleles.size());
            final String refKey = window.substring(site.start, site.end);
            for (int a = 0; a < alleles.size(); a++) {
                bases[a] = alleles.get(a).getBytes(StandardCharsets.US_ASCII);
                events.add(normalize(site.start, refKey, alleles.get(a), site.start));
            }
            siteAlleleBases.add(bases);
            siteAlleleEvents.add(events);
        }

        // ---- fragment patterns -------------------------------------------------------------------------------------
        spanningReadAlleles = 0;
        alignedReadAlleles = 0;
        final Map<String, long[]> fragmentMasks = new LinkedHashMap<>();
        final Map<String, boolean[]> fragmentSeen = new HashMap<>();
        for (int r = 0; r < nReads; r++) {
            final Anchors flat = anchors.get(r);
            if (flat.size() == 0) {
                continue;
            }
            final PositionalRead read = reads.get(r);
            final int low = flat.referenceOffsets[0] - flat.readOffsets[0];
            final int high = flat.referenceOffsets[flat.size() - 1] + (read.length() - flat.readOffsets[flat.size() - 1]);
            final long[] masks = fragmentMasks.computeIfAbsent(read.fragment, f -> new long[nSites]);
            final boolean[] seen = fragmentSeen.computeIfAbsent(read.fragment, f -> new boolean[nSites]);
            for (int si = 0; si < nSites; si++) {
                final Site site = sites.get(si);
                if (high < site.start || low > site.end) {
                    continue;
                }
                // a read spanning the site carries the allele its supported events make there, whatever errors it has
                final int carried = carriedAlleles.get(si)[r];
                if (carried >= 0) {
                    spanningReadAlleles++;
                }
                final long consistent = carried >= 0 ? 1L << carried : consistency(read, flat, site, siteAlleleBases.get(si));
                if (consistent != 0) {
                    // mates that disagree leave an empty set, which counts as uninformative
                    masks[si] = seen[si] ? masks[si] & consistent : consistent;
                    seen[si] = true;
                }
            }
        }
        final Map<Pattern, Integer> patternCountMap = new LinkedHashMap<>();
        for (final long[] masks : fragmentMasks.values()) {
            if (Arrays.stream(masks).anyMatch(m -> m != 0)) {
                patternCountMap.merge(new Pattern(masks), 1, Integer::sum);
            }
        }
        patternCount = patternCountMap.size();
        final List<long[]> patterns = new ArrayList<>(patternCount);
        final int[] counts = new int[patternCount];
        for (final Map.Entry<Pattern, Integer> entry : patternCountMap.entrySet()) {
            counts[patterns.size()] = entry.getValue();
            patterns.add(entry.getKey().masks);
        }

        // ---- haplotype cover -------------------------------------------------------------------------------------
        final Cover cover = cover(nSites, siteAlleles, patterns, counts);
        final List<AssembledHaplotype> result = new ArrayList<>(cover.alleles.size());
        final Set<String> seen = new HashSet<>();
        for (int h = 0; h < cover.alleles.size(); h++) {
            final int[] alleles = cover.alleles.get(h);
            final StringBuilder sequence = new StringBuilder(window.length() + 16);
            final List<PositionalEvent> events = new ArrayList<>();
            int next = 0;
            for (int si = 0; si < nSites; si++) {
                final Site site = sites.get(si);
                sequence.append(window, next, site.start).append(siteAlleles.get(si).get(alleles[si]));
                next = site.end;
                events.addAll(siteAlleleEvents.get(si).get(alleles[si]));
            }
            sequence.append(window, next, window.length());
            if (seen.add(sequence.toString())) {
                Collections.sort(events);
                result.add(new AssembledHaplotype(events, cover.supports[h]));
            }
        }
        return result;
    }

    private static int eventsEnd(final List<PositionalEvent> events, final int floor) {
        int end = floor;
        for (final PositionalEvent event : events) {
            end = Math.max(end, event.position() + Math.max(1, event.ref().length()));
        }
        return end;
    }

    private static boolean anyClipped(final PositionalRead read, final int from, final int to) {
        for (int i = from; i < to; i++) {
            if (read.clipped[i]) {
                return true;
            }
        }
        return false;
    }

    /** are the bases at offsets [from, to] aligned to consecutive reference positions? */
    private static boolean unitSteps(final PositionalRead read, final int from, final int to) {
        for (int i = from; i < to; i++) {
            if (read.projection[i + 1] - read.projection[i] != 1) {
                return false;
            }
        }
        return true;
    }

    /**
     * Read depths as a de Bruijn graph's branch vertices see them.  left[p] counts the reads with usable bases at window
     * offsets p - k through p (the kmer before an event at p, and the event's first base); right[e] counts the reads with
     * usable bases at e - 1 through e + k - 1 (the last base of an event ending at e, and the kmer after it).
     */
    private int[][] kmerDepths() {
        final int length = windowBases.length;
        final int[] left = new int[length + 1], right = new int[length + 1];
        for (final PositionalRead read : reads) {
            for (final int[] segment : read.segments) {
                final int low = (int) Math.ceil(read.projection[segment[0]]), high = (int) Math.floor(read.projection[segment[1] - 1]);
                addRange(left, low + k, high);
                addRange(right, low + 1, high - k + 1);
            }
        }
        for (int i = 1; i <= length; i++) {
            left[i] += left[i - 1];
            right[i] += right[i - 1];
        }
        return new int[][] {left, right};
    }

    /** add 1 at offsets from through to, inclusive, of a difference array */
    private static void addRange(final int[] difference, final int from, final int to) {
        final int a = Math.max(0, from), b = Math.min(difference.length - 2, to);
        if (a <= b) {
            difference[a]++;
            difference[b + 1]--;
        }
    }

    /**
     * The events with at least minSupport reads and, with adaptive pruning, the log odds against sequencing error that
     * {@link AdaptiveChainPruner} requires of a chain: the event's reads against the reads through the kmer before the
     * event or, if that is better, the kmer after it, with an error rate estimated, as the pruner does, from the events
     * that fail the test at the initial error rate.
     */
    private List<PositionalEvent> supportedEvents(final Map<PositionalEvent, Integer> support) {
        final List<PositionalEvent> result = new ArrayList<>();
        if (!params.adaptivePruning()) {
            support.forEach((event, count) -> {
                if (count >= params.minSupport()) {
                    result.add(event);
                }
            });
            return result;
        }
        final int[][] depths = kmerDepths();
        long errors = 0;
        for (final Map.Entry<PositionalEvent, Integer> entry : support.entrySet()) {
            if (!passesPruning(entry.getKey(), entry.getValue(), depths, params.initialErrorRate())) {
                errors += entry.getValue();
            }
        }
        // the pruner divides by the total multiplicity of the graph's edges, that is, by the reads' kmer-to-kmer steps
        long steps = 0;
        for (final PositionalRead read : reads) {
            for (final int[] segment : read.segments) {
                steps += segment[1] - segment[0] - k;
            }
        }
        final double errorRate = Math.min(1.0, (double) errors / Math.max(1, steps));
        support.forEach((event, count) -> {
            if (count >= params.minSupport() && passesPruning(event, count, depths, errorRate)) {
                result.add(event);
            }
        });
        return result;
    }

    private boolean passesPruning(final PositionalEvent event, final int count, final int[][] depths, final double errorRate) {
        final int last = windowBases.length - 1;
        final int start = Math.min(last, Math.max(0, event.position()));
        final int end = Math.min(last, Math.max(0, event.position() + event.ref().length()));
        final int left = Math.max(count, depths[0][start]), right = Math.max(count, depths[1][end]);
        return Math.max(Mutect2Engine.logLikelihoodRatio(left - count, count, errorRate),
                Mutect2Engine.logLikelihoodRatio(right - count, count, errorRate)) >= params.pruningLogOddsThreshold();
    }

    /** merge supported events whose footprints (widened over tandem repeats) are closer than a kmer */
    private List<Site> findSites(final List<PositionalEvent> supported) {
        final int n = supported.size();
        final int[][] footprints = new int[n][];
        final Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) {
            footprints[i] = footprint(supported.get(i));
            order[i] = i;
        }
        Arrays.sort(order, Comparator.<Integer>comparingInt(i -> footprints[i][0]).thenComparingInt(i -> footprints[i][1]));
        final List<int[]> bounds = new ArrayList<>();
        final List<List<PositionalEvent>> events = new ArrayList<>();
        for (final int i : order) {
            final int s = footprints[i][0], t = footprints[i][1];
            if (!bounds.isEmpty() && s < bounds.get(bounds.size() - 1)[1] + k) {
                final int[] last = bounds.get(bounds.size() - 1);
                last[1] = Math.max(last[1], t);
                events.get(events.size() - 1).add(supported.get(i));
            } else {
                bounds.add(new int[] {s, t});
                events.add(new ArrayList<>(List.of(supported.get(i))));
            }
        }
        final List<Site> sites = new ArrayList<>(bounds.size());
        for (int i = 0; i < bounds.size(); i++) {
            sites.add(new Site(Math.max(0, bounds.get(i)[0]), Math.min(windowBases.length, bounds.get(i)[1]), events.get(i)));
        }
        return sites;
    }

    /** the window interval an event occupies: a substitution its bases, an indel its tandem repeat tract plus a base on each side */
    private int[] footprint(final PositionalEvent event) {
        final int p = event.position();
        if (event.isSubstitution()) {
            return new int[] {p, p + event.ref().length()};
        }
        final String indel = event.ref().isEmpty() ? event.alt() : event.ref();
        final int end = repeatTractEnd(p, indel, event.end());
        return new int[] {p - 1, end + 1};
    }

    /** extend a tandem repeat of unitSequence's minimal period, phased at position, from start while the window continues it */
    private int repeatTractEnd(final int position, final String unitSequence, final int start) {
        final int length = unitSequence.length();
        int period = length;
        for (int p = 1; p <= length; p++) {
            if (length % p == 0 && unitSequence.substring(0, p).repeat(length / p).equals(unitSequence)) {
                period = p;
                break;
            }
        }
        int j = start;
        while (j < windowBases.length && windowBases[j] == unitSequence.charAt((j - position) % period)) {
            j++;
        }
        return j;
    }

    /**
     * The alleles of a site: the reference, then the sequences that at least minSupport fragments anchored on both
     * sides carry across it (as the reference edited by the read's supported events at the site), most frequent first,
     * then each event in the context of the site's events that most of its reads carry, and each event no allele
     * carries yet.
     *
     * @param readKeys filled with the sequence each read spanning the site carries across it
     */
    private List<String> siteAlleles(final Site site, final List<Set<PositionalEvent>> onPath, final List<Anchors> anchors,
                                     final String[] readKeys) {
        final int s = site.start, t = site.end;
        final Set<PositionalEvent> siteEvents = new HashSet<>(site.events);
        final List<PositionalEvent> sortedSiteEvents = new ArrayList<>(siteEvents);
        Collections.sort(sortedSiteEvents);
        final String refKey = window.substring(s, t);
        final Map<String, Set<String>> fragmentsByKey = new HashMap<>();
        final Map<String, Set<PositionalEvent>> carried = new HashMap<>();
        for (int r = 0; r < reads.size(); r++) {
            final Anchors flat = anchors.get(r);
            int leftOffset = -1, leftReference = -1, rightOffset = -1, rightReference = -1;
            for (int x = 0; x < flat.size(); x++) {
                final int ref = flat.referenceOffsets[x];
                if (ref < s) {
                    leftOffset = flat.readOffsets[x];
                    leftReference = ref;
                } else if (ref >= t) {
                    rightOffset = flat.readOffsets[x];
                    rightReference = ref;
                    break;
                }
            }
            if (leftOffset < 0 || rightOffset < 0) {
                continue;
            }
            final PositionalRead read = reads.get(r);
            final int iS = leftOffset + (s - leftReference), iT = rightOffset - (rightReference - t);
            if (!(0 <= iS && iS <= iT && iT <= read.length())) {
                continue;
            }
            final List<PositionalEvent> readSiteEvents = intersection(sortedSiteEvents, onPath.get(r));
            final String projected = PositionalEvent.apply(refKey, s, readSiteEvents);
            final String key = projected != null ? projected : read.basesString(iS, iT);
            readKeys[r] = key;
            fragmentsByKey.computeIfAbsent(key, x -> new HashSet<>()).add(read.fragment);
            carried.computeIfAbsent(key, x -> new HashSet<>()).addAll(readSiteEvents);
        }
        final List<String> keys = new ArrayList<>();
        keys.add(refKey);
        fragmentsByKey.entrySet().stream()
                .filter(e -> !e.getKey().equals(refKey) && e.getValue().size() >= params.minSupport())
                .sorted(Comparator.<Map.Entry<String, Set<String>>>comparingInt(e -> -e.getValue().size()).thenComparing(Map.Entry::getKey))
                .forEach(e -> keys.add(e.getKey()));
        final Set<PositionalEvent> covered = new HashSet<>();
        for (final String key : keys) {
            covered.addAll(carried.getOrDefault(key, Collections.emptySet()));
        }

        for (final PositionalEvent event : sortedSiteEvents) {
            // the event in the context of the site's events that most of its reads carry, so that it keeps its phase
            // with nearby events when the reads' full patterns are too fragmented (e.g. by stutter) to recur
            int carriers = 0;
            final Map<PositionalEvent, Integer> shared = new HashMap<>();
            for (final Set<PositionalEvent> readEvents : onPath) {
                if (readEvents.contains(event)) {
                    carriers++;
                    for (final PositionalEvent other : sortedSiteEvents) {
                        if (!other.equals(event) && readEvents.contains(other)) {
                            shared.merge(other, 1, Integer::sum);
                        }
                    }
                }
            }
            final Set<PositionalEvent> context = new TreeSet<>();
            context.add(event);
            for (final Map.Entry<PositionalEvent, Integer> entry : shared.entrySet()) {
                if (2 * entry.getValue() > carriers) {
                    context.add(entry.getKey());
                }
            }
            final String inContext = PositionalEvent.apply(refKey, s, context);
            if (inContext != null && !keys.contains(inContext)) {
                keys.add(inContext);
                covered.addAll(context);
            }
            if (!covered.contains(event)) {
                final String alone = PositionalEvent.apply(refKey, s, List.of(event));
                if (alone != null && !keys.contains(alone)) {
                    keys.add(alone);
                }
            }
        }
        return keys;
    }

    private static List<PositionalEvent> intersection(final List<PositionalEvent> sorted, final Set<PositionalEvent> set) {
        final List<PositionalEvent> result = new ArrayList<>();
        for (final PositionalEvent event : sorted) {
            if (set.contains(event)) {
                result.add(event);
            }
        }
        return result;
    }

    /** gap costs, in mismatches, when a read's bases are aligned to a site's alleles */
    private static final int GAP_OPEN = 2, GAP_EXTEND = 1;

    /**
     * The alleles of a site that best match a read's bases next to its nearest reference anchor, as a bit mask, or 0
     * if the read is uninformative: it doesn't reach the site, matches every allele equally, or mismatches all of them.
     * The bases are compared without gaps unless that leaves every allele with too many mismatches, as when the read
     * carries an indel no allele explains (an error, or a pruned event); then they are aligned, allowing one more gap.
     */
    private long consistency(final PositionalRead read, final Anchors flat, final Site site, final byte[][] alleles) {
        final int s = site.start, t = site.end, n = read.length();
        int left = -1, right = -1;
        for (int x = 0; x < flat.size(); x++) {
            if (flat.referenceOffsets[x] < s) {
                left = x;
            } else if (flat.referenceOffsets[x] >= t && right < 0) {
                right = x;
            }
        }
        // each allele's haplotype bases from the anchor through the site, in the direction away from the anchor
        final byte[][] haplotypes = new byte[alleles.length][];
        final boolean forward = left >= 0;
        final int ia;
        if (forward) {
            ia = flat.readOffsets[left];
            final int ra = flat.referenceOffsets[left];
            if (n - ia < s - ra + 1) {
                return 0;
            }
            for (int a = 0; a < alleles.length; a++) {
                final int dl = alleles[a].length - (t - s), end = Math.min(t + dl + 2, windowBases.length + dl);
                haplotypes[a] = new byte[Math.max(0, end - ra)];
                for (int x = 0; x < haplotypes[a].length; x++) {
                    haplotypes[a][x] = haplotypeBase(site, alleles[a], dl, ra + x);
                }
            }
        } else if (right >= 0) {
            ia = flat.readOffsets[right];
            final int ra = flat.referenceOffsets[right];
            if (ia < ra - t + 1) {
                return 0;
            }
            final int lowest = Math.max(0, s - 2);
            for (int a = 0; a < alleles.length; a++) {
                final int dl = alleles[a].length - (t - s), hp = ra + dl;
                haplotypes[a] = new byte[Math.max(0, hp - lowest)];
                for (int x = 0; x < haplotypes[a].length; x++) {
                    haplotypes[a][x] = haplotypeBase(site, alleles[a], dl, hp - 1 - x);
                }
            }
        } else {
            return 0;
        }
        int longest = 0, band = 0;
        for (int a = 0; a < alleles.length; a++) {
            longest = Math.max(longest, haplotypes[a].length);
            band = Math.max(band, Math.abs(alleles[a].length - (t - s)));
        }
        band += 3;
        // the read's bases from the anchor, in the same direction
        final int available = forward ? n - ia : ia, length = Math.min(available, longest + band);
        final byte[] bases = new byte[length], qualities = new byte[length];
        for (int x = 0; x < length; x++) {
            final int i = forward ? ia + x : ia - 1 - x;
            bases[x] = read.bases[i];
            qualities[x] = read.qualities[i];
        }

        final int[] costs = new int[alleles.length];
        for (int a = 0; a < alleles.length; a++) {
            costs[a] = mismatches(haplotypes[a], bases, qualities, params.minBaseQuality());
        }
        int best = Arrays.stream(costs).min().orElse(Integer.MAX_VALUE);
        if (best > params.maxSiteMismatches()) {
            for (int a = 0; a < alleles.length; a++) {
                costs[a] = alignmentCost(haplotypes[a], bases, qualities, length == available, band, params.minBaseQuality());
            }
            best = Arrays.stream(costs).min().orElse(Integer.MAX_VALUE);
            if (best > params.maxSiteMismatches() + GAP_OPEN) {
                return 0;
            }
            alignedReadAlleles++;
        }
        long consistent = 0;
        for (int a = 0; a < alleles.length; a++) {
            if (costs[a] == best) {
                consistent |= 1L << a;
            }
        }
        return consistent == allAlleles(alleles.length) ? 0 : consistent;
    }

    /** mismatches between bases and haplotype, compared without gaps, at bases of at least minBaseQuality */
    private static int mismatches(final byte[] haplotype, final byte[] bases, final byte[] qualities, final int minBaseQuality) {
        int result = 0;
        for (int x = 0; x < Math.min(haplotype.length, bases.length); x++) {
            if (bases[x] != haplotype[x] && qualities[x] >= minBaseQuality) {
                result++;
            }
        }
        return result;
    }

    /**
     * Cost of the best alignment of read bases to haplotype bases that both start next to an anchor: a mismatch at a
     * base of at least minBaseQuality costs 1, and a gap of g bases GAP_OPEN + (g - 1) * GAP_EXTEND.  Read bases after
     * the haplotype's end are free, as are haplotype bases after the read's end if readEnded.  Only cells within band
     * of the diagonal are filled.
     */
    @VisibleForTesting
    static int alignmentCost(final byte[] haplotype, final byte[] bases, final byte[] qualities, final boolean readEnded,
                             final int band, final int minBaseQuality) {
        final int rows = haplotype.length, columns = bases.length, width = columns + 1, infinity = Integer.MAX_VALUE / 4;
        // costs of alignments ending in a base pair, a haplotype base the read lacks, and a read base the haplotype lacks
        final int[] pair = new int[(rows + 1) * width], readGap = new int[(rows + 1) * width], haplotypeGap = new int[(rows + 1) * width];
        Arrays.fill(pair, infinity);
        Arrays.fill(readGap, infinity);
        Arrays.fill(haplotypeGap, infinity);
        pair[0] = 0;
        for (int i = 0; i <= rows; i++) {
            for (int j = Math.max(0, i - band); j <= Math.min(columns, i + band); j++) {
                final int cell = i * width + j;
                if (i > 0 && j > 0) {
                    final int diagonal = cell - width - 1;
                    final int cost = bases[j - 1] != haplotype[i - 1] && qualities[j - 1] >= minBaseQuality ? 1 : 0;
                    pair[cell] = Math.min(pair[diagonal], Math.min(readGap[diagonal], haplotypeGap[diagonal])) + cost;
                }
                if (i > 0) {
                    final int up = cell - width;
                    readGap[cell] = Math.min(Math.min(pair[up], haplotypeGap[up]) + GAP_OPEN, readGap[up] + GAP_EXTEND);
                }
                if (j > 0) {
                    final int back = cell - 1;
                    haplotypeGap[cell] = Math.min(Math.min(pair[back], readGap[back]) + GAP_OPEN, haplotypeGap[back] + GAP_EXTEND);
                }
            }
        }
        int best = infinity;
        for (int j = 0; j <= columns; j++) {
            final int cell = rows * width + j;
            best = Math.min(best, Math.min(pair[cell], Math.min(readGap[cell], haplotypeGap[cell])));
        }
        if (readEnded) {
            for (int i = 0; i <= rows; i++) {
                final int cell = i * width + columns;
                best = Math.min(best, Math.min(pair[cell], Math.min(readGap[cell], haplotypeGap[cell])));
            }
        }
        return best;
    }

    /** base x of the window with the site replaced by the allele, which changes the length by dl */
    private byte haplotypeBase(final Site site, final byte[] allele, final int dl, final int x) {
        if (x < site.start) {
            return windowBases[x];
        } else if (x < site.start + allele.length) {
            return allele[x - site.start];
        } else {
            return windowBases[x - dl];
        }
    }

    private static long allAlleles(final int n) {
        return n >= 64 ? -1L : (1L << n) - 1;
    }

    // ---- haplotype cover ---------------------------------------------------------------------------------------------

    /** a fragment pattern: per site the mask of consistent alleles, 0 where uninformative */
    private record Pattern(long[] masks) {
        @Override
        public boolean equals(final Object other) {
            return other instanceof Pattern p && Arrays.equals(masks, p.masks);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(masks);
        }
    }

    /** haplotypes as an allele per site, and their supports */
    record Cover(List<int[]> alleles, int[] supports) { }

    private record Choice(int[] alleles) {
        @Override
        public boolean equals(final Object other) {
            return other instanceof Choice c && Arrays.equals(alleles, c.alleles);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(alleles);
        }
    }

    private static boolean compatible(final long[] h, final long[] p) {
        for (int i = 0; i < h.length; i++) {
            if (h[i] != 0 && p[i] != 0 && (h[i] & p[i]) == 0) {
                return false;
            }
        }
        return true;
    }

    private static long[] merge(final long[] h, final long[] p) {
        final long[] result = new long[h.length];
        for (int i = 0; i < h.length; i++) {
            result[i] = h[i] == 0 ? p[i] : (p[i] == 0 ? h[i] : h[i] & p[i]);
        }
        return result;
    }

    private static int singleAllele(final long mask) {
        return mask != 0 && Long.bitCount(mask) == 1 ? Long.numberOfTrailingZeros(mask) : -1;
    }

    /** fragment counts per allele */
    private static final class AlleleCounter {
        private final Map<Integer, Integer> counts = new HashMap<>();

        void add(final int key, final int count) {
            counts.merge(key, count, Integer::sum);
        }

        boolean isEmpty() {
            return counts.isEmpty();
        }

        /** alleles by decreasing count, ties to the lower allele (the reference first) */
        List<Integer> mostCommon() {
            final List<Map.Entry<Integer, Integer>> entries = new ArrayList<>(counts.entrySet());
            entries.sort(Comparator.<Map.Entry<Integer, Integer>>comparingInt(e -> -e.getValue()).thenComparingInt(Map.Entry::getKey));
            final List<Integer> keys = new ArrayList<>(entries.size());
            entries.forEach(e -> keys.add(e.getKey()));
            return keys;
        }
    }

    /** Python-style representation of a pattern, the final tie-breaker in the order of patterns */
    private static String representation(final long[] pattern) {
        final StringBuilder result = new StringBuilder("(");
        for (int si = 0; si < pattern.length; si++) {
            if (si > 0) {
                result.append(", ");
            }
            if (pattern[si] == 0) {
                result.append("None");
            } else {
                result.append("frozenset({");
                long mask = pattern[si];
                boolean first = true;
                while (mask != 0) {
                    if (!first) {
                        result.append(", ");
                    }
                    result.append(Long.numberOfTrailingZeros(mask));
                    mask &= mask - 1;
                    first = false;
                }
                result.append("})");
            }
        }
        return result.append(pattern.length == 1 ? ",)" : ")").toString();
    }

    /**
     * Choose haplotypes, as an allele per site, so that every recurring pattern is consistent with one and every allele
     * is on one.  Most frequent patterns seed haplotypes first; a pattern merges into the compatible haplotype it
     * overlaps most, unless that would pin a site of the haplotype to a minority allele, in which case a pattern seen
     * in at least minSupport fragments branches into a new haplotype.
     * @return the haplotypes, sorted by decreasing support, then by their alleles (so the reference first among ties)
     */
    @VisibleForTesting
    Cover cover(final int nSites, final List<List<String>> siteAlleles, final List<long[]> patterns, final int[] counts) {
        final int nPatterns = patterns.size();
        final AlleleCounter[] siteFrequencies = new AlleleCounter[nSites];
        final int[] majority = new int[nSites];
        for (int si = 0; si < nSites; si++) {
            siteFrequencies[si] = new AlleleCounter();
            for (int p = 0; p < nPatterns; p++) {
                final int allele = singleAllele(patterns.get(p)[si]);
                if (allele >= 0) {
                    siteFrequencies[si].add(allele, counts[p]);
                }
            }
            majority[si] = siteFrequencies[si].isEmpty() ? 0 : siteFrequencies[si].mostCommon().get(0);
        }

        final Integer[] ordered = new Integer[nPatterns];
        final int[] singles = new int[nPatterns];
        final String[] representations = new String[nPatterns];
        for (int p = 0; p < nPatterns; p++) {
            ordered[p] = p;
            singles[p] = (int) Arrays.stream(patterns.get(p)).filter(m -> singleAllele(m) >= 0).count();
            representations[p] = representation(patterns.get(p));
        }
        Arrays.sort(ordered, Comparator.<Integer>comparingInt(p -> -counts[p]).thenComparingInt(p -> -singles[p])
                .thenComparing(p -> representations[p]));

        final List<long[]> haplotypes = new ArrayList<>();
        final IntArrayList supports = new IntArrayList();
        for (final int p : ordered) {
            final long[] pattern = patterns.get(p);
            int best = -1, bestOverlap = -1, bestSupport = -1;
            for (int h = 0; h < haplotypes.size(); h++) {
                final long[] haplotype = haplotypes.get(h);
                if (!compatible(haplotype, pattern)) {
                    continue;
                }
                int overlap = 0;
                for (int si = 0; si < nSites; si++) {
                    if (haplotype[si] != 0 && pattern[si] != 0) {
                        overlap++;
                    }
                }
                if (best < 0 || overlap > bestOverlap || (overlap == bestOverlap && supports.getInt(h) > bestSupport)) {
                    best = h;
                    bestOverlap = overlap;
                    bestSupport = supports.getInt(h);
                }
            }
            if (best >= 0) {
                if (!introducesMinority(haplotypes.get(best), pattern, majority)) {
                    haplotypes.set(best, merge(haplotypes.get(best), pattern));
                    supports.set(best, supports.getInt(best) + counts[p]);
                } else if (counts[p] >= params.minSupport()) {
                    haplotypes.add(merge(haplotypes.get(best), pattern));
                    supports.add(counts[p]);
                }
            } else if (counts[p] >= params.minSupport()) {
                haplotypes.add(pattern.clone());
                supports.add(counts[p]);
            }
        }

        final Map<Choice, Integer> finals = new LinkedHashMap<>();
        for (int h = 0; h < haplotypes.size(); h++) {
            finals.merge(new Choice(finalize(haplotypes.get(h), siteAlleles, siteFrequencies)), supports.getInt(h), Integer::sum);
        }
        finals.putIfAbsent(new Choice(new int[nSites]), 0);
        // every allele must be on some haplotype: merge the patterns that determine it, and complete the rest
        for (int si = 0; si < nSites; si++) {
            for (int a = 0; a < siteAlleles.get(si).size(); a++) {
                final int site = si, allele = a;
                if (finals.keySet().stream().anyMatch(c -> c.alleles[site] == allele)) {
                    continue;
                }
                long[] haplotype = new long[nSites];
                haplotype[si] = 1L << a;
                int support = 0;
                for (final int p : ordered) {
                    if (patterns.get(p)[si] == 1L << a && compatible(haplotype, patterns.get(p))) {
                        haplotype = merge(haplotype, patterns.get(p));
                        support += counts[p];
                    }
                }
                final int[] choice = finalize(haplotype, siteAlleles, siteFrequencies);
                choice[si] = a;
                finals.merge(new Choice(choice), support, Integer::sum);
            }
        }

        final List<Map.Entry<Choice, Integer>> entries = new ArrayList<>(finals.entrySet());
        entries.sort((x, y) -> {
            final int bySupport = Integer.compare(y.getValue(), x.getValue());
            return bySupport != 0 ? bySupport : Arrays.compare(x.getKey().alleles, y.getKey().alleles);
        });
        final List<int[]> alleles = new ArrayList<>(entries.size());
        entries.forEach(e -> alleles.add(e.getKey().alleles));
        return new Cover(alleles, entries.stream().mapToInt(Map.Entry::getValue).toArray());
    }

    /** would merging the pattern into the haplotype newly pin one of its sites to alleles excluding the site's majority allele? */
    private static boolean introducesMinority(final long[] haplotype, final long[] pattern, final int[] majority) {
        for (int si = 0; si < haplotype.length; si++) {
            final long a = haplotype[si], b = pattern[si];
            if (b != 0 && (a == 0 || (a & ~b) != 0) && (b & (1L << majority[si])) == 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * An allele for every site of a partially determined haplotype: the determined allele, or else the site's most
     * frequent allele among those allowed.  Fragments that would pin an open site to a minority allele have made their
     * own haplotype instead (see {@link #cover}), and fragments that reach only that site carry no phase, so neither
     * decides it.
     */
    private static int[] finalize(final long[] haplotype, final List<List<String>> siteAlleles, final AlleleCounter[] siteFrequencies) {
        final int[] result = new int[haplotype.length];
        for (int si = 0; si < haplotype.length; si++) {
            final int determined = singleAllele(haplotype[si]);
            if (determined >= 0) {
                result[si] = determined;
                continue;
            }
            final long options = haplotype[si] != 0 ? haplotype[si] : allAlleles(siteAlleles.get(si).size());
            final int choice = firstAllowed(siteFrequencies[si], options);
            result[si] = choice >= 0 ? choice : Long.numberOfTrailingZeros(options);
        }
        return result;
    }

    private static int firstAllowed(final AlleleCounter counter, final long options) {
        for (final int allele : counter.mostCommon()) {
            if ((options & (1L << allele)) != 0) {
                return allele;
            }
        }
        return -1;
    }
}
