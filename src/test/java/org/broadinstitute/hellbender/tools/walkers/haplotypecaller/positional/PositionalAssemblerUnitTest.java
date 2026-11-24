package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.positional;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Reads are simulated from known haplotypes and aligned by the truth, soft-clipped where an event is near a read end as
 * BWA tends to do; the assembled haplotypes must include every true haplotype.
 */
public class PositionalAssemblerUnitTest extends GATKBaseTest {
    private static final int READ_LENGTH = 101, FRAGMENT_LENGTH = 300, CLIP_IF_WITHIN = 12, PADDING = 100, MIN_DANGLING_BRANCH_LENGTH = 4;
    private static final double ERROR_RATE = 0.002;
    private static final String SAMPLE = "tumor", READ_GROUP = "rg", CONTIG = "1";

    private record TrueHaplotype(List<PositionalEvent> events, double weight) { }

    private static String randomBases(final Random random, final int n) {
        final StringBuilder result = new StringBuilder(n);
        for (int i = 0; i < n; i++) {
            result.append("ACGT".charAt(random.nextInt(4)));
        }
        return result.toString();
    }

    @DataProvider(name = "cases")
    public Object[][] cases() {
        final Random random = new Random(7);
        final String left = randomBases(random, 600), right = randomBases(random, 600);
        final List<Object[]> cases = new ArrayList<>();

        // short tandem repeat expansion and contraction: cycles in a de Bruijn graph with k = 10
        final String str = left + "CA".repeat(12) + right;
        final int p = left.length();
        cases.add(new Object[] {"(CA)12 -> (CA)14", str, List.of(new TrueHaplotype(List.of(), 1),
                new TrueHaplotype(List.of(new PositionalEvent(p, "", "CACA")), 1)), p - 20, p + 45, 80});
        cases.add(new Object[] {"(CA)12 -> (CA)10", str, List.of(new TrueHaplotype(List.of(), 1),
                new TrueHaplotype(List.of(new PositionalEvent(p, "CACA", "")), 1)), p - 20, p + 45, 80});

        // tandem duplication of 30 bases
        final String dupRef = left + randomBases(random, 60) + right;
        final int dupPosition = left.length() + 30;
        cases.add(new Object[] {"30-base tandem duplication", dupRef, List.of(new TrueHaplotype(List.of(), 1),
                new TrueHaplotype(List.of(new PositionalEvent(dupPosition, "", dupRef.substring(dupPosition - 30, dupPosition))), 1)),
                dupPosition - 40, dupPosition + 40, 80});

        // a 50-base insertion, seen mostly through soft clips
        final String plain = left + right;
        cases.add(new Object[] {"50-base insertion", plain, List.of(new TrueHaplotype(List.of(), 1),
                new TrueHaplotype(List.of(new PositionalEvent(p, "", randomBases(random, 50))), 1)), p - 30, p + 30, 80});

        // phasing: a germline het SNV and a somatic SNV 60 bases away on the alt haplotype only, in 10% of reads
        final int p1 = p - 30, p2 = p + 30;
        final PositionalEvent germline = new PositionalEvent(p1, plain.substring(p1, p1 + 1), plain.charAt(p1) != 'A' ? "A" : "C");
        final PositionalEvent somatic = new PositionalEvent(p2, plain.substring(p2, p2 + 1), plain.charAt(p2) != 'G' ? "G" : "T");
        cases.add(new Object[] {"somatic SNV phased with a germline het", plain, List.of(new TrueHaplotype(List.of(), 0.5),
                new TrueHaplotype(List.of(germline), 0.4), new TrueHaplotype(List.of(germline, somatic), 0.1)), p1 - 10, p2 + 10, 120});
        return cases.toArray(new Object[0][]);
    }

    @Test(dataProvider = "cases")
    public void testTrueHaplotypesAreAssembled(final String name, final String reference, final List<TrueHaplotype> truth,
                                               final int spanStart, final int spanEnd, final int depth) {
        final SAMFileHeader header = header(reference.length());

        // 0-based half-open span -> 1-based assembly region
        final AssemblyRegion region = new AssemblyRegion(new SimpleInterval(CONTIG, spanStart + 1, spanEnd), PADDING, header);
        simulate(header, reference, truth, depth, new Random(1)).stream()
                .filter(read -> region.getPaddedSpan().overlaps(read))
                .forEach(region::add);
        AssemblyBasedCallerUtils.finalizeRegion(region, false, false, (byte) 9, header, new IndexedSampleList(SAMPLE),
                false, false, false, false);

        final PositionalAssembler assembler = assembler();
        final byte[] window = reference.getBytes(StandardCharsets.US_ASCII);
        final int low = spanStart - PADDING, high = spanEnd + PADDING;
        final Set<List<PositionalEvent>> found = assembler.assemble(region.getReads(), window, 1, header).stream()
                .map(h -> h.events().stream().filter(e -> e.position() >= low && e.position() < high).collect(Collectors.toList()))
                .collect(Collectors.toSet());
        final Set<List<PositionalEvent>> expected = new HashSet<>();
        for (final TrueHaplotype haplotype : truth) {
            final List<PositionalEvent> normalized = new ArrayList<>();
            int floor = 0;
            for (final PositionalEvent event : haplotype.events.stream().sorted().toList()) {
                for (final PositionalEvent e : PositionalEvent.normalize(reference, event.position(), event.ref(), event.alt(), floor,
                        SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS)) {
                    normalized.add(e);
                    floor = e.position() + Math.max(1, e.ref().length());
                }
            }
            normalized.sort(Comparator.naturalOrder());
            expected.add(normalized);
        }
        for (final List<PositionalEvent> haplotype : expected) {
            Assert.assertTrue(found.contains(haplotype), name + ": true haplotype " + haplotype + " not among assembled " + found);
        }
    }

    /**
     * As in a de Bruijn graph, a mismatch in a read's first or last few bases (a dangling end shorter than
     * --min-dangling-branch-length kmers) can't make an event by itself, but it supports an event that other reads make.
     */
    @Test
    public void testShortDanglingEndsOnlySupportEvents() {
        final String reference = randomBases(new Random(3), 400);
        final int p = 200, length = 60;
        final PositionalEvent snv = new PositionalEvent(p, reference.substring(p, p + 1), reference.charAt(p) == 'A' ? "C" : "A");
        final SAMFileHeader header = header(reference.length());
        for (final boolean midRead : new boolean[] {false, true}) {
            final List<GATKRead> reads = new ArrayList<>();
            for (int i = 0; i < 30; i++) {
                reads.add(read(header, reference, p - 30 - i, length, null, "ref" + i));
            }
            // reads that start one base before the SNV or end one base after it: dangling ends of 2 kmers
            for (int i = 0; i < 5; i++) {
                reads.add(read(header, reference, p - 1, length, snv, "head" + i));
                reads.add(read(header, reference, p + 2 - length, length, snv, "tail" + i));
            }
            if (midRead) {
                reads.add(read(header, reference, p - length / 2, length, snv, "middle"));
            }
            reads.sort(Comparator.comparingInt(GATKRead::getStart).thenComparing(GATKRead::getName));
            final boolean assembled = assembler().assemble(reads, reference.getBytes(StandardCharsets.US_ASCII), 1, header).stream()
                    .anyMatch(h -> h.events().contains(snv));
            Assert.assertEquals(assembled, midRead);
        }
    }

    /**
     * Two sites: a germline het at site 1 and a recurrent error, allele 1 at site 2, in 2 fragments that reach only
     * site 2.  The error branches off the alt haplotype, whose fragments never reach site 2; that open site must take
     * the majority allele 0, or the germline haplotype (1, 0) would be lost.
     */
    @Test
    public void testOpenSitesTakeTheMajorityAllele() {
        final PositionalHaplotypeFinder finder = new PositionalHaplotypeFinder(PositionalAssemblyParameters.defaults(10, 2, true, 0.001,
                ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD, 10, MIN_DANGLING_BRANCH_LENGTH,
                SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS), new byte[0], List.of());
        final List<long[]> patterns = List.of(new long[] {1, 1}, new long[] {2, 0}, new long[] {0, 2});
        final PositionalHaplotypeFinder.Cover cover = finder.cover(2, List.of(List.of("A", "C"), List.of("G", "T")), patterns, new int[] {40, 30, 2});
        final List<List<Integer>> haplotypes = cover.alleles().stream().map(a -> Arrays.stream(a).boxed().toList()).toList();
        Assert.assertEquals(haplotypes, List.of(List.of(0, 0), List.of(1, 0), List.of(1, 1)));
        Assert.assertEquals(cover.supports(), new int[] {40, 30, 2});
    }

    @Test
    public void testAlignmentCost() {
        final byte[] haplotype = "ACGTTGCAAGTC".getBytes(StandardCharsets.US_ASCII);
        Assert.assertEquals(cost(haplotype, "ACGTTGCAAGTC", -1, false), 0);
        Assert.assertEquals(cost(haplotype, "ACGTAGCAAGTC", -1, false), 1);
        // a mismatch at a base below the minimum quality is free
        Assert.assertEquals(cost(haplotype, "ACGTAGCAAGTC", 4, false), 0);
        // a base missing from the read, or an extra one, is a gap; read bases past the haplotype are free
        Assert.assertEquals(cost(haplotype, "ACGTGCAAGTCA", -1, false), 2);
        Assert.assertEquals(cost(haplotype, "ACGTTTGCAAGTC", -1, false), 2);
        // haplotype bases past the read's end are free only if the read really ends there
        Assert.assertEquals(cost(haplotype, "ACGTTG", -1, true), 0);
        Assert.assertTrue(cost(haplotype, "ACGTTG", -1, false) > 2);
    }

    private static PositionalAssembler assembler() {
        return new PositionalAssembler(PositionalAssembler.DEFAULT_KMER_SIZE, 2, true, 0.001, ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD,
                10, MIN_DANGLING_BRANCH_LENGTH, 128, SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS, false);
    }

    private static SAMFileHeader header(final int referenceLength) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, referenceLength);
        final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(READ_GROUP);
        readGroup.setSample(SAMPLE);
        header.addReadGroup(readGroup);
        return header;
    }

    /** an unpaired read of the reference from 0-based offset start, with an SNV if one is given */
    private static GATKRead read(final SAMFileHeader header, final String reference, final int start, final int length,
                                 final PositionalEvent snv, final String name) {
        final byte[] bases = reference.substring(start, start + length).getBytes(StandardCharsets.US_ASCII);
        if (snv != null) {
            bases[snv.position() - start] = (byte) snv.alt().charAt(0);
        }
        final byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) 30);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, start + 1, bases, qualities, length + "M");
        read.setReadGroup(READ_GROUP);
        read.setMappingQuality(60);
        return read;
    }

    private static int cost(final byte[] haplotype, final String read, final int lowQualityBase, final boolean readEnded) {
        final byte[] qualities = new byte[read.length()];
        Arrays.fill(qualities, (byte) 30);
        if (lowQualityBase >= 0) {
            qualities[lowQualityBase] = 5;
        }
        return PositionalHaplotypeFinder.alignmentCost(haplotype, read.getBytes(StandardCharsets.US_ASCII), qualities, readEnded, 4, 10);
    }

    /** the haplotype's bases and, for each, its reference offset (-1 for inserted bases) */
    private static List<int[]> haplotype(final String reference, final List<PositionalEvent> events, final StringBuilder bases) {
        final List<int[]> positions = new ArrayList<>();
        int next = 0;
        for (final PositionalEvent event : events.stream().sorted().toList()) {
            for (int i = next; i < event.position(); i++) {
                bases.append(reference.charAt(i));
                positions.add(new int[] {i});
            }
            for (int i = 0; i < event.alt().length(); i++) {
                bases.append(event.alt().charAt(i));
                positions.add(new int[] {-1});
            }
            next = event.end();
        }
        for (int i = next; i < reference.length(); i++) {
            bases.append(reference.charAt(i));
            positions.add(new int[] {i});
        }
        return positions;
    }

    /**
     * The CIGAR operators and 0-based reference start of haplotype bases [a, b), as an aligner might report them:
     * inserted bases at either end, or any event within {@link #CLIP_IF_WITHIN} bases of an end, become soft clips.
     */
    private static String cigar(final List<int[]> positions, final int a, final int b, final int[] start) {
        final List<Character> ops = new ArrayList<>();
        int previous = -1;
        for (int i = a; i < b; i++) {
            final int p = positions.get(i)[0];
            if (p < 0) {
                ops.add('I');
            } else {
                if (previous >= 0 && p > previous + 1) {
                    for (int d = 0; d < p - previous - 1; d++) {
                        ops.add('D');
                    }
                }
                ops.add('M');
                previous = p;
            }
        }
        List<Character> clipped = clipNearEnd(ops);
        final List<Character> reversed = new ArrayList<>(clipped);
        java.util.Collections.reverse(reversed);
        clipped = clipNearEnd(reversed);
        java.util.Collections.reverse(clipped);
        while (!clipped.isEmpty() && (clipped.get(0) == 'I' || clipped.get(0) == 'D')) {
            if (clipped.get(0) == 'I') {
                clipped.set(0, 'S');
                break;
            }
            clipped.remove(0);
        }
        while (!clipped.isEmpty() && (clipped.get(clipped.size() - 1) == 'I' || clipped.get(clipped.size() - 1) == 'D')) {
            if (clipped.get(clipped.size() - 1) == 'I') {
                clipped.set(clipped.size() - 1, 'S');
                break;
            }
            clipped.remove(clipped.size() - 1);
        }
        // leading soft clips, including inserted bases next to them, then the first aligned base
        int leadingReadBases = 0;
        while (clipped.get(leadingReadBases) == 'S') {
            leadingReadBases++;
        }
        int index = a + leadingReadBases;
        while (positions.get(index)[0] < 0) {
            index++;
        }
        start[0] = positions.get(index)[0];
        final StringBuilder cigar = new StringBuilder();
        for (int i = 0; i < clipped.size(); ) {
            int j = i;
            while (j < clipped.size() && clipped.get(j).equals(clipped.get(i))) {
                j++;
            }
            cigar.append(j - i).append(clipped.get(i));
            i = j;
        }
        return cigar.toString();
    }

    /** soft-clip the start of the read up to the first aligned base after its first event, if that event is near the start */
    private static List<Character> clipNearEnd(final List<Character> ops) {
        int firstEvent = -1;
        for (int i = 0; i < ops.size(); i++) {
            if (ops.get(i) != 'M') {
                firstEvent = i;
                break;
            }
        }
        if (firstEvent < 0 || ops.get(firstEvent) == 'S') {
            return new ArrayList<>(ops);
        }
        final long readBasesBefore = ops.subList(0, firstEvent).stream().filter(op -> op != 'D').count();
        if (readBasesBefore >= CLIP_IF_WITHIN) {
            return new ArrayList<>(ops);
        }
        int k = firstEvent;
        while (k < ops.size() && ops.get(k) != 'M') {
            k++;
        }
        final long clippedReadBases = ops.subList(0, k).stream().filter(op -> op != 'D').count();
        final List<Character> result = new ArrayList<>();
        for (int i = 0; i < clippedReadBases; i++) {
            result.add('S');
        }
        result.addAll(ops.subList(k, ops.size()));
        return result;
    }

    private static List<GATKRead> simulate(final SAMFileHeader header, final String reference, final List<TrueHaplotype> haplotypes,
                                           final int depth, final Random random) {
        final List<GATKRead> reads = new ArrayList<>();
        final int fragments = (int) ((double) depth * reference.length() / (2 * READ_LENGTH));
        final double total = haplotypes.stream().mapToDouble(TrueHaplotype::weight).sum();
        final byte[] qualities = new byte[READ_LENGTH];
        Arrays.fill(qualities, (byte) 40);
        for (int f = 0; f < fragments; f++) {
            double x = random.nextDouble() * total;
            TrueHaplotype chosen = haplotypes.get(haplotypes.size() - 1);
            for (final TrueHaplotype haplotype : haplotypes) {
                if (x < haplotype.weight) {
                    chosen = haplotype;
                    break;
                }
                x -= haplotype.weight;
            }
            final StringBuilder sequence = new StringBuilder();
            final List<int[]> positions = haplotype(reference, chosen.events, sequence);
            final int s = random.nextInt(sequence.length() - FRAGMENT_LENGTH);
            final int[][] mates = {{s, s + READ_LENGTH}, {s + FRAGMENT_LENGTH - READ_LENGTH, s + FRAGMENT_LENGTH}};
            final GATKRead[] pair = new GATKRead[2];
            for (int mate = 0; mate < 2; mate++) {
                final char[] bases = sequence.substring(mates[mate][0], mates[mate][1]).toCharArray();
                for (int i = 0; i < bases.length; i++) {
                    if (random.nextDouble() < ERROR_RATE) {
                        final char original = bases[i];
                        do {
                            bases[i] = "ACGT".charAt(random.nextInt(4));
                        } while (bases[i] == original);
                    }
                }
                final int[] start = new int[1];
                final String cigar = cigar(positions, mates[mate][0], mates[mate][1], start);
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "fragment" + f, 0, start[0] + 1,
                        new String(bases).getBytes(StandardCharsets.US_ASCII), qualities.clone(), cigar);
                read.setReadGroup(READ_GROUP);
                read.setMappingQuality(60);
                read.setIsPaired(true);
                read.setIsProperlyPaired(true);
                read.setIsReverseStrand(mate == 1);
                read.setMateIsReverseStrand(mate == 0);
                if (mate == 0) {
                    read.setIsFirstOfPair();
                } else {
                    read.setIsSecondOfPair();
                }
                pair[mate] = read;
            }
            for (int mate = 0; mate < 2; mate++) {
                pair[mate].setMatePosition(CONTIG, pair[1 - mate].getStart());
                pair[mate].setFragmentLength(mate == 0 ? FRAGMENT_LENGTH : -FRAGMENT_LENGTH);
                reads.add(pair[mate]);
            }
        }
        reads.sort(Comparator.comparingInt(GATKRead::getStart).thenComparing(GATKRead::getName).thenComparing(GATKRead::isReverseStrand));
        return reads;
    }
}
