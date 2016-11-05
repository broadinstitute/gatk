package org.broadinstitute.hellbender.utils.pileup;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;

/**
 * Represents a pileup of reads at a given position.
 */
public class ReadPileup implements Iterable<PileupElement> {
    private final Locatable loc;
    private PileupElementTracker pileupElements;

    /**
     * Create a new pileup at loc, using the reads and their corresponding
     * offsets.
     * Note: This constructor keeps an alias to the given list.
     */
    public ReadPileup(final Locatable loc, final List<PileupElement> pileup) {
        this(loc, pileup, false);
    }

    /**
     * Create a new pileup at loc, using the reads and their corresponding
     * offsets.
     * Note: This constructor keeps an alias to the given list.
     */
    public ReadPileup(final Locatable loc, final List<PileupElement> pileup, final boolean preSorted) {
        Utils.nonNull(loc, "loc is null");
        Utils.nonNull(pileup, "element list is null");
        this.loc = loc;
        this.pileupElements = new UnifiedPileupElementTracker(pileup, preSorted);
    }

    /**
     * Create a new pileup at loc, using an stratified pileup
     */
    public ReadPileup(final Locatable loc, final Map<String, ReadPileup> stratifiedPileup) {
        this.loc = loc;
        this.pileupElements = new PerSamplePileupElementTracker(stratifiedPileup.entrySet().stream()
                .collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue().pileupElements)));
    }

    /**
     * Create a new pileup without any aligned reads
     */
    public ReadPileup(final Locatable loc) {
        this(loc, Collections.emptyList(), true);
    }

    /**
     * Create a new pileup with the given reads.
     */
    public ReadPileup(final Locatable loc, final List<GATKRead> reads, final int offset) {
        this(loc, readsOffsetsToPileup(reads, offset));
    }

    /**
     * Create a new pileup with the given reads.
     */
    public ReadPileup(final Locatable loc, final List<GATKRead> reads, final List<Integer> offsets) {
        this(loc, readsOffsetsToPileup(reads, offsets));
    }

    // helper constructor
    private ReadPileup(final Locatable loc, final PileupElementTracker tracker) {
        this.loc = loc;
        this.pileupElements = tracker;
    }

    /**
     * Returns the first element corresponding to the given read or null there is no such element.
     *
     * @param read or null if elements with no reads are to be retrieved from this pileup.
     * @return the first element corresponding to the given read or null there is no such element.
     */
    @VisibleForTesting
    PileupElement getElementForRead(final GATKRead read) {
        return pileupElements.getElementStream().filter(el -> Objects.equals(el.getRead(), read)).findAny().orElse(null);
    }

    /**
     * Helper routine for converting reads and offset lists to a PileupElement list.
     */
    private static List<PileupElement> readsOffsetsToPileup(final List<GATKRead> reads, final List<Integer> offsets) {
        Utils.nonNull(reads, "Illegal null read list in UnifiedReadBackedPileup");
        Utils.nonNull(offsets, "Illegal null offsets list in UnifiedReadBackedPileup");
        Utils.validateArg(reads.size() == offsets.size(), "Reads and offset lists have different sizes!");

        final List<PileupElement> pileup = new ArrayList<>(reads.size());
        for (int i = 0; i < reads.size(); i++) {
            pileup.add(PileupElement.createPileupForReadAndOffset(reads.get(i), offsets.get(i)));
        }

        return pileup;
    }

    /**
     * Helper routine for converting reads and a single offset to a PileupElement list.
     */
    private static List<PileupElement> readsOffsetsToPileup(final List<GATKRead> reads, final int offset) {
        Utils.nonNull(reads, "Illegal null read list");
        return reads.stream().map(r -> PileupElement.createPileupForReadAndOffset(r, offset)).collect(Collectors.toList());
    }

    /**
     * Make a new pileup consisting of elements of this pileup that satisfy the predicate.
     * NOTE: the new pileup will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public ReadPileup makeFilteredPileup(final Predicate<PileupElement> filter) {
        Utils.nonNull(filter);
        return new ReadPileup(loc, pileupElements.makeFilteredTracker(filter));
    }

    /**
     * Make a new pileup from elements whose reads have read groups that agree with the given lane ID.
     * (if they have a name equal to the ID or starting with ID followed by a period ".").
     * Also, if both laneID and read group ID are {@code null}, the read will be included.
     * Returns empty pileup if no suitable elements are found.
     * NOTE: the new pileup will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public ReadPileup getPileupForLane(final String laneID) {
        return makeFilteredPileup(p -> {
            final GATKRead read = p.getRead();
            final String readGroupID = read.getReadGroup();
            if (laneID == null && readGroupID == null) {
                return true;
            }
            if (laneID != null && readGroupID != null) {
                final boolean laneSame = readGroupID.startsWith(laneID + "."); // lane is the same, but sample identifier is different
                final boolean exactlySame = readGroupID.equals(laneID);        // in case there is no sample identifier, they have to be exactly the same
                if (laneSame || exactlySame) {
                    return true;
                }
            }
            return false;
        });
    }

    /**
     * Make a new pileup from elements whose reads belong to the given sample.
     * Passing null sample as an argument retrieves reads whose read group or sample name is {@code null}
     * NOTE: the new pileup will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public ReadPileup getPileupForSample(final String sample, final SAMFileHeader header) {
        return new ReadPileup(loc, pileupElements.getTrackerForSample(sample, header));
    }

    /**
     * Gets a set of the read groups represented in this pileup.
     * Note: contains null if a read has a null read group
     */
    public Set<String> getReadGroupIDs() {
        return pileupElements.getElementStream().map(pe -> pe.getRead().getReadGroup()).collect(Collectors.toSet());
    }

    /**
     * Gets a set of the samples represented in this pileup.
     * Note: contains null if a read has a null read group or a null sample name.
     */
    public Set<String> getSamples(final SAMFileHeader header) {
        return pileupElements.getSamples(header);
    }

    /**
     * Splits the ReadPileup by sample
     *
     * @param header            the header to retrieve the samples from
     * @param unknownSampleName the sample name if there is no read group/sample name; {@code null} if lack of RG is not expected
     * @return a Map of sample name to ReadPileup (including empty pileups for a sample)
     * @throws org.broadinstitute.hellbender.exceptions.UserException.ReadMissingReadGroup if unknownSampleName is {@code null} and there are reads without RG/sample name
     */
    public Map<String, ReadPileup> splitBySample(final SAMFileHeader header, final String unknownSampleName) {
        // keep the splitting
        pileupElements = pileupElements.splitBySample(header);
        final Set<String> samples = getSamples(header);
        if (unknownSampleName == null && samples.contains(null)) {
            throw new UserException.ReadMissingReadGroup(getPileupForSample(null, header).iterator().next().getRead());
        }
        return getSamples(header).stream().collect(Collectors.toMap(
                s -> (s == null) ? unknownSampleName : s,
                s -> getPileupForSample(s, header)
        ));
    }

    /**
     * The best way to access PileupElements where you only care about the bases and quals in the pileup.
     * <p>
     * for (PileupElement p : this) { doSomething(p); }
     * <p>
     * Provides efficient iteration of the data.
     */
    @Override
    public Iterator<PileupElement> iterator() {
        return pileupElements.iterator();
    }

    /**
     * Iterator over sorted by read start PileupElements.
     */
    public Iterator<PileupElement> sortedIterator() {
        return pileupElements.sortedIterator();
    }

    /**
     * The number of elements in this pileup.
     */
    public int size() {
        return pileupElements.size();
    }

    /**
     * @return true if there are 0 elements in the pileup, false otherwise
     */
    public boolean isEmpty() {
        return size() == 0;
    }

    /**
     * @return the location of this pileup.
     */
    public Locatable getLocation() {
        return loc;
    }

    /**
     * Get counts of A, C, G, T in order, which returns a int[4] vector with counts according
     * to BaseUtils.simpleBaseToBaseIndex for each base.
     * Deletions are not counted.
     */
    public int[] getBaseCounts() {
        final int[] counts = new int[4];

        for (final PileupElement pile : this) {
            // skip deletion sites
            if (!pile.isDeletion()) {
                final int index = BaseUtils.simpleBaseToBaseIndex(pile.getBase());
                if (index != -1) {
                    counts[index]++;
                }
            }
        }

        return counts;
    }

    /**
     * Fixes the quality of all the elements that come from an overlapping pair in the same way as
     * samtools does {@see tweak_overlap_quality function in
     * <a href="https://github.com/samtools/htslib/blob/master/sam.c">samtools</a>}.
     * <p>
     * Setting the quality of one of the bases to 0 effectively removes the redundant base for
     * calling. In addition, if the bases overlap we have increased confidence if they agree (or
     * reduced if they don't). Thus, the algorithm proceeds as following:
     * <p>
     * 1. If the bases are the same, the quality of the first element is the sum of both qualities
     * and the quality of the second is reduced to 0.
     * 2. If the bases are different, the base with the highest quality is reduced with a factor of
     * 0.8, and the quality of the lowest is reduced to 0.
     * <p>
     * Note: Resulting qualities higher than {@link QualityUtils#MAX_SAM_QUAL_SCORE} are capped.
     */
    public void fixOverlaps() {
        pileupElements.fixOverlaps();
    }

    public final static Comparator<PileupElement> baseQualTieBreaker = new Comparator<PileupElement>() {
        @Override
        public int compare(PileupElement o1, PileupElement o2) {
            return Byte.compare(o1.getQual(), o2.getQual());
        }
    };

    public final static Comparator<PileupElement> mapQualTieBreaker = new Comparator<PileupElement>() {
        @Override
        public int compare(PileupElement o1, PileupElement o2) {
            return Integer.compare(o1.getMappingQual(), o2.getMappingQual());
        }
    };

    /**
     * Returns a new ReadPileup where only one read from an overlapping read
     * pair is retained.  If the two reads in question disagree to their basecall,
     * neither read is retained.  If they agree on the base, the read with the higher
     * base quality observation is retained
     *
     * @return the newly filtered pileup
     */
    public ReadPileup getOverlappingFragmentFilteredPileup(SAMFileHeader header) {
        return getOverlappingFragmentFilteredPileup(true, baseQualTieBreaker, header);
    }

    /**
     * Returns a new ReadPileup where only one read from an overlapping read
     * pair is retained.  If discardDiscordant and the two reads in question disagree to their basecall,
     * neither read is retained.  Otherwise, the read with the higher
     * quality (base or mapping, depending on baseQualNotMapQual) observation is retained
     *
     * @return the newly filtered pileup
     */
    public ReadPileup getOverlappingFragmentFilteredPileup(boolean discardDiscordant, Comparator<PileupElement> tieBreaker, SAMFileHeader header) {
        List<PileupElement> filteredPileupList = new ArrayList<PileupElement>();

        for (ReadPileup pileup : this.splitBySample(header, null).values()) {
            Collection<PileupElement> elements = filterSingleSampleForOverlaps(pileup, tieBreaker, discardDiscordant);
            filteredPileupList.addAll(elements);
        }

        return new ReadPileup(loc, filteredPileupList);
    }

    private Collection<PileupElement> filterSingleSampleForOverlaps(ReadPileup pileup, Comparator<PileupElement> tieBreaker, boolean discardDiscordant) {
        Map<String, PileupElement> filteredPileup = new HashMap<String, PileupElement>();
        Set<String> readNamesDeleted = new HashSet<>();

        for (PileupElement p : pileup) {
            String readName = p.getRead().getName();

            // if we've never seen this read before, life is good
            if (!filteredPileup.containsKey(readName)) {
                if(!readNamesDeleted.contains(readName)) {
                    filteredPileup.put(readName, p);
                }
            } else {
                PileupElement existing = filteredPileup.get(readName);

                // if the reads disagree at this position, throw them all out.  Otherwise
                // keep the element with the highest quality score
                if (discardDiscordant && existing.getBase() != p.getBase()) {
                    filteredPileup.remove(readName);
                    readNamesDeleted.add(readName);
                } else if (tieBreaker.compare(existing, p) < 0) {
                    filteredPileup.put(readName, p);
                }
            }
        }
        return(filteredPileup.values());
    }

    @Override
    public String toString() {
        return String.format("%s %s %s %s",
                loc.getContig(),
                loc.getStart(),
                new String(getBases()),
                getQualsString());
    }

    /**
     * Format, assuming a single-sample, in a samtools-like string.
     * Each line represents a genomic position, consisting of chromosome name, coordinate,
     * reference base, read bases and read qualities
     *
     * @param ref the reference base
     * @return pileup line
     */
    public String getPileupString(final char ref) {
        // In the pileup format,
        return String.format("%s %s %c %s %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                ref,                                                     // reference base
                new String(getBases()),
                getQualsString());
    }

    /**
     * Returns a list of the reads in this pileup. Note this call costs O(n) and allocates fresh lists each time
     */
    public List<GATKRead> getReads() {
        return pileupElements.getElementStream().map(pe -> pe.getRead()).collect(Collectors.toList());
    }

    /**
     * Returns the number of elements that satisfy the predicate.
     */
    public int getNumberOfElements(final Predicate<PileupElement> peFilter) {
        Utils.nonNull(peFilter);
        //Note: pileups are small so int is fine.
        return (int) pileupElements.getElementStream().filter(peFilter).count();
    }

    /**
     * Returns a list of the offsets in this pileup.
     * Note: this call costs O(n) and allocates fresh lists each time
     */
    public List<Integer> getOffsets() {
        return pileupElements.getElementStream().map(pe -> pe.getOffset()).collect(Collectors.toList());
    }

    //Extracts an int array by mapping each element in the pileup to int.
    private int[] extractIntArray(final ToIntFunction<PileupElement> map) {
        return pileupElements.getElementStream().mapToInt(map).toArray();
    }

    /**
     * Returns an array of the bases in this pileup.
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public byte[] getBases() {
        return toByteArray(extractIntArray(pe -> pe.getBase()));
    }

    /**
     * Returns an array of the quals in this pileup.
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public byte[] getBaseQuals() {
        return toByteArray(extractIntArray(pe -> pe.getQual()));
    }

    //Converts array of ints to array of bytes by hard casting (loses precision if ints are large).
    private byte[] toByteArray(final int[] ints) {
        final byte[] bytes = new byte[ints.length];
        for (int i = 0; i < ints.length; i++) {
            bytes[i] = (byte) ints[i];
        }
        return bytes;
    }

    /**
     * Get an array of the mapping qualities.
     */
    public int[] getMappingQuals() {
        return extractIntArray(pe -> pe.getMappingQual());
    }

    private String getQualsString() {
        final byte[] quals = getBaseQuals();
        for (int i = 0; i < quals.length; i++) {
            quals[i] = (byte) (33 + quals[i]);  //as per SAM spec
        }
        return new String(quals);
    }
}
