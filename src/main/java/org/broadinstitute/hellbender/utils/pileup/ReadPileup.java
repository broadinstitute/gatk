package org.broadinstitute.hellbender.utils.pileup;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.qc.Pileup;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Represents a pileup of reads at a given position.
 */
public class ReadPileup implements Iterable<PileupElement> {
    private final Locatable loc;
    private final List<PileupElement> pileupElements;

    // Optional index of the elements by the samples of their reads, for pileups whose elements are grouped by sample, such
    // as those from LocusIteratorByState.  The elements of indexedSamples.get(n) (null for reads without a read group or
    // sample name) are those from sampleStarts[n] (inclusive) to sampleStarts[n + 1] (exclusive).  Both are null for a
    // pileup without an index, whose per-sample methods look up the sample of each read.
    private final SampleOrder indexedSamples;
    private final int[] sampleStarts;

    /** Constant used by samtools to downgrade a quality for overlapping reads that disagrees in their base. */
    public static final double SAMTOOLS_OVERLAP_LOW_CONFIDENCE = 0.8;

    /**
     * Create a new pileup at loc, using the reads and their corresponding
     * offsets.
     * Note: This constructor keeps an alias to the given list.
     */
    public ReadPileup(final Locatable loc, final List<PileupElement> pileup) {
        this.loc = loc;
        this.pileupElements = pileup;
        this.indexedSamples = null;
        this.sampleStarts = null;
    }

    /**
     * Create a new pileup at loc whose elements are grouped by the samples of their reads, so that splitting it by sample
     * takes no per-read work.  The elements of samples.get(n) are those from index sampleStarts[n] (inclusive) to
     * sampleStarts[n + 1] (exclusive).  The caller must ensure that every element's read belongs to its sample according to
     * {@link ReadUtils#getSampleName}, with a null sample for reads without a read group or sample name.
     * Note: This constructor keeps an alias to the given list and array.
     */
    public ReadPileup(final Locatable loc, final List<PileupElement> pileup, final SampleOrder samples, final int[] sampleStarts) {
        Utils.nonNull(pileup);
        Utils.nonNull(samples);
        Utils.nonNull(sampleStarts);
        Utils.validateArg(sampleStarts.length == samples.size() + 1, "there must be one more sample start than samples");
        Utils.validateArg(sampleStarts[0] == 0 && sampleStarts[samples.size()] == pileup.size(),
                "sample starts must begin at 0 and end at the size of the pileup");
        for (int n = 0; n < samples.size(); n++) {
            Utils.validateArg(sampleStarts[n] <= sampleStarts[n + 1], "sample starts must be non-decreasing");
        }
        this.loc = loc;
        this.pileupElements = pileup;
        this.indexedSamples = samples;
        this.sampleStarts = sampleStarts;
    }

    /**
     * Create a new pileup at loc, using an stratified pileup
     * Note: the current implementation of ReadPileup does not efficiently retrieve the stratified pileup
     */
    public ReadPileup(final Locatable loc, final Map<String, ReadPileup> stratifiedPileup) {
        // NOTE: this bit of code is a potential performance hotspot for the HaplotypeCaller.
        // This straightforward loop outperforms the equivalent streaming expression by over 2x.
        List<PileupElement> allElements = new ArrayList<>(stratifiedPileup.size() * 1000);
        for ( final Map.Entry<String, ReadPileup> pileupEntry : stratifiedPileup.entrySet() ) {
            allElements.addAll(pileupEntry.getValue().pileupElements);
        }

        this.loc = loc;
        this.pileupElements = allElements;
        this.indexedSamples = null;
        this.sampleStarts = null;
    }

    /**
     * Create a new pileup without any aligned reads
     */
    public ReadPileup(final Locatable loc) {
        this(loc, new ArrayList<PileupElement>());
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

    /**
     * Get the pileup of reads covering a locus.  This is useful, for example, in VariantWalkers, which work on
     * ReadsContexts and not AlignmentContexts.
     */
    public ReadPileup(final Locatable loc, final Iterable<GATKRead> reads) {
        final List<PileupElement> pile = StreamSupport.stream(reads.spliterator(), false)
                .filter(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK.and(ReadFilterLibrary.NOT_DUPLICATE))
                .map(AlignmentStateMachine::new)
                .map(asm -> {
                    while ( asm.stepForwardOnGenome() != null && asm.getGenomePosition() < loc.getStart()) { }
                    return asm.getGenomePosition() == loc.getStart() ? asm.makePileupElement() : null;
                }).filter(Objects::nonNull).collect(Collectors.toList());

        this.loc = loc;
        pileupElements = pile;
        indexedSamples = null;
        sampleStarts = null;
    }

    /**
     * Returns the first element corresponding to the given read or null there is no such element.
     *
     * @param read or null if elements with no reads are to be retrieved from this pileup.
     * @return the first element corresponding to the given read or null there is no such element.
     */
    @VisibleForTesting
    PileupElement getElementForRead(final GATKRead read) {
        return getElementStream().filter(el -> Objects.equals(el.getRead(), read)).findAny().orElse(null);
    }

    /**
     * Helper routine for converting reads and offset lists to a PileupElement list.
     */
    private static List<PileupElement> readsOffsetsToPileup(final List<GATKRead> reads, final List<Integer> offsets) {
        if (reads.size() != offsets.size()) {
            throw new IllegalArgumentException("Reads and offset lists have different sizes!");
        }

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
        return reads.stream().map(r -> PileupElement.createPileupForReadAndOffset(r, offset)).collect(Collectors.toList());
    }

    /**
     * Make a new pileup consisting of elements of this pileup that satisfy the predicate.
     * NOTE: the new pileup will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public ReadPileup makeFilteredPileup(final Predicate<PileupElement> filter) {
        return new ReadPileup(loc, getElementStream().filter(filter).collect(Collectors.toList()));
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
     * For a pileup indexed by sample, such as one from LocusIteratorByState, this takes no per-read work.
     * NOTE: the new pileup will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public ReadPileup getPileupForSample(final String sample, final SAMFileHeader header) {
        if (indexedSamples != null) {
            final int n = indexedSamples.indexOf(sample);
            return n == -1 ? new ReadPileup(loc) : getPileupForIndexedSample(n);
        }
        return makeFilteredPileup(pe -> Objects.equals(ReadUtils.getSampleName(pe.getRead(), header), sample));
    }

    /**
     * Make a new pileup from elements whose reads belong to samples that pass the filter, keeping their order.  A null
     * sample stands for reads without a read group or sample name, and the filter must depend only on the sample.
     * For a pileup indexed by sample, such as one from LocusIteratorByState, this takes no per-read work, and when every
     * element passes it returns this pileup itself.
     * NOTE: the new pileup will not be independent of the old one (no deep copy of the underlying data is performed).
     */
    public ReadPileup getPileupForSamples(final Predicate<String> sampleFilter, final SAMFileHeader header) {
        Utils.nonNull(sampleFilter);
        if (indexedSamples == null) {
            return makeFilteredPileup(pe -> sampleFilter.test(ReadUtils.getSampleName(pe.getRead(), header)));
        }

        final int nSamples = indexedSamples.size();
        final boolean[] kept = new boolean[nSamples];
        int nKeptSamples = 0;
        int nKeptElements = 0;
        int lastKeptSample = -1;
        for (int n = 0; n < nSamples; n++) {
            // a sample without elements has nothing to keep or drop
            if (sampleStarts[n] < sampleStarts[n + 1] && sampleFilter.test(indexedSamples.get(n))) {
                kept[n] = true;
                nKeptSamples++;
                nKeptElements += sampleStarts[n + 1] - sampleStarts[n];
                lastKeptSample = n;
            }
        }

        if (nKeptElements == size()) {
            return this;
        } else if (nKeptSamples == 0) {
            return new ReadPileup(loc);
        } else if (nKeptSamples == 1) {
            return getPileupForIndexedSample(lastKeptSample);
        }

        // copy the elements of the kept samples, in order, and index the copy by those samples
        final List<PileupElement> elements = new ArrayList<>(nKeptElements);
        final List<String> samples = new ArrayList<>(nKeptSamples);
        final int[] starts = new int[nKeptSamples + 1];
        for (int n = 0; n < nSamples; n++) {
            if (kept[n]) {
                starts[samples.size()] = elements.size();
                samples.add(indexedSamples.get(n));
                elements.addAll(pileupElements.subList(sampleStarts[n], sampleStarts[n + 1]));
            }
        }
        starts[nKeptSamples] = elements.size();
        return new ReadPileup(loc, elements, new SampleOrder(samples), starts);
    }

    // the elements of indexedSamples.get(n), as a view of this pileup's elements indexed by that sample
    private ReadPileup getPileupForIndexedSample(final int n) {
        final int start = sampleStarts[n];
        final int end = sampleStarts[n + 1];
        return new ReadPileup(loc, pileupElements.subList(start, end), SampleOrder.ofSingleSample(indexedSamples.get(n)), new int[] {0, end - start});
    }

    /**
     * Gets a set of the read groups represented in this pileup.
     * Note: contains null if a read has a null read group
     */
    public Set<String> getReadGroupIDs() {
        return getElementStream().map(pe -> pe.getRead().getReadGroup()).collect(Collectors.toSet());
    }

    /**
     * Gets a set of the samples represented in this pileup.
     * Note: contains null if a read has a null read group or a null sample name.
     */
    public Set<String> getSamples(final SAMFileHeader header) {
        if (indexedSamples != null) {
            // add the samples in order of first appearance, as below, so that the set iterates in the same order
            final Set<String> samples = new HashSet<>();
            for (int n = 0; n < indexedSamples.size(); n++) {
                if (sampleStarts[n] < sampleStarts[n + 1]) {
                    samples.add(indexedSamples.get(n));
                }
            }
            return samples;
        }
        return getElementStream().map(pe -> pe.getRead()).map(r -> ReadUtils.getSampleName(r, header)).collect(Collectors.toSet());
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
        final Map<String, ReadPileup> toReturn = new HashMap<>();
        for (final String sample : getSamples(header)) {
            final ReadPileup pileupBySample = getPileupForSample(sample, header);
            if (sample != null) {
                toReturn.put(sample, pileupBySample);
            } else {
                if (unknownSampleName == null) {
                    throw new UserException.ReadMissingReadGroup(pileupBySample.iterator().next().getRead());
                }
                toReturn.put(unknownSampleName, pileupBySample);
            }
        }
        return toReturn;
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
        // Profiling has determined that returning a custom unmodifiable iterator is faster than
        // Collections.unmodifiableList(pileupElements).iterator()
        return new Iterator<PileupElement>() {
            private final int len = pileupElements.size();
            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < len;
            }

            @Override
            public PileupElement next() {
                return pileupElements.get(i++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException("Cannot remove from a pileup element iterator");
            }
        };
    }

    /**
     * Iterator over sorted by read start PileupElements.
     */
    public Iterator<PileupElement> sortedIterator() {
        return getElementStream()
                .sorted((l, r) -> Integer.compare(l.getRead().getStart(), r.getRead().getStart()))
                .iterator();
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
        final FragmentCollection<PileupElement> fragments = FragmentCollection.create(this);
        fragments.getOverlappingPairs().stream()
                .forEach(
                        elements -> fixPairOverlappingQualities(elements.getLeft(), elements.getRight())
                );
    }

    /**
     * Fixes the quality of two elements that come from an overlapping pair in the same way as
     * samtools does {@see tweak_overlap_quality function in
     * <a href="https://github.com/samtools/htslib/blob/master/sam.c">samtools</a>}.
     * The only difference with the samtools API is the cap for high values ({@link QualityUtils#MAX_SAM_QUAL_SCORE}).
     */
    @VisibleForTesting
    static void fixPairOverlappingQualities(final PileupElement firstElement,
                                            final PileupElement secondElement) {
        // only if they do not represent deletions
        if (!secondElement.isDeletion() && !firstElement.isDeletion()) {
            final byte[] firstQuals = firstElement.getRead().getBaseQualities();
            final byte[] secondQuals = secondElement.getRead().getBaseQualities();
            if (firstElement.getBase() == secondElement.getBase()) {
                // if both have the same base, extra confidence in the firts of them
                firstQuals[firstElement.getOffset()] =
                        (byte) (firstQuals[firstElement.getOffset()] + secondQuals[secondElement
                                .getOffset()]);
                // cap to maximum byte value
                if (firstQuals[firstElement.getOffset()] < 0
                        || firstQuals[firstElement.getOffset()] > QualityUtils.MAX_SAM_QUAL_SCORE) {
                    firstQuals[firstElement.getOffset()] = QualityUtils.MAX_SAM_QUAL_SCORE;
                }
                secondQuals[secondElement.getOffset()] = 0;
            } else {
                // if not, we lost confidence in the one with higher quality
                if (firstElement.getQual() >= secondElement.getQual()) {
                    firstQuals[firstElement.getOffset()] =
                            (byte) (SAMTOOLS_OVERLAP_LOW_CONFIDENCE * firstQuals[firstElement.getOffset()]);
                    secondQuals[secondElement.getOffset()] = 0;
                } else {
                    secondQuals[secondElement.getOffset()] =
                            (byte) (SAMTOOLS_OVERLAP_LOW_CONFIDENCE * secondQuals[secondElement.getOffset()]);
                    firstQuals[firstElement.getOffset()] = 0;
                }
            }
            firstElement.getRead().setBaseQualities(firstQuals);
            secondElement.getRead().setBaseQualities(secondQuals);
        }
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
        return getElementStream().map(pe -> pe.getRead()).collect(Collectors.toList());
    }

    private Stream<PileupElement> getElementStream() {
        return pileupElements.stream();
    }

    /**
     * Returns the number of elements that satisfy the predicate.
     */
    public int getNumberOfElements(final Predicate<PileupElement> peFilter) {
        Utils.nonNull(peFilter);
        //Note: pileups are small so int is fine.
        return (int) getElementStream().filter(peFilter).count();
    }

    /**
     * Returns a list of the offsets in this pileup.
     * Note: this call costs O(n) and allocates fresh lists each time
     */
    public List<Integer> getOffsets() {
        return getElementStream().map(pe -> pe.getOffset()).collect(Collectors.toList());
    }

    //Extracts an int array by mapping each element in the pileup to int.
    private int[] extractIntArray(final ToIntFunction<PileupElement> map) {
        return getElementStream().mapToInt(map).toArray();
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

    /**
     * An immutable list of distinct samples with constant-time lookup of each sample's position, for indexing pileups by
     * sample.  A null sample stands for reads without a read group or sample name.  Build one for many pileups, as
     * LocusIteratorByState does, since construction takes time proportional to the number of samples.
     */
    public static final class SampleOrder {
        private final List<String> samples;
        private final Map<String, Integer> positions;

        public SampleOrder(final List<String> samples) {
            Utils.nonNull(samples);
            final Map<String, Integer> positions = new HashMap<>();
            for (int n = 0; n < samples.size(); n++) {
                final String sample = samples.get(n);
                Utils.validateArg(positions.put(sample, n) == null, () -> "sample " + sample + " appears more than once");
            }
            this.samples = Collections.unmodifiableList(new ArrayList<>(samples));
            this.positions = Collections.unmodifiableMap(positions);
        }

        private SampleOrder(final String sample) {
            this.samples = Collections.singletonList(sample);
            this.positions = Collections.singletonMap(sample, 0);
        }

        private static SampleOrder ofSingleSample(final String sample) {
            return new SampleOrder(sample);
        }

        public int size() {
            return samples.size();
        }

        public String get(final int n) {
            return samples.get(n);
        }

        /**
         * @return the position of the sample, or -1 if it isn't one of these samples
         */
        public int indexOf(final String sample) {
            final Integer position = positions.get(sample);
            return position == null ? -1 : position;
        }
    }
}
