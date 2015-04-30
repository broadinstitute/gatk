package org.broadinstitute.hellbender.utils.pileup;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * Represents a pileup of reads at a given position.
 */
public final class ReadPileup implements Iterable<PileupElement>{
    private final Locatable loc;
    private final List<PileupElement> pileupElements;

    /**
     * Create a new version of a read backed pileup at loc, using the reads and their corresponding
     * offsets.
     * Note: This internal-use constructor keeps an alias pointer to the given list.
     */
    @VisibleForTesting
    ReadPileup(final Locatable loc, final List<PileupElement> pileup) {
        this.loc = loc;
        this.pileupElements = pileup;
    }

    /**
     * Create a new pileup without any aligned reads
     */
    public ReadPileup(final Locatable loc) {
        this(loc, Collections.emptyList());
    }

    /**
     * Create a new pileup with the given reads.
     */
    public ReadPileup(final Locatable loc, final List<SAMRecord> reads, final int offset) {
        this(loc, readsOffsets2Pileup(reads, offset));
    }


    /**
     * Create a new pileup with the given reads.
     */
    public ReadPileup(final Locatable loc, final List<SAMRecord> reads, final List<Integer> offsets) {
        this(loc, readsOffsets2Pileup(reads, offsets));
    }

    /**
     * Helper routine for converting reads and offset lists to a PileupElement list.
     */
    private static List<PileupElement> readsOffsets2Pileup(final List<SAMRecord> reads, final List<Integer> offsets) {
        if (reads == null) {
            throw new GATKException("Illegal null read list in UnifiedReadBackedPileup");
        }
        if (offsets == null) {
            throw new GATKException("Illegal null offsets list in UnifiedReadBackedPileup");
        }
        if (reads.size() != offsets.size()) {
            throw new GATKException("Reads and offset lists have different sizes!");
        }

        final List<PileupElement> pileup = new LinkedList<>();
        for (int i = 0; i < reads.size(); i++) {
            pileup.add(PileupElement.createPileupForReadAndOffset(reads.get(i), offsets.get(i)));
        }

        return pileup;
    }

    /**
     * Helper routine for converting reads and a single offset to a PileupElement list.
     */
    private static List<PileupElement> readsOffsets2Pileup(final List<SAMRecord> reads, final int offset) {
        if (reads == null) {
            throw new GATKException("Illegal null read list in UnifiedReadBackedPileup");
        }
        if (offset < 0) {
            throw new GATKException("Illegal offset < 0 UnifiedReadBackedPileup");
        }

        final List<PileupElement> pileup = new LinkedList<>();
        for (SAMRecord read : reads) {
            pileup.add(PileupElement.createPileupForReadAndOffset(read, offset));
        }

        return pileup;
    }

    /**
     * Make a new pileup consisting of elements of this pileup that satisfy the predicate.
     */
    public ReadPileup makeFilteredPileup(final Predicate<PileupElement> filter){
        return new ReadPileup(loc, pileupElements.stream().filter(filter).collect(Collectors.toList()));
    }

    /**
     * Make a new pileup consisting of only reads on the negative strand.
     */
    public ReadPileup getPositiveStrandPileup() {
        return makeFilteredPileup(p -> !p.getRead().getReadNegativeStrandFlag());
    }

    /**
     * Make a new pileup consisting of only reads on the negative strand.
     */
    public ReadPileup getNegativeStrandPileup() {
        return makeFilteredPileup(p -> p.getRead().getReadNegativeStrandFlag());
    }

    /**
     * Make a new pileup that contains of only bases with quality >= minBaseQ, coming from
     * reads with mapping qualities >= minMapQ.
     */
    public ReadPileup getBaseAndMappingFilteredPileup(int minBaseQ, int minMapQ) {
        return makeFilteredPileup(p -> p.getRead().getMappingQuality() >= minMapQ && (p.isDeletion() || p.getQual() >= minBaseQ));
    }

    /**
     * Make a new pileup that contains of only bases with quality >= minBaseQ.
     */
    public ReadPileup getBaseFilteredPileup(int minBaseQ) {
        return makeFilteredPileup(p -> (p.isDeletion() || p.getQual() >= minBaseQ));
    }

    /**
     * Make a new pileup that contains of only bases coming from reads with mapping quality >= minMapQ.
     */
    public ReadPileup getMappingFilteredPileup(int minMapQ) {
        return makeFilteredPileup(p -> p.getRead().getMappingQuality() >= minMapQ);
    }

    /**
     * Makes a new pileup by subsetting reads from a given read group.
     * Returns null if no reads are from the given read group.
     * @param targetReadGroupId Identifier for the read group.
     * @return A read-backed pileup containing only the reads in the given read group.
     */
    public ReadPileup getPileupForReadGroup(String targetReadGroupId) {
        final ReadPileup pu = makeFilteredPileup(p -> {
            SAMRecord read = p.getRead();
            if (targetReadGroupId != null) {
                if (read.getReadGroup() != null && targetReadGroupId.equals(read.getReadGroup().getReadGroupId())) {
                    return true;
                }
            } else {
                if (read.getReadGroup() == null || read.getReadGroup().getReadGroupId() == null) {
                    return true;
                }
            }
            return false;
        });
        return pu.isEmpty() ? null: pu;
    }

    /**
     * Make a new pileup from elements whose reads have read groups that agree with the given ID.
     * (it they have a name equal to the ID or starting with ID followed by a period ".").
     * Retrurns null if no suitable elements are found.
     */
    public ReadPileup getPileupForLane(String laneID) {
        List<PileupElement> filteredTracker = new LinkedList<>();
        for (PileupElement p : pileupElements) {
            SAMRecord read = p.getRead();
            final SAMReadGroupRecord readGroup = read.getReadGroup();
            if (laneID != null) {
                if (readGroup != null &&
                        (readGroup.getReadGroupId().startsWith(laneID + ".") ||   // lane is the same, but sample identifier is different
                         readGroup.getReadGroupId().equals(laneID)))               // in case there is no sample identifier, they have to be exactly the same
                {
                    filteredTracker.add(p);
                }
            } else {
                if (readGroup == null || readGroup.getReadGroupId() == null) {
                    filteredTracker.add(p);
                }
            }
        }
        return filteredTracker.size() > 0 ? new ReadPileup(loc, filteredTracker) : null;
    }

    /**
     * Gets a set of the read groups represented in this pileup.
     * Note: contains null if a read has a null read group.
     */
    public Set<SAMReadGroupRecord> getReadGroupIDs() {
        return pileupElements.stream().map(pe -> pe.getRead().getReadGroup()).collect(Collectors.toSet());
    }

    /**
     * Gets a set of the samples represented in this pileup.
     * Note: contains null if a read has a null read group or a null sample name.
     */
    public Set<String> getSamples() {
        return pileupElements.stream().map(pe -> pe.getRead()).map(r -> r.getReadGroup() == null ? null: r.getReadGroup().getSample()).collect(Collectors.toSet());
    }


    /**
     * The best way to access PileupElements where you only care about the bases and quals in the pileup.
     * <p/>
     * for (PileupElement p : this) { doSomething(p); }
     * <p/>
     * Provides efficient iteration of the data.
     */
    public Iterator<PileupElement> iterator() {
        return Collections.unmodifiableList(pileupElements).iterator();
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
     * @return the location of this pileup
     */
    public Locatable getLocation() {
        return loc;
    }

    /**
     * Get counts of A, C, G, T in order, which returns a int[4] vector with counts according
     * to BaseUtils.simpleBaseToBaseIndex for each base.
     */
    public int[] getBaseCounts() {
        final int[] counts = new int[4];

        for (final PileupElement pile : this) {
            // skip deletion sites
            if (!pile.isDeletion()) {
                int index = BaseUtils.simpleBaseToBaseIndex(pile.getBase());
                if (index != -1) {
                    counts[index]++;
                }
            }
        }

        return counts;
    }

    /**
     * Creates a pileup formar string.
     * In the pileup format, each line represents a genomic position, consisting of chromosome name,
     * coordinate, reference base, read bases, read qualities and alignment mapping qualities.
     */
    public String getPileupString(final Character ref) {
        return String.format("%s %s %c %s %s",
                getLocation().getContig(), getLocation().getStart(),    // chromosome name and coordinate
                ref,                                                     // reference base
                new String(getBases()),
                getQualsString());
    }

    /**
     * Returns a list of the reads in this pileup. Note this call costs O(n) and allocates fresh lists each time
     */
    public List<SAMRecord> getReads() {
        return pileupElements.stream().map(pe -> pe.getRead()).collect(Collectors.toList());
    }

    /**
     * The number of deletion bases in this pileup.
     */
    public int getNumberOfDeletions() {
        //Note: pileups are small so int is fine.
        return (int)pileupElements.stream().filter(p -> p.isDeletion()).count();
    }

    /**
     * The number of reads with mapping quality 0 in this pileup.
     */
    public int getNumberOfMappingQualityZeroReads() {
        //Note: pileups are small so int is fine.
        return (int)pileupElements.stream().filter(p -> p.getRead().getMappingQuality() == 0).count();
    }

    /**
     * The number of reads with deletions after this element.
     */
    public int getNumberOfDeletionsAfterThisElement() {
        //Note: pileups are small so int is fine.
        return (int)pileupElements.stream().filter(pe -> pe.isBeforeDeletionStart()).count();
    }

    /**
     * The number of reads with insertions after this element.
     */
    public int getNumberOfInsertionsAfterThisElement() {
        //Note: pileups are small so int is fine.
        return (int)pileupElements.stream().filter(pe -> pe.isBeforeInsertion()).count();
    }

    /**
     * Returns a list of the offsets in this pileup.
     * Note: this call costs O(n) and allocates fresh lists each time
     */
    public List<Integer> getOffsets() {
        return pileupElements.stream().map(pe -> pe.getOffset()).collect(Collectors.toList());
    }

    /**
     * Returns an array of the bases in this pileup.
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public byte[] getBases() {
        final byte[] v = new byte[size()];
        int pos = 0;
        for (final PileupElement pile : pileupElements) {
            v[pos++] = pile.getBase();
        }
        return v;
    }

    /**
     * Returns an array of the quals in this pileup.
     * Note: this call costs O(n) and allocates fresh array each time
     */
    public byte[] getQuals() {
        final byte[] v = new byte[size()];
        int pos = 0;
        for (final PileupElement pile : pileupElements) {
            v[pos++] = pile.getQual();
        }
        return v;
    }

    /**
     * Get an array of the mapping qualities.
     */
    public int[] getMappingQuals() {
        final int[] v = new int[size()];
        int pos = 0;
        for (final PileupElement pile : pileupElements) {
            v[pos++] = pile.getRead().getMappingQuality();
        }
        return v;
    }

    private static String quals2String(byte[] quals) {
        final StringBuilder qualStr = new StringBuilder();
        for (int qual : quals) {
            qual = Math.min(qual, 63);              // todo: fixme, this isn't a good idea
            char qualChar = (char) (33 + qual);     // todo: warning, this is illegal for qual > 63
            qualStr.append(qualChar);
        }

        return qualStr.toString();
    }

    private String getQualsString() {
        return quals2String(getQuals());
    }

    /**
     * Returns a new ReadBackedPileup that is sorted by start coordinate of the reads.
     */
    public ReadPileup getStartSortedPileup() {

        final TreeSet<PileupElement> sortedElements = new TreeSet<>(new Comparator<PileupElement>() {
            @Override
            public int compare(PileupElement element1, PileupElement element2) {
                final int difference = element1.getRead().getAlignmentStart() - element2.getRead().getAlignmentStart();
                return difference != 0 ? difference : element1.getRead().getReadName().compareTo(element2.getRead().getReadName());
            }
        });

        sortedElements.addAll(pileupElements);
        return new ReadPileup(loc, new LinkedList<>(sortedElements));
    }

    public FragmentCollection<PileupElement> toFragments() {
        return FragmentCollection.create(this);
    }
}
