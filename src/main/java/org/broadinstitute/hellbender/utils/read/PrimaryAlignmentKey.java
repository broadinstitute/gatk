package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMRecord;

/**
 * It is useful to define a key such that the key will occur at most once among the primary alignments in a given file
 * (assuming the file is valid). The read name + pairing status should be sufficient for this.
 */
public final class PrimaryAlignmentKey implements Comparable<PrimaryAlignmentKey> {

    private enum PairStatus {UNPAIRED, FIRST, SECOND} // note the order here; it should correspond to that of SAMRecordQueryNameComparator!

    private final PairStatus pairStatus;
    private final String readName;

    public PrimaryAlignmentKey(final SAMRecord rec) {
        this.pairStatus = rec.getReadPairedFlag() ?
                (rec.getSecondOfPairFlag() ? PairStatus.SECOND : PairStatus.FIRST) :
                PairStatus.UNPAIRED;
        this.readName = rec.getReadName();
    }

    @Override
    public int hashCode() {
        return 31 * this.readName.hashCode() + this.pairStatus.hashCode();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof PrimaryAlignmentKey)) return false;
        final PrimaryAlignmentKey that = (PrimaryAlignmentKey) o;
        return (this.compareTo(that) == 0);
    }


    @Override
    public int compareTo(final PrimaryAlignmentKey that) {
        int comp = this.readName.compareTo(that.readName);
        if (comp == 0) comp = this.pairStatus.compareTo(that.pairStatus);
        return comp;
    }
}
