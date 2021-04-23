package org.broadinstitute.hellbender.utils;

import com.google.inject.ImplementedBy;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.builder.EqualsExclude;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import javax.validation.OverridesAttribute;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public interface IntervalPileup {

    byte NO_BQ = (byte) -1;
    byte GAP = (byte) '-';
    byte NO_BASE = (byte) -1;

    List<GATKRead> reads();
    ReferenceBases reference();
    int width();
    int height();
    byte baseAt(final int row, final int column); // if overlapping and different calls, it returns the call for the "first" read.
    byte qualAt(final int row, final int column); // if overlapping and different quals, it return the qual for the "first" read.
    Insert insertAt(final int row, final int column);

    GATKRead readAt(final int row, final int column); // if overlapping, it returns the "first" read (possibly random).
    List<GATKRead> readsAt(final int row, final int column); // if mates are overlapping and we force mates on the same row we may have more than one read.

    Element element(final GATKRead read);
    default Element element(int row) {
        final GATKRead read = reads().get(row);
        return element(read);
    }

    boolean hasInsertAt(int i, int j);

    byte[] insertBasesAt(int i, int j);

    byte[] insertQualsAt(int i, int j);

    interface Element {
        GATKRead read();
        int row();
        int minColumn();
        int maxColumn();
        Insert insertAt(final int column);
        List<Insert> inserts();
        default List<Insert> inserts(final int firstColumn, final int lastColumn) {
            if (lastColumn < firstColumn) {
                return Collections.emptyList();
            }
            final List<Insert> all = inserts();
            if (all.isEmpty()) {
                return Collections.emptyList();
            } else if (all.size() == 1) {
                final int column = all.get(0).column();
                return column >= firstColumn && column <= lastColumn ? all : Collections.emptyList();
            } else {
                int i;
                for (i = 0; i < all.size(); i++) {
                    if (all.get(i).column() >= firstColumn) {
                        break;
                    }
                }
                int j, k;
                for (j = i, k = 0; j < all.size(); j++, k++) {
                    if (all.get(j).column() > lastColumn) {
                        break;
                    }
                }
                if (k == 0) {
                    return Collections.emptyList();
                } else if (k == 1) {
                    return Collections.singletonList(all.get(i));
                } else {
                    return all.subList(i, i + k);
                }
            }
        }

        boolean hasInsertAt(final int column);
        int insertSize(final int column);
        int copyInsertBases(int column, final byte[] dest, final int offset, final int length);
        int copyInsertQuals(int column, byte[] dest, int offset, int maxLength);
        byte[] insertQualsAt(int column);
        byte[] insertBasesAt(final int column);
        byte baseAt(final int column);
        byte qualAt(final int column);
    }

    interface Insert {
        int column();
        int length();
        byte[] bases();
        byte[] quals();
        /**
         * Returns true iff and only iff the {@code other} object is also an insert and
         * has exactly the same bases and qualities.
         * @param other
         * @return
         * @see #hashCode()
         */
        @Override
        boolean equals(final Object other);

        /**
         * Must be overrided in agreement with equals.
         * @return
         */
        @Override
        int hashCode();

        int copyBases(int offset, byte[] dest, int destOffset, final int maxLength);

        default int copyBases(int offset, byte[] dest, int destOffset) {
            return copyBases(offset, dest, destOffset, Integer.MAX_VALUE);
        }
        default int copyBases(byte[] dest) {
            return copyBases(0, dest, 0, Integer.MAX_VALUE);
        }

        int copyQuals(int offset, byte[] dest, int destOffset, final int maxLength);

        default int copyQuals(int offset, byte[] dest, int destOffset) {
            return copyQuals(offset, dest, destOffset, Integer.MAX_VALUE);
        }
        default int copyQuals(byte[] dest) {
            return copyQuals(0, dest, 0, Integer.MAX_VALUE);
        }
    }


    static IntervalPileup of(final Locatable loc, final ReadsDataSource aln, final ReferenceDataSource ref) {
        Utils.nonNull(loc);
        Utils.nonNull(aln);
        Utils.nonNull(ref);
        final String contig = loc.getContig();
        final int start = loc.getStart();
        final int end = loc.getEnd();
        Utils.nonNull(contig);
        Utils.validate(start <= end, "start must be less than end");
        final SAMSequenceRecord sequence = ref.getSequenceDictionary().getSequence(contig);
        if (sequence == null) {
            throw new IllegalArgumentException("unknown contig id " + contig);
        } else if (start <= 0) {
            throw new IllegalArgumentException("the start must be greater than 0");
        } else if (end > sequence.getSequenceLength()) {
            throw new IllegalArgumentException("the end must be less than the sequence length or equal to");
        }
        final SimpleInterval interval = new SimpleInterval(loc);
        final ReferenceBases referenceBases = new ReferenceBases(ref.queryAndPrefetch(interval).getBases(), interval);
        final List<GATKRead> reads = Utils.stream(aln.query(interval)).collect(Collectors.toList());
        return of(reads, referenceBases);
    }

    static IntervalPileup of(final List<GATKRead> reads, final ReferenceBases referenceBases) {
        Utils.nonNull(reads);
        Utils.nonNull(referenceBases);
        if (reads.isEmpty()) {
            return new EmptyIntervalPileup(referenceBases);
        } else {
            return new ByteMapIntervalPileup(referenceBases, reads);
        }
    }
}
