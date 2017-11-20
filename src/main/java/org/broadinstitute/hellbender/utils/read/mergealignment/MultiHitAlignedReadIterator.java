package org.broadinstitute.hellbender.utils.read.mergealignment;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Iterate over queryname-sorted SAM, and return each group of reads with the same queryname.  Unmapped reads
 * are filtered out, as are alignments that don't seem to match any part of the reference.
 * If there are multiple hits for the same read, and the first and second ends need to be correlated,
 * then they are sorted by hit index. Supplemental alignments are discarded, with a logged message.
 * A set of hits for a single query may then be filtered with a caller-supplied filter, which will remove any
 * alignments that do not pass the filter.  If the primary alignment is removed, the best-mapping secondary alignment
 * or alignment pair will be marked as primary.
 *
 */
public final class MultiHitAlignedReadIterator implements CloseableIterator<HitsForInsert> {
    private final PeekableIterator<SAMRecord> peekIterator;
    private final SAMRecordQueryNameComparator queryNameComparator = new SAMRecordQueryNameComparator();
    private final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy;

    private HitsForInsert theNext = null;

    /**
     *
     * @param querynameOrderIterator
     * @param primaryAlignmentSelectionStrategy Algorithm for selecting primary alignment when it is not clear from
     *                                          the input what should be primary.
     */
    MultiHitAlignedReadIterator(final CloseableIterator<SAMRecord> querynameOrderIterator,
                                final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy) {
        this.primaryAlignmentSelectionStrategy = primaryAlignmentSelectionStrategy;
        peekIterator = new PeekableIterator<>(new FilteringSamIterator(querynameOrderIterator,
                new SamRecordFilter() {
                    // Filter unmapped reads.
                    @Override
                    public boolean filterOut(final SAMRecord record) {
                        return record.getReadUnmappedFlag() || SAMUtils.cigarMapsNoBasesToRef(record.getCigar());
                    }
                    @Override
                    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
                        return ((first.getReadUnmappedFlag() || SAMUtils.cigarMapsNoBasesToRef(first.getCigar()))
                                && (second.getReadUnmappedFlag() || SAMUtils.cigarMapsNoBasesToRef(second.getCigar())));
                    }
                }));


        advance();
    }

    @Override
    public void close() {
        peekIterator.close();
    }

    @Override
    public boolean hasNext() {
        return theNext != null;
    }

    /**
     * @throws IllegalStateException if the input is not queryname-sorted.
     */
    @Override
    public HitsForInsert next() {
        if (!hasNext()) throw new NoSuchElementException();
        final HitsForInsert ret = theNext;
        advance();
        return ret;
    }

    private void advance() {
        while (peekIterator.hasNext()) {
            theNext = nextMaybeEmpty();
            if (theNext.numHits() > 0) return;
        }
        theNext = null;
    }

    private HitsForInsert nextMaybeEmpty() {
        Utils.validate(peekIterator.hasNext(), "iterator has no next");
        final String readName = peekIterator.peek().getReadName();
        final HitsForInsert hits = new HitsForInsert();

        Boolean isPaired = null;

        // Accumulate the alignments matching readName.
        do {
            final SAMRecord rec = peekIterator.next();
            replaceHardWithSoftClips(rec);
            // It is critical to do this here, because SamAlignmentMerger uses this exception to determine
            // if the aligned input needs to be sorted.
            if (peekIterator.hasNext() && queryNameComparator.fileOrderCompare(rec, peekIterator.peek()) > 0) {
                throw new IllegalStateException("Underlying iterator is not queryname sorted: " +
                        rec + " > " + peekIterator.peek());
            }

            if (isPaired == null) {
                isPaired = rec.getReadPairedFlag();
            } else if (isPaired != rec.getReadPairedFlag()) {
                throw new GATKException("Got a mix of paired and unpaired alignments for read " + readName);
            }

            // Records w/ a supplemental flag are stashed to the side until the primary alignment has
            // been determined, and then re-added into the process later
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                if (rec.getSupplementaryAlignmentFlag()) {
                    hits.addSupplementalFirstOfPairOrFragment(rec);
                } else {
                    hits.addFirstOfPairOrFragment(rec);
                }
            } else if (rec.getSecondOfPairFlag()) {
                if (rec.getSupplementaryAlignmentFlag()) {
                    hits.addSupplementalSecondOfPair(rec);
                } else {
                    hits.addSecondOfPair(rec);
                }
            } else throw new GATKException("Read is marked as pair but neither first or second: " + readName);
        } while (peekIterator.hasNext() && peekIterator.peek().getReadName().equals(readName));

        // If there is no more than one alignment for each end, no need to do any coordination.
        if (hits.numHits() <= 1) {
            // No HI tags needed if only a single hit
            if (hits.getFirstOfPair(0) != null) {
                hits.getFirstOfPair(0).setAttribute(SAMTag.HI.name(), null);
                hits.getFirstOfPair(0).setSecondaryAlignment(false);
            }
            if (hits.getSecondOfPair(0) != null) {
                hits.getSecondOfPair(0).setAttribute(SAMTag.HI.name(), null);
                hits.getSecondOfPair(0).setSecondaryAlignment(false);
            }
        } else {
            primaryAlignmentSelectionStrategy.pickPrimaryAlignment(hits);
        }

        // Used to check that alignments for first and second were correlated, but this is no longer required.
        return hits;
    }

    /** Replaces hard clips with soft clips and fills in bases and qualities with dummy values as needed. */
    private void replaceHardWithSoftClips(final SAMRecord rec) {
        if (rec.getReadUnmappedFlag()) return;
        if (rec.getCigar().isEmpty()) return;

        List<CigarElement> elements = rec.getCigar().getCigarElements();
        final CigarElement first = elements.get(0);
        final CigarElement last  = elements.size() == 1 ? null : elements.get(elements.size()-1);
        final int startHardClip = first.getOperator() == CigarOperator.H ? first.getLength() : 0;
        final int endHardClip   = (last != null && last.getOperator() == CigarOperator.H) ? last.getLength() : 0;

        if (startHardClip + endHardClip > 0) {
            final int len = rec.getReadBases().length + startHardClip + endHardClip;

            // Fix the basecalls
            final byte[] bases = new byte[len];
            Arrays.fill(bases, (byte) 'N');
            System.arraycopy(rec.getReadBases(), 0, bases, startHardClip, rec.getReadBases().length);

            // Fix the quality scores
            final byte[] quals = new byte[len];
            Arrays.fill(quals, (byte) 2  );
            System.arraycopy(rec.getBaseQualities(), 0, quals, startHardClip, rec.getBaseQualities().length);

            // Fix the cigar!
            elements = new ArrayList<>(elements); // make it modifiable
            if (startHardClip > 0) elements.set(0, new CigarElement(first.getLength(), CigarOperator.S));
            if (endHardClip   > 0) elements.set(elements.size()-1, new CigarElement(last.getLength(), CigarOperator.S));

            // Set the update structures on the new record
            rec.setReadBases(bases);
            rec.setBaseQualities(quals);
            rec.setCigar(new Cigar(elements));
        }
    }

    /** Unsupported operation. */
    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
