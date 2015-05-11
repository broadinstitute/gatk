package org.broadinstitute.hellbender.tools.picard.vcf.concordance;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.util.*;

/**
 * Takes a VCFFileReader and an IntervalList and provides a single iterator over all variants in all the intervals.
 *
 * @TODO Currently this uses the VCFFileReader.query method - could be useful to make a version of this iterator that uses the .iterator method
 *
 * @author Tim Fennell
 * @author George Grant
 */
public final class ByIntervalListVariantContextIterator implements Iterator<VariantContext> {
    private final VCFFileReader reader;
    private final Iterator<Interval> intervals;
    private CloseableIterator<VariantContext> currentIterator;

    public ByIntervalListVariantContextIterator(final VCFFileReader reader, final IntervalList intervals) {
        this.reader = reader;
        this.intervals = intervals.uniqued().iterator();
    }

    /** If the current iterator is null or exhausted, move to the next interval. */
    private void advance() {
        while ((currentIterator == null || !currentIterator.hasNext()) && this.intervals.hasNext()) {
            if (currentIterator != null) currentIterator.close();
            final Interval i = this.intervals.next();
            this.currentIterator = this.reader.query(i.getContig(), i.getStart(), i.getEnd());
        }
    }

    @Override
    public boolean hasNext() {
        advance();
        return this.currentIterator.hasNext();
    }

    @Override
    public VariantContext next() {
        advance();
        return this.currentIterator.next();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
