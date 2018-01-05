package org.broadinstitute.hellbender.utils.iterators;


import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Iterator;

/**
 * Wraps a SAMRecord iterator within an iterator of GATKReads.
 */
public final class SAMRecordToReadIterator implements Iterator<GATKRead>, Iterable<GATKRead> {
    private final Iterator<SAMRecord> samIterator;

    public SAMRecordToReadIterator( final Iterator<SAMRecord> samIterator ) {
        this.samIterator = samIterator;
    }

    @Override
    public boolean hasNext() {
        return samIterator.hasNext();
    }

    @Override
    public GATKRead next() {
        return new SAMRecordToGATKReadAdapter(samIterator.next());
    }

    @Override
    public Iterator<GATKRead> iterator() {
        return this;
    }
}