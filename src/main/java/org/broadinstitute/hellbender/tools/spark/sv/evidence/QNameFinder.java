package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.apache.commons.collections4.iterators.SingletonIterator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

/**
 * Class to find the template names associated with reads in specified intervals.
 */
public final class QNameFinder implements Function<GATKRead, Iterator<QNameAndInterval>> {
    private final ReadMetadata metadata;
    private final List<SVInterval> intervals;
    private static final Iterator<QNameAndInterval> noName = Collections.emptyIterator();
    private int intervalsIndex = 0;

    public QNameFinder( final ReadMetadata metadata,
                        final List<SVInterval> intervals ) {
        this.metadata = metadata;
        this.intervals = intervals;
    }

    @Override
    public Iterator<QNameAndInterval> apply( final GATKRead read ) {
        final int readContigId = metadata.getContigID(read.getContig());
        final int readStart = read.getUnclippedStart();
        final int intervalsSize = intervals.size();
        while ( intervalsIndex < intervalsSize ) {
            final SVInterval interval = intervals.get(intervalsIndex);
            if ( interval.getContig() > readContigId ) break;
            if ( interval.getContig() == readContigId && interval.getEnd() > read.getStart() ) break;
            intervalsIndex += 1;
        }
        if ( intervalsIndex >= intervalsSize ) return noName;
        final SVInterval indexedInterval = intervals.get(intervalsIndex);
        final SVInterval readInterval = new SVInterval(readContigId, readStart, read.getUnclippedEnd()+1);
        if ( indexedInterval.isDisjointFrom(readInterval) ) return noName;
        return new SingletonIterator<>(new QNameAndInterval(read.getName(), intervalsIndex));
    }
}
