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
    private final SVReadFilter filter;
    private static final Iterator<QNameAndInterval> noName = Collections.emptyIterator();
    private int intervalsIndex = 0;

    public QNameFinder( final ReadMetadata metadata,
                        final List<SVInterval> intervals,
                        final SVReadFilter filter ) {
        this.metadata = metadata;
        this.intervals = intervals;
        this.filter = filter;
    }

    @Override
    public Iterator<QNameAndInterval> apply( final GATKRead read ) {
        if ( !filter.isMapped(read) ) return Collections.emptyIterator();

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
