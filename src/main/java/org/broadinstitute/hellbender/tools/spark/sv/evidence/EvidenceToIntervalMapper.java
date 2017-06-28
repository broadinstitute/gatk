package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.apache.commons.collections4.iterators.SingletonIterator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;

import java.util.Collections;
import java.util.Iterator;
import java.util.function.Function;

/**
 * A class to examine a stream of BreakpointEvidence, and group it into Intervals.
 */
public final class EvidenceToIntervalMapper implements Function<BreakpointEvidence, Iterator<SVInterval>> {
    private final int gapSize;
    private SVInterval curInterval;
    private static final Iterator<SVInterval> EMPTY_ITERATOR = Collections.emptyIterator();

    public EvidenceToIntervalMapper( final int gapSize ) {
        this.gapSize = gapSize;
    }

    @Override
    public Iterator<SVInterval> apply( final BreakpointEvidence evidence ) {
        Iterator<SVInterval> result = EMPTY_ITERATOR;
        final SVInterval interval = evidence.getLocation();
        if ( curInterval == null ) {
            curInterval = interval;
        } else if ( curInterval.gapLen(interval) < gapSize ) {
            curInterval = curInterval.join(interval);
        } else {
            result = new SingletonIterator<>(curInterval);
            curInterval = interval;
        }
        return result;
    }
}
