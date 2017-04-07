package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.commons.collections4.iterators.SingletonIterator;

import java.util.Collections;
import java.util.Iterator;
import java.util.function.Function;

/**
 * A class to examine a stream of BreakpointEvidence, and group it into Intervals.
 */
public final class EvidenceToIntervalMapper implements Function<BreakpointEvidence, Iterator<SVInterval>> {
    private static final Iterator<SVInterval> noInterval = Collections.emptyIterator();
    private final int gapSize;
    private int contig = -1;
    private int start;
    private int end;

    public EvidenceToIntervalMapper( final int gapSize ) {
        this.gapSize = gapSize;
    }

    @Override
    public Iterator<SVInterval> apply( final BreakpointEvidence evidence ) {
        Iterator<SVInterval> result = noInterval;
        if ( evidence.getContigIndex() != contig ) {
            if ( contig != -1 ) {
                result = new SingletonIterator<>(new SVInterval(contig, start, end));
            }
            contig = evidence.getContigIndex();
            start = evidence.getEventStartPosition();
            end = evidence.getContigEnd();
        } else if ( evidence.getEventStartPosition() >= end + gapSize ) {
            result = new SingletonIterator<>(new SVInterval(contig, start, end));
            start = evidence.getEventStartPosition();
            end = evidence.getContigEnd();
        } else {
            end = Math.max(end, evidence.getContigEnd());
        }
        return result;
    }
}
