package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.Iterator;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Class to find the coverage of the intervals.
 */
public final class IntervalCoverageFinder implements Iterable<Tuple2<Integer, Integer>> {
    private final int[] basesInInterval;

    public IntervalCoverageFinder( final ReadMetadata metadata,
                                   final List<SVInterval> intervals,
                                   final Iterator<GATKRead> readItr ) {
        basesInInterval = new int[intervals.size()];
        int intervalsIndex = 0;
        while ( readItr.hasNext() ) {
            final GATKRead read = readItr.next();
            final int readContigId = metadata.getContigID(read.getContig());
            final int readStart = read.getUnclippedStart();
            final int intervalsSize = intervals.size();
            while ( intervalsIndex < intervalsSize ) {
                final SVInterval interval = intervals.get(intervalsIndex);
                if ( interval.getContig() > readContigId ) break;
                if ( interval.getContig() == readContigId && interval.getEnd() > read.getStart() ) break;
                intervalsIndex += 1;
            }
            if ( intervalsIndex >= intervalsSize ) break;
            final SVInterval indexedInterval = intervals.get(intervalsIndex);
            final SVInterval readInterval = new SVInterval(readContigId, readStart, read.getUnclippedEnd()+1);
            basesInInterval[intervalsIndex] += indexedInterval.overlapLen(readInterval);
        }
    }

    @Override
    public Iterator<Tuple2<Integer, Integer>> iterator() {
        return IntStream
                .range(0, basesInInterval.length)
                .filter(idx -> basesInInterval[idx] > 0)
                .mapToObj(idx -> new Tuple2<>(idx, basesInInterval[idx]))
                .iterator();
    }
}
