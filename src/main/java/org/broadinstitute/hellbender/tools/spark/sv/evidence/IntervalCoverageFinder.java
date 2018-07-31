package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Class to find the coverage of the intervals.
 */
public final class IntervalCoverageFinder implements Iterable<Tuple2<Integer, int[]>> {
    private final int[][] intervalCoverage;

    public IntervalCoverageFinder( final ReadMetadata metadata,
                                   final List<SVInterval> intervals,
                                   final Iterator<GATKRead> unfilteredReadItr,
                                   final SVReadFilter filter ) {
        intervalCoverage = new int[intervals.size()][];

        int intervalsIndex = 0;
        final int intervalsSize = intervals.size();
        final Iterator<GATKRead> readItr = filter.applyFilter(unfilteredReadItr,
                (svReadFilter, read) -> svReadFilter.isMappedToPrimaryContig(read, metadata));
        while ( readItr.hasNext() ) {
            final GATKRead read = readItr.next();
            final int readContigId = metadata.getContigID(read.getContig());
            final int readStart = read.getStart();
            while ( intervalsIndex < intervalsSize ) {
                final SVInterval interval = intervals.get(intervalsIndex);
                if ( interval.getContig() > readContigId ) break;
                if ( interval.getContig() == readContigId && interval.getEnd() > read.getStart() ) break;
                intervalsIndex += 1;
            }
            if ( intervalsIndex >= intervalsSize ) break;
            final SVInterval readInterval = new SVInterval(readContigId, readStart, read.getEnd()+1);

            int intervalsContainingReadIndex = intervalsIndex;
            while (intervalsContainingReadIndex < intervals.size()) {
                final SVInterval indexedInterval = intervals.get(intervalsContainingReadIndex);
                if ( !indexedInterval.overlaps(readInterval) ) break;

                final SVInterval overlapInterval = readInterval.intersect(indexedInterval);
                if (overlapInterval == null) throw new GATKException.ShouldNeverReachHereException("read was supposed to intersect interval but overlap is null");
                if (intervalCoverage[intervalsContainingReadIndex] == null) { 
                    intervalCoverage[intervalsContainingReadIndex] = new int[indexedInterval.getLength()];
                }
                final int[] coverageArray = intervalCoverage[intervalsContainingReadIndex];
                for (int i = overlapInterval.getStart(); i < overlapInterval.getEnd(); i++) {
                    coverageArray[i - indexedInterval.getStart()] += 1;
                }
                intervalsContainingReadIndex++;
            }
        }
    }

    static Iterator<CandidateCoverageInterval> getHighCoverageIntervalsInWindow(final int minFlankingHighCoverageValue,
                                                                                final int minPeakHighCoverageValue,
                                                                                final SVInterval interval,
                                                                                final int[] coverageArray) {
        boolean inRegion = false;
        int regionStart = -1;
        boolean foundHighRegion = false;
        final List<CandidateCoverageInterval> highCovIntervals = new ArrayList<>(2);
        for (int i = 0; i < coverageArray.length; i++) {
            if (coverageArray[i] >= minFlankingHighCoverageValue) {
                if (! inRegion) {
                    inRegion = true;
                    foundHighRegion = false;
                    regionStart = i;
                }
                if (coverageArray[i] >= minPeakHighCoverageValue) {
                    foundHighRegion = true;
                }
            } else {
                if (inRegion) {
                    if (foundHighRegion || regionStart == 0) {
                        final SVInterval highCovInterval = new SVInterval(interval.getContig(),
                                interval.getStart() + regionStart,
                                interval.getStart() + i);
                        highCovIntervals.add(new CandidateCoverageInterval(highCovInterval, foundHighRegion));
                    }
                }
                inRegion = false;
            }
        }
        if (inRegion) {
            final SVInterval highCovInterval = new SVInterval(interval.getContig(),
                    interval.getStart() + regionStart,
                    interval.getEnd());
            highCovIntervals.add(new CandidateCoverageInterval(highCovInterval, foundHighRegion));
        }
        return highCovIntervals.iterator();
    }

    @Override
    public Iterator<Tuple2<Integer, int[]>> iterator() {
        return IntStream
                .range(0, intervalCoverage.length)
                .filter(idx -> intervalCoverage[idx] != null)
                .mapToObj(idx -> new Tuple2<>(idx, intervalCoverage[idx]))
                .iterator();
    }

    @DefaultSerializer(CandidateCoverageInterval.Serializer.class)
    static class CandidateCoverageInterval {
        private static final SVInterval.Serializer intervalSerializer = new SVInterval.Serializer();
        private boolean containsMaxCoveragePeak;
        private SVInterval interval;

        public CandidateCoverageInterval(final SVInterval interval, final boolean containsMaxCoveragePeak) {
            this.containsMaxCoveragePeak = containsMaxCoveragePeak;
            this.interval = interval;
        }

        public CandidateCoverageInterval(final Kryo kryo, final Input input) {
            this.interval = intervalSerializer.read(kryo, input, SVInterval.class);
            this.containsMaxCoveragePeak = input.readBoolean();
        }

        public boolean containsMaxCoveragePeak() {
            return containsMaxCoveragePeak;
        }

        public SVInterval getInterval() {
            return interval;
        }

        protected void serialize( final Kryo kryo, final Output output ) {
            intervalSerializer.write(kryo, output, interval);
            output.writeBoolean(containsMaxCoveragePeak);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<CandidateCoverageInterval> {
            @Override
            public void write(final Kryo kryo, final Output output, final CandidateCoverageInterval candidateCoverageInterval ) {
                candidateCoverageInterval.serialize(kryo, output);
            }

            @Override
            public CandidateCoverageInterval read(final Kryo kryo, final Input input, final Class<CandidateCoverageInterval> klass ) {
                return new CandidateCoverageInterval(kryo, input);
            }
        }

    }
}
