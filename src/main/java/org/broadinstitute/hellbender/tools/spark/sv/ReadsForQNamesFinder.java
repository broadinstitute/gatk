package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Find <intervalId,list<fastqread>> pairs for interesting template names.
 */
public final class ReadsForQNamesFinder implements Iterable<Tuple2<Integer, List<SVFastqUtils.FastqRead>>> {
    private final List<Tuple2<Integer, List<SVFastqUtils.FastqRead>>> fastQRecords;

    public ReadsForQNamesFinder( final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
                                 final int nIntervals, final boolean includeMappingLocation, final boolean dumpFASTQs,
                                 final Iterator<GATKRead> unfilteredReadItr, final SVReadFilter filter ) {
        final int nReadsPerInterval = 2 * qNamesMultiMap.size() / nIntervals;
        @SuppressWarnings({ "unchecked", "rawtypes" })
        final List<SVFastqUtils.FastqRead>[] intervalReads = new List[nIntervals];
        int nPopulatedIntervals = 0;
        final Iterator<GATKRead> readItr = filter.applyFilter(unfilteredReadItr, SVReadFilter::isPrimaryLine);
        while ( readItr.hasNext() ) {
            final GATKRead read = readItr.next();
            final Iterator<QNameAndInterval> namesItr = qNamesMultiMap.findEach(read.getName());
            SVFastqUtils.FastqRead FastqRead = null;
            while ( namesItr.hasNext() ) {
                final int intervalId = namesItr.next().getIntervalId();
                if ( intervalReads[intervalId] == null ) {
                    intervalReads[intervalId] = new ArrayList<>(nReadsPerInterval);
                    nPopulatedIntervals += 1;
                }
                if ( FastqRead == null ) {
                    final String readName =
                            dumpFASTQs ? SVFastqUtils.readToFastqSeqId(read, includeMappingLocation) : null;
                    FastqRead = new SVFastqUtils.FastqRead(readName, read.getBases(), read.getBaseQualities());
                }
                intervalReads[intervalId].add(FastqRead);
            }
        }
        fastQRecords = new ArrayList<>(nPopulatedIntervals);
        if ( nPopulatedIntervals > 0 ) {
            for ( int idx = 0; idx != nIntervals; ++idx ) {
                final List<SVFastqUtils.FastqRead> readList = intervalReads[idx];
                if ( readList != null ) fastQRecords.add(new Tuple2<>(idx, readList));
            }
        }
    }

    @Override
    public Iterator<Tuple2<Integer, List<SVFastqUtils.FastqRead>>> iterator() { return fastQRecords.iterator(); }
}
