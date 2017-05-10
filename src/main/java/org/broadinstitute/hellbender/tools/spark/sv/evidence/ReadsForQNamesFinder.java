package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Find <intervalId,list<fastqread>> pairs for interesting template names.
 */
public final class ReadsForQNamesFinder {
    private final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap;
    private final int nIntervals;
    private final int nReadsPerInterval;
    private final boolean includeMappingLocation;
    private final boolean dumpFASTQs;

    public ReadsForQNamesFinder( final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
                                 final int nIntervals, final boolean includeMappingLocation, final boolean dumpFASTQs ) {
        this.qNamesMultiMap = qNamesMultiMap;
        this.nIntervals = nIntervals;
        this.nReadsPerInterval = 2 * qNamesMultiMap.size() / nIntervals;
        this.includeMappingLocation = includeMappingLocation;
        this.dumpFASTQs = dumpFASTQs;
    }

    public Iterable<Tuple2<Integer, List<SVFastqUtils.FastqRead>>> call(final Iterator<GATKRead> readsItr ) {
        @SuppressWarnings({ "unchecked", "rawtypes" })
        final List<SVFastqUtils.FastqRead>[] intervalReads = new List[nIntervals];
        int nPopulatedIntervals = 0;
        while ( readsItr.hasNext() ) {
            final GATKRead read = readsItr.next();
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
        final List<Tuple2<Integer, List<SVFastqUtils.FastqRead>>> fastQRecords = new ArrayList<>(nPopulatedIntervals);
        if ( nPopulatedIntervals > 0 ) {
            for ( int idx = 0; idx != nIntervals; ++idx ) {
                final List<SVFastqUtils.FastqRead> readList = intervalReads[idx];
                if ( readList != null ) fastQRecords.add(new Tuple2<>(idx, readList));
            }
        }
        return fastQRecords;
    }
}
