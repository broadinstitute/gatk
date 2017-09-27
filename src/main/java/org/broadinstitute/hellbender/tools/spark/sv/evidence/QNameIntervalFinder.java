package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerLong;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

/**
 * Iterates over reads, kmerizing them, and checking the kmers against a set of KmerAndIntervals
 * to figure out which intervals (if any) a read belongs in.
 * Results are returned as a QNameAndInterval iterator.
 */
public final class QNameIntervalFinder implements Function<GATKRead,Iterator<QNameAndInterval>> {
    private final int kSize;
    private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMap;

    public QNameIntervalFinder( final int kSize, final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMap ) {
        this.kSize = kSize;
        this.kmerMap = kmerMap;
    }

    @Override
    public Iterator<QNameAndInterval> apply( final GATKRead read ) {
        final List<Integer> intervals = new ArrayList<>();
        SVKmerizer.canonicalStream(read.getBases(), kSize, new SVKmerLong())
                .forEach(kmer -> {
                    final Iterator<KmerAndInterval> kmerAndIntervalIterator = kmerMap.findEach(kmer);
                    while ( kmerAndIntervalIterator.hasNext() ) {
                        final Integer intervalId = kmerAndIntervalIterator.next().getValue();
                        if ( !intervals.contains(intervalId) ) {
                            intervals.add(intervalId);
                        }
                    }
                });
        final String qName = read.getName();
        return intervals.stream().map(intervalId -> new QNameAndInterval(qName, intervalId)).iterator();
    }
}
