package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSetSpark;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMapSpark;
import scala.Tuple2;

import java.util.Iterator;

/**
 * Eliminates dups, and removes over-represented kmers.
 */
public final class KmerCleaner implements Iterable<KmerAndInterval> {

    private final HopscotchUniqueMultiMapSpark<SVKmer, Integer, KmerAndInterval> kmerMultiMap;

    public KmerCleaner( final Iterator<Tuple2<KmerAndInterval, Integer>> kmerCountItr,
                        final int kmersPerPartitionGuess,
                        final int minKmerCount,
                        final int maxKmerCount,
                        final int maxIntervalsPerKmer ) {
        kmerMultiMap = new HopscotchUniqueMultiMapSpark<>(kmersPerPartitionGuess);

        // remove kmers with extreme counts that won't help in building a local assembly
        while ( kmerCountItr.hasNext() ) {
            final Tuple2<KmerAndInterval, Integer> kmerCount = kmerCountItr.next();
            final int count = kmerCount._2;
            if ( count >= minKmerCount && count <= maxKmerCount ) kmerMultiMap.add(kmerCount._1);
        }

        final HopscotchSetSpark<SVKmer> uniqueKmers = new HopscotchSetSpark<>(kmerMultiMap.size());
        kmerMultiMap.stream().map(KmerAndInterval::getKey).forEach(uniqueKmers::add);
        uniqueKmers.stream()
                .filter(kmer -> SVUtils.iteratorSize(kmerMultiMap.findEach(kmer)) > maxIntervalsPerKmer)
                .forEach(kmerMultiMap::removeEach);
    }

    @Override
    public Iterator<KmerAndInterval> iterator() {
        return kmerMultiMap.iterator();
    }
}
