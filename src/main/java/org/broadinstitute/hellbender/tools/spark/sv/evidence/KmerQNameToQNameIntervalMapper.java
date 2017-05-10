package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import scala.Tuple2;

import java.util.*;

/**
 * Class that maps a stream of <kmer,qname> pairs into a stream of QNameAndIntervals.
 * A multimap of kmers onto intervalIds is given to the constructor.
 * Kmers that have too many (defined by constructor param) associated qnames are discarded.
 */
public final class KmerQNameToQNameIntervalMapper {
    private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap;
    private final int maxQNamesPerKmer;
    private final int kmerMapSize;

    public KmerQNameToQNameIntervalMapper( final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap,
                                           final int maxQNamesPerKmer,
                                           final int kmerMapSize ) {
        this.kmerMultiMap = kmerMultiMap;
        this.maxQNamesPerKmer = maxQNamesPerKmer;
        this.kmerMapSize = kmerMapSize;
    }

    public Iterable<QNameAndInterval> call( final Iterator<Tuple2<SVKmer, String>> pairItr ) {
        final HopscotchMap<SVKmer, List<String>, Map.Entry<SVKmer, List<String>>> kmerQNamesMap =
                new HopscotchMap<>(kmerMapSize);
        while ( pairItr.hasNext() ) {
            final Tuple2<SVKmer, String> pair = pairItr.next();
            final SVKmer kmer = pair._1();
            Map.Entry<SVKmer, List<String>> entry = kmerQNamesMap.find(kmer);
            if ( entry == null ) {
                // new entries are created with an empty list of qnames as their value,
                // but if the list becomes too long we destroy it (by setting the value to null).
                entry = new AbstractMap.SimpleEntry<>(kmer, new ArrayList<>());
                kmerQNamesMap.add(entry);
            }
            final List<String> qNames = entry.getValue();
            // if we're still growing the list
            if ( qNames != null ) {
                // if the list becomes too long, discard it
                if ( qNames.size() >= maxQNamesPerKmer ) entry.setValue(null);
                else qNames.add(pair._2());
            }
        }

        final int qNameCount =
                kmerQNamesMap.stream().mapToInt(entry -> entry.getValue() == null ? 0 : entry.getValue().size()).sum();
        final HopscotchSet<QNameAndInterval> qNameAndIntervals = new HopscotchSet<>(qNameCount);
        for ( final Map.Entry<SVKmer, List<String>> entry : kmerQNamesMap ) {
            final List<String> qNames = entry.getValue();
            // if the list hasn't been discarded for having grown too big
            if ( qNames != null ) {
                final Iterator<KmerAndInterval> intervalItr = kmerMultiMap.findEach(entry.getKey());
                while ( intervalItr.hasNext() ) {
                    final int intervalId = intervalItr.next().getIntervalId();
                    for ( final String qName : qNames ) {
                        qNameAndIntervals.add(new QNameAndInterval(qName, intervalId));
                    }
                }
            }
        }
        return qNameAndIntervals;
    }
}
