package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVDUSTFilteredKmerizer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerLong;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;
import java.util.function.Function;

/**
 * Class that acts as a mapper from a stream of reads to a stream of KmerAndIntervals.
 * The template names of reads to kmerize, along with a set of kmers to ignore are passed in (by broadcast).
 */
public final class QNameKmerizer implements Function<GATKRead, Iterator<Tuple2<KmerAndInterval, Integer>>> {
    private final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNameAndIntervalMultiMap;
    private final Set<SVKmer> kmersToIgnore;
    private final int kSize;
    private final int maxDUSTScore;
    private final ArrayList<Tuple2<KmerAndInterval, Integer>> tupleList = new ArrayList<>();

    public QNameKmerizer( final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNameAndIntervalMultiMap,
                          final Set<SVKmer> kmersToIgnore, final int kSize, final int maxDUSTScore ) {
        this.qNameAndIntervalMultiMap = qNameAndIntervalMultiMap;
        this.kmersToIgnore = kmersToIgnore;
        this.kSize = kSize;
        this.maxDUSTScore = maxDUSTScore;
    }

    @Override
    public Iterator<Tuple2<KmerAndInterval, Integer>> apply( final GATKRead read ) {
        final String qName = read.getName();
        final Iterator<QNameAndInterval> names = qNameAndIntervalMultiMap.findEach(qName);
        tupleList.clear();
        while ( names.hasNext() ) {
            final int intervalId = names.next().getIntervalId();
            SVDUSTFilteredKmerizer.stream(read.getBases(), kSize, maxDUSTScore, new SVKmerLong())
                    .map(kmer -> kmer.canonical(kSize))
                    .filter(kmer -> !kmersToIgnore.contains(kmer))
                    .map(kmer -> new KmerAndInterval(kmer, intervalId))
                    .forEach(kmerCountAndInterval -> tupleList.add(new Tuple2<>(kmerCountAndInterval, 1)));
        }
        return tupleList.iterator();
    }
}
