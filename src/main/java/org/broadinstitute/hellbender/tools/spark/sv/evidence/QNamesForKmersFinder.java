package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerLong;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;

/**
 * Class that acts as a mapper from a stream of reads to a stream of <kmer,qname> pairs for a set of interesting kmers.
 * A multimap of interesting kmers is given to the constructor (by broadcast).
 */
public final class QNamesForKmersFinder implements Function<GATKRead, Iterator<Tuple2<SVKmer, String>>> {
    private final int kSize;
    private final SVReadFilter filter;
    private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap;

    public QNamesForKmersFinder( final int kSize,
                                 final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap,
                                 final SVReadFilter filter ) {
        this.kSize = kSize;
        this.kmerMultiMap = kmerMultiMap;
        this.filter = filter;
    }

    @Override
    public Iterator<Tuple2<SVKmer, String>> apply( final GATKRead read ) {
        if ( !filter.notJunk(read) || !filter.isPrimaryLine(read) ) return Collections.emptyIterator();

        final List<Tuple2<SVKmer, String>> results = new ArrayList<>();
        SVKmerizer.canonicalStream(read.getBases(), kSize, new SVKmerLong())
                .forEach(kmer -> {
                    if ( kmerMultiMap.contains(kmer) ) {
                        results.add(new Tuple2<>(kmer, read.getName()));
                    }
                });
        return results.iterator();
    }
}
