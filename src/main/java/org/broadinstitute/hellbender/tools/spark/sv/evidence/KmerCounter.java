package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.KmerAndCount;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerLong;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;

/**
 * Iterates over reads, kmerizing them, and counting up just the kmers that appear in a passed-in set.
 * The counts are returned as a KmerAndCount iterator.
 */
public final class KmerCounter {
    private final int kSize;
    private final int kmersPerPartitionGuess;
    private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMap;

    public KmerCounter( final int kSize, final int kmersPerPartitionGuess,
                        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMap ) {
        this.kSize = kSize;
        this.kmerMap = kmerMap;
        this.kmersPerPartitionGuess = kmersPerPartitionGuess;
    }

    public Iterator<KmerAndCount> apply( final Iterator<GATKRead> readItr ) {
        final HopscotchMap<SVKmer, Integer, KmerAndCount> counts = new HopscotchMap<>(kmersPerPartitionGuess);
        while ( readItr.hasNext() ) {
            final GATKRead read = readItr.next();
            SVKmerizer.canonicalStream(read.getBases(), kSize, new SVKmerLong())
                    .forEach(kmer -> {
                        if ( kmerMap.contains(kmer) ) {
                            final KmerAndCount kmerAndCount = counts.find(kmer);
                            if ( kmerAndCount != null ) kmerAndCount.bumpCount();
                            else counts.add(new KmerAndCount((SVKmerLong)kmer));
                        }
                    });
        }
        return counts.iterator();
    }
}
