package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.apache.spark.api.java.function.Function;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.spark.pathseq.ContainsKmerReadFilter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Set;

/**
 * Spark version of ContainsKmerReadFilter that efficiently handles kmer Set broadcasting
 */
public class ContainsKmerReadFilterSpark implements Function<GATKRead, Boolean> {
    private static final long serialVersionUID = 1L;
    final Broadcast<Set<SVKmer>> broadcastSet;
    final int kmerSize;
    transient ContainsKmerReadFilter filter;

    public ContainsKmerReadFilterSpark(final Broadcast<Set<SVKmer>> broadcastSet, final int kmerSize ) {
        this.broadcastSet = broadcastSet;
        this.kmerSize = kmerSize;
    }

    @Override
    public Boolean call( final GATKRead read ) {
        if ( filter == null ) filter = new ContainsKmerReadFilter(broadcastSet.value(), kmerSize);
        return filter.test(read);
    }
}