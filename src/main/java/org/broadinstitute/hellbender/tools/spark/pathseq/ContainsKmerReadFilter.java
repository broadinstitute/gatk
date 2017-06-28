package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmer;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerizer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Set;

/**
 * Keep reads that DO NOT contain at least one kmer from a Set of SVKmerShorts
 */
public class ContainsKmerReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;
    private Set<SVKmer> kmerLib;
    private int kSize;

    public ContainsKmerReadFilter(final Set<SVKmer> kmer_lib, int kmer_size) {
        kmerLib = kmer_lib;
        kSize = kmer_size;
    }

    @Override
    public boolean test( final GATKRead read ) {
        final SVKmerizer kmers = new SVKmerizer(read.getBases(),kSize,new SVKmerShort(kSize));
        while (kmers.hasNext()) {
            if (kmerLib.contains(kmers.next())) {return false;}
        }
        return true;
    }
}
