package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Spark version of ContainsKmerReadFilter
 */
public class ContainsKmerReadFilterSpark implements Function<GATKRead, Boolean> {
    private static final long serialVersionUID = 1L;
    private final String kmerSetPath;
    private final int kmerCountThreshold;
    private transient ContainsKmerReadFilter filter;

    public ContainsKmerReadFilterSpark(final String kmerSetPath, final int kmerCountThreshold) {
        this.kmerSetPath = kmerSetPath;
        this.kmerCountThreshold = kmerCountThreshold;
    }

    @Override
    public Boolean call(final GATKRead read) {
        if (filter == null) filter = new ContainsKmerReadFilter(kmerSetPath, kmerCountThreshold);
        return filter.test(read);
    }
}
