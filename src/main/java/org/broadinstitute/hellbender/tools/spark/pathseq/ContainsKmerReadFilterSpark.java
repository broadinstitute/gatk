package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Spark version of ContainsKmerReadFilter
 */
public class ContainsKmerReadFilterSpark implements Function<GATKRead, Boolean> {
    private static final long serialVersionUID = 1L;
    private final String kmerSetPath;
    private final PipelineOptions options;
    private final int kmerSize, kmerCountThreshold;
    private final byte[] mask;
    private transient ContainsKmerReadFilter filter;

    public ContainsKmerReadFilterSpark(final String kmerSetPath, final int kmerSize, final byte[] mask,
                                       final int kmerCountThreshold, final PipelineOptions options) {
        this.kmerSetPath = kmerSetPath;
        this.kmerSize = kmerSize;
        this.kmerCountThreshold = kmerCountThreshold;
        this.mask = mask;
        this.options = options;
    }

    @Override
    public Boolean call(final GATKRead read) {
        if (filter == null) filter = new ContainsKmerReadFilter(kmerSetPath, kmerSize, mask, kmerCountThreshold, options);
        return filter.test(read);
    }
}
