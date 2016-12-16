package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.tools.spark.utils.QueryableLongSet;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Keep reads that do NOT contain one or more kmers from a set of SVKmerShorts
 */
public class ContainsKmerReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;
    private static QueryableLongSet kmerLib = null;
    private final int kSize, kmerCountThreshold;
    private final byte[] mask;

    public ContainsKmerReadFilter(final String kmerLibPath, final int kSize, final byte[] mask, final int kmerCountThreshold,
                                  final PipelineOptions options) {
        this.kSize = kSize;
        this.mask = mask;
        this.kmerCountThreshold = kmerCountThreshold;
        if (kmerLib == null) {
            synchronized (ContainsKmerReadFilter.class) {
                if (kmerLib == null) {
                    kmerLib = PSKmerLibUtils.readLargeQueryableSet(kmerLibPath, options);
                }
            }
        }
    }

    @Override
    public boolean test(final GATKRead read) {
        final SVKmerizer kmers = new SVKmerizer(read.getBases(), kSize, 1, new SVKmerShort(kSize));
        int numKmersFound = 0;
        while (kmers.hasNext()) {
            if (kmerLib.contains(kmers.next().mask(mask,kSize).canonical(kSize-mask.length).getLong())) {
                if (++numKmersFound >= kmerCountThreshold) {
                    return false;
                }
            }
        }
        return true;
    }

    //Static variables can't be garbage collected until the object is unloaded
    public static void closeKmerLib() {
        kmerLib = null;
    }
}
