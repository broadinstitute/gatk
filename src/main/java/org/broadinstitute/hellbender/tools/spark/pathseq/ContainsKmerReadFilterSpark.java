package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.io.ObjectInputStream;

/**
 * Wrapper for ContainsKmerReadFilter to avoid serializing the kmer filter in Spark
 */
public class ContainsKmerReadFilterSpark implements Function<GATKRead, Boolean> {
    private static final long serialVersionUID = 1L;
    private final String kmerSetPath;
    private final int kmerCountThreshold;
    private transient ContainsKmerReadFilter filter; //Load lazily to avoid its serialization

    public ContainsKmerReadFilterSpark(final String kmerSetPath, final int kmerCountThreshold) {
        this.kmerSetPath = kmerSetPath;
        this.kmerCountThreshold = kmerCountThreshold;
    }

    public void readObject(final ObjectInputStream in) throws IOException, ClassNotFoundException {
        in.defaultReadObject();
        //Load in executors upon deserialization
        filter = new ContainsKmerReadFilter(kmerSetPath, kmerCountThreshold);
    }

    @Override
    public Boolean call(final GATKRead read) {
        //Still need this check for local mode
        if (filter == null) filter = new ContainsKmerReadFilter(kmerSetPath, kmerCountThreshold);
        return filter.test(read);
    }
}
