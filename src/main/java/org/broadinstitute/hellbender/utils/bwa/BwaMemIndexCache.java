package org.broadinstitute.hellbender.utils.bwa;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Manage a global collection of {@link BwaMemIndex} instances.
 */
public class BwaMemIndexCache {

    private final static Map<String, BwaMemIndex> instances = new HashMap<>();

    /**
     * Returns a {@link BwaMemIndex} instance that corresponds to  given index image file.
     * @param indexImageFile the target image file.
     * @return never {@code null}.
     */
    public static synchronized BwaMemIndex getInstance( final String indexImageFile ) {
        Utils.nonNull(indexImageFile, "the index image file name provided cannot be null");
        if (!instances.containsKey(indexImageFile)) {
            instances.put(indexImageFile, new BwaMemIndex(indexImageFile));
        }
        return instances.get(indexImageFile);
    }

    /**
     * Closes all instances in the VM.
     */
    public static synchronized void closeInstances() {
        instances.values().forEach(BwaMemIndex::close);
        instances.clear();
    }

    /**
     * Closes all instances in all the VMs involved in the spark context provided.
     * @param ctx the spark context.
     */
    public static void closeAllDistributedInstances( final JavaSparkContext ctx ) {
        Utils.nonNull(ctx, "the context provided cannot be null");
        int nJobs = ctx.defaultParallelism();
        final List<Integer> jobList = new ArrayList<>(nJobs);
        for ( int idx = 0; idx != nJobs; ++idx ) jobList.add(idx);
        ctx.parallelize(jobList, nJobs).foreach(idx -> closeInstances());
    }
}
