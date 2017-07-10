package org.broadinstitute.hellbender.utils.bwa;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

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
     * @param ctx 
     */
    public static void closeAllDistributedInstances( final JavaSparkContext ctx ) {
        int nJobs = ctx.defaultParallelism();
        final List<Integer> jobList = new ArrayList<>(nJobs);
        for ( int idx = 0; idx != nJobs; ++idx ) jobList.add(idx);
        ctx.parallelize(jobList, nJobs).foreach(idx -> closeInstances());
    }
}
