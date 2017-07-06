package org.broadinstitute.hellbender.utils.bwa;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

/**
 * Manage a BwaMemIndex singleton.
 */
public class BwaMemIndexCache {
    private final static Map<String, BwaMemIndex> instances = new HashMap<>();
    private static String globalIndexImageFile;
    private static BwaMemIndex globalInstance;

    public static synchronized BwaMemIndex getInstance( final String indexImageFile ) {

        if (!instances.containsKey(indexImageFile)) {
            instances.put(indexImageFile, new BwaMemIndex(indexImageFile));
        }
        return instances.get(indexImageFile);
    }

    public static synchronized void closeInstances() {
        instances.values().stream().forEach(instance -> instance.close());
        instances.clear();
    }

    public static void closeAllDistributedInstances( final JavaSparkContext ctx ) {
        int nJobs = ctx.defaultParallelism();
        final List<Integer> jobList = new ArrayList<>(nJobs);
        for ( int idx = 0; idx != nJobs; ++idx ) jobList.add(idx);
        ctx.parallelize(jobList, nJobs).foreach(idx -> closeInstances());
    }
}
