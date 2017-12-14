package org.broadinstitute.hellbender.utils.bwa;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

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
     * Closes an index instance in the cache given its index file name.
     * <p>
     *     Notice that you need to pass in exactly the same file name that was used when invoking {@link #getInstance}.
     * </p>
     * <p>
     *     An attempt to close a missing instance, won't have any effect.
     * </p>
     *
     * @param indexImageFile the index file name of the instance to close.
     */
    public static synchronized void closeInstance(final String indexImageFile) {
        Utils.nonNull(indexImageFile, "the input image file cannot be null");
        if (instances.containsKey(indexImageFile)) {
            instances.get(indexImageFile).close();
            instances.remove(indexImageFile);
        }
    }

    /**
     * Closes an index instance.
     *<p>
     *     An attempt to close a instance that is not present in the cache, won't have any effect.
     *     Thus if the input instance is not part of the cache an is not closed, will remind unclosed.
     * </p>
     * @param instance the instance ot close.
     */
    public static synchronized void closeInstance(final BwaMemIndex instance) {
        Utils.nonNull(instance, "the input index cannot be null");
        if (instances.values().contains(instance)) {
            instance.close();
            instances.values().remove(instance);
        }
    }

    /**
     * Closes all instances in the VM.
     */
    public static synchronized void closeInstances() {
        final Iterator<BwaMemIndex> it = instances.values().iterator();
        while (it.hasNext()) {
            it.next().close();
            it.remove();
        }
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
