package org.broadinstitute.hellbender.utils.bwa;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.List;

/**
 * Manage a BwaMemIndex singleton.
 */
public class BwaMemIndexSingleton {
    private static String globalIndexImageFile;
    private static BwaMemIndex globalInstance;

    public static synchronized BwaMemIndex getInstance( final String indexImageFile ) {
        if ( globalIndexImageFile != null && !globalIndexImageFile.equals(indexImageFile) ) {
            throw new GATKException("Can't load bwa index image file "+indexImageFile+
                    " becase "+globalIndexImageFile+" is already loaded.");
        }
        if ( globalIndexImageFile == null ) {
            globalIndexImageFile = indexImageFile;
            globalInstance = new BwaMemIndex(indexImageFile);
        }
        return globalInstance;
    }

    public static synchronized void closeInstance() {
        if ( globalInstance != null ) {
            globalInstance.close();
            globalIndexImageFile = null;
            globalInstance = null;
        }
    }

    public static void closeAllDistributedInstances( final JavaSparkContext ctx ) {
        int nJobs = ctx.defaultParallelism();
        final List<Integer> jobList = new ArrayList<>(nJobs);
        for ( int idx = 0; idx != nJobs; ++idx ) jobList.add(idx);
        ctx.parallelize(jobList, nJobs).foreach(idx -> closeInstance());
    }
}
